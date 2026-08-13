"""
Paper collection service for PubMed data retrieval.
"""

import time
import logging
from datetime import datetime
from typing import List, Dict, Optional

from Bio import Entrez

from ..models import Paper
from ..utils.error_handler import APIError

logger = logging.getLogger(__name__)


class PaperCollectionService:
    """Service for collecting paper data from PubMed."""

    def __init__(
        self, email: str = "pubmed.miner@example.com", rate_limit: float = 3.0
    ):
        """Initialize the paper collection service.

        Args:
            email: Email address for Entrez API (required by NCBI)
            rate_limit: Maximum requests per second to PubMed
        """
        self.email = email
        self.rate_limit = rate_limit
        self.last_request_time = 0.0

        self._total_searches = 0
        self._total_papers_fetched = 0
        self._total_api_calls = 0
        self._average_response_time = 0.0
        self._error_count = 0
        self._rate_limit_hits = 0

        Entrez.email = self.email

        logger.info(f"Initialized PaperCollectionService with email: {email}")

    def search_papers(
        self,
        query: str,
        max_results: int = 1000,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None,
    ) -> List[str]:
        """Search PubMed for papers matching the query.

        Args:
            query: PubMed search query string
            max_results: Maximum number of papers to retrieve

        Returns:
            List of PMIDs as strings

        Raises:
            ValueError: If query is empty or max_results is invalid
            APIError: If PubMed API call fails
        """
        if not query.strip():
            raise ValueError("Query cannot be empty")

        if max_results <= 0:
            raise ValueError("Max results must be positive")

        self._rate_limit()

        self._total_searches += 1

        max_retries = 3
        for attempt in range(max_retries):
            try:
                logger.info(f"Searching PubMed for: '{query}' (max: {max_results})")

                term_parts = [query]
                if date_from:
                    term_parts.append(f"({date_from}:{date_to}[Date - Publication])")
                full_term = " AND ".join(term_parts)

                handle = Entrez.esearch(
                    db="pubmed",
                    term=full_term,
                    retmax=max_results,
                    sort="pub+date",
                    retmode="xml",
                )

                records = Entrez.read(handle)
                handle.close()

                pmids = records.get("IdList", [])
                logger.info(f"Found {len(pmids)} papers for query: '{query}'")

                return pmids
            except Exception as e:
                if attempt == max_retries - 1:
                    logger.error(f"Error searching PubMed for query '{query}': {e}")
                    raise APIError(f"PubMed search failed: {e}")
                logger.warning(
                    f"Retry {attempt + 1}/{max_retries} for query '{query}': {e}"
                )
                continue

        # Should not reach here
        raise APIError(f"PubMed search failed: max retries exceeded")

    def get_paper_details(
        self, pmids: List[str], topic: Optional[str] = None
    ) -> List[Paper]:
        """Retrieve detailed information for a list of PMIDs.

        Args:
            pmids: List of PubMed IDs
            topic: Optional topic name to associate with papers

        Returns:
            List of Paper objects with detailed metadata

        Raises:
            ValueError: If PMID list contains invalid values
            APIError: If PubMed API call fails
        """
        if not pmids:
            return []

        for pmid in pmids:
            if not pmid or not isinstance(pmid, str):
                raise ValueError("All PMIDs must be non-empty strings")

        batches: List[List[str]] = []
        for i in range(0, len(pmids), 100):
            batches.append(pmids[i:i + 100])

        all_papers: List[Paper] = []
        for batch in batches:
            batch_papers = self._fetch_paper_batch(batch, topic=topic)
            all_papers.extend(batch_papers)

        return all_papers

    def _fetch_paper_batch(
        self, pmids: List[str], topic: Optional[str] = None
    ) -> List[Paper]:
        """Fetch a batch of paper details from PubMed.

        Args:
            pmids: List of PubMed IDs (max 100)
            topic: Optional topic name to associate with papers

        Returns:
            List of Paper objects
        """
        self._rate_limit()

        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(pmids),
                rettype="medline",
                retmode="xml",
            )

            records = Entrez.read(handle)
            handle.close()

            papers = []
            for record in records["PubmedArticle"]:
                try:
                    paper = self._parse_paper_record(record, topic=topic)
                    if paper:
                        papers.append(paper)
                except Exception as e:
                    pmid = record.get("MedlineCitation", {}).get("PMID", "Unknown")
                    logger.warning(f"Error parsing paper {pmid}: {e}")
                    continue

            self._total_papers_fetched += len(papers)
            return papers

        except Exception as e:
            logger.error(f"Error fetching paper batch {pmids[:3]}...: {e}")
            raise APIError(f"Failed to fetch paper details: {e}")

    def _parse_paper_record(
        self, record: Dict, topic: Optional[str] = None
    ) -> Optional[Paper]:
        """Parse a PubMed record into a Paper object.

        Args:
            record: PubMed XML record as dictionary
            topic: Optional topic name to associate with the paper

        Returns:
            Paper object or None if parsing fails
        """
        try:
            citation = record["MedlineCitation"]
            article = citation["Article"]

            pmid = str(citation["PMID"])

            title = article.get("ArticleTitle", "No Title")

            authors = self._parse_authors(
                article.get("AuthorList", [])
            )

            journal = article.get("Journal", {}).get("Title", "Unknown Journal")

            pub_date = self._parse_date(article.get("ArticleDate", {}))

            abstract = self._extract_abstract(record)

            doi = self._extract_doi(record)

            return Paper(
                pmid=pmid,
                title=title,
                authors=authors,
                journal=journal,
                publication_date=pub_date,
                abstract=abstract,
                doi=doi,
                topic=topic,
            )

        except Exception as e:
            logger.warning(f"Error parsing paper record: {e}")
            return None

    def _parse_authors(self, author_list) -> List[str]:
        """Parse author list from XML.

        Args:
            author_list: Author list from PubMed record (dict list, XML string, or XML element)

        Returns:
            List of author full names
        """
        # If it's an XML string, parse it first
        if isinstance(author_list, str):
            from xml.etree import ElementTree as ET
            try:
                root = ET.fromstring(author_list)
                author_list = []
                for author in root.findall("Author"):
                    author_dict = {}
                    last = author.findtext("LastName")
                    fore = author.findtext("ForeName")
                    if last:
                        author_dict["LastName"] = last
                    if fore:
                        author_dict["ForeName"] = fore
                    if author_dict:
                        author_list.append(author_dict)
            except Exception:
                return []

        authors = []
        for author in author_list:
            last_name = author.get("LastName", "")
            fore_name = author.get("ForeName", "")
            if last_name and fore_name:
                authors.append(f"{fore_name} {last_name}")
            elif last_name:
                authors.append(last_name)
        return authors

    def _extract_abstract(self, record: Dict) -> Optional[str]:
        """Extract abstract from PubMed record.

        Args:
            record: PubMed XML record as dictionary

        Returns:
            Abstract text or None
        """
        try:
            citation = record["MedlineCitation"]
            article = citation["Article"]
            abstract_text = article.get("Abstract", {}).get("AbstractText", "")
            if abstract_text:
                return str(abstract_text)
        except (KeyError, AttributeError):
            pass
        return None

    def _extract_doi(self, record: Dict) -> Optional[str]:
        """Extract DOI from PubMed record.

        Args:
            record: PubMed XML record as dictionary

        Returns:
            DOI string or None
        """
        try:
            if "PubmedData" in record:
                article_ids = record["PubmedData"].get("ArticleIdList", [])
                for article_id in article_ids:
                    if isinstance(article_id, dict):
                        if article_id.get("IdType") == "doi":
                            return article_id.get("value")
                    elif hasattr(article_id, "attributes"):
                        if article_id.attributes.get("IdType") == "doi":
                            return str(article_id)
        except (KeyError, AttributeError):
            pass
        return None

    def _rate_limit(self) -> None:
        """Apply rate limiting between API calls."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        min_interval = 1.0 / self.rate_limit

        if time_since_last < min_interval:
            sleep_time = min_interval - time_since_last
            logger.debug(f"Rate limiting: sleeping {sleep_time:.2f}s")
            time.sleep(sleep_time)

        self.last_request_time = time.time()

    def _validate_pmid(self, pmid: Optional[str]) -> bool:
        """Validate PMID format.

        Args:
            pmid: PubMed ID to validate

        Returns:
            True if valid, False otherwise
        """
        if not pmid:
            return False
        if not isinstance(pmid, str):
            return False
        if not pmid.isdigit():
            return False
        if len(pmid) > 8:
            return False
        return True

    def reset_statistics(self) -> None:
        """Reset all statistics counters to zero."""
        self._total_searches = 0
        self._total_papers_fetched = 0
        self._total_api_calls = 0
        self._average_response_time = 0.0
        self._error_count = 0
        self._rate_limit_hits = 0
        self.last_request_time = 0.0

    def get_statistics(self) -> Dict[str, Any]:
        """Get current statistics for the service.

        Returns:
            Dictionary with statistic names and values.
        """
        return {
            "total_searches": self._total_searches,
            "total_papers_fetched": self._total_papers_fetched,
            "total_api_calls": self._total_api_calls,
            "average_response_time": self._average_response_time,
            "error_count": self._error_count,
            "rate_limit_hits": self._rate_limit_hits,
        }

    def _parse_authors_from_xml(self, author_data) -> str:
        """Parse author information from XML data.

        Args:
            author_data: XML element or data structure containing author information.

        Returns:
            String containing formatted author names.
        """
        if isinstance(author_data, str):
            return author_data

        if hasattr(author_data, "findall") and callable(getattr(author_data, "findall", None)):
            authors = []
            for author in author_data.findall("Author"):
                last_name = author.findtext("LastName", "")
                fore_name = author.findtext("ForeName", "")
                if last_name and fore_name:
                    authors.append(f"{fore_name} {last_name}")
                elif last_name:
                    authors.append(last_name)
            return ", ".join(authors)

        if isinstance(author_data, dict):
            last_name = author_data.get("LastName", "")
            fore_name = author_data.get("ForeName", "")
            if last_name and fore_name:
                return f"{fore_name} {last_name}"
            return last_name or fore_name or ""

        return str(author_data)

    def _parse_date(self, date_dict) -> datetime:
        """Parse publication date from XML.

        Args:
            date_dict: Date dictionary from PubMed record (dict, XML string, or XML element)

        Returns:
            datetime object (defaults to epoch if parsing fails)
        """
        # Handle XML string
        if isinstance(date_dict, str):
            try:
                from xml.etree import ElementTree as ET
                root = ET.fromstring(date_dict)
                year_text = root.findtext("Year", "2000")
                month_text = root.findtext("Month", "1")
                day_text = root.findtext("Day", "1")
                year = int(year_text) if year_text else 2000
                month = int(month_text) if month_text else 1
                day = int(day_text) if day_text else 1
                return datetime(year, month, day)
            except Exception:
                pass
            # Fallback for ISO format string
            try:
                return datetime.fromisoformat(date_dict)
            except (ValueError, TypeError):
                return datetime(2000, 1, 1)

        # Handle XML element
        if hasattr(date_dict, "findtext") and callable(getattr(date_dict, "findtext", None)):
            try:
                year_text = date_dict.findtext("Year", "2000")
                month_text = date_dict.findtext("Month", "1")
                day_text = date_dict.findtext("Day", "1")
                year = int(year_text) if year_text else 2000
                month = int(month_text) if month_text else 1
                day = int(day_text) if day_text else 1
                return datetime(year, month, day)
            except (ValueError, TypeError):
                return datetime(2000, 1, 1)

        # Handle datetime
        if isinstance(date_dict, datetime):
            return date_dict

        # Handle dict
        try:
            year = int(date_dict.get("Year", "2000"))
            month = int(date_dict.get("Month", "1"))
            day = int(date_dict.get("Day", "1"))
            return datetime(year, month, day)
        except (ValueError, KeyError, TypeError, AttributeError):
            return datetime(2000, 1, 1)

    def _parse_date_from_xml(self, date_dict) -> datetime:
        """Parse publication date from XML string or element.

        Args:
            date_dict: XML string or element with date info

        Returns:
            datetime object
        """
        return self._parse_date(date_dict)

    def _parse_date_from_element(self, date_element) -> datetime:
        """Parse date from XML element."""
        return self._parse_date(date_element)

    def _clean_text(self, text: str) -> str:
        """Clean text by unescaping HTML entities and collapsing whitespace.
        
        Args:
            text: Text to clean
            
        Returns:
            Cleaned text
        """
        import html as _html
        cleaned = _html.unescape(text)
        cleaned = ' '.join(cleaned.split())
        return cleaned
