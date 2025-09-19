"""
Paper collection service for PubMed data retrieval.
"""

from Bio import Entrez
import time
import logging
from datetime import datetime
from typing import List, Dict, Optional

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

        # Set Entrez email
        Entrez.email = self.email

        logger.info(f"Initialized PaperCollectionService with email: {email}")

    def search_papers(self, query: str, max_results: int = 1000) -> List[str]:
        """Search PubMed for papers matching the query.

        Args:
            query: PubMed search query string
            max_results: Maximum number of papers to retrieve

        Returns:
            List of PMIDs as strings

        Raises:
            APIError: If PubMed API call fails
            ValueError: If query is empty or max_results is invalid
        """
        if not query.strip():
            raise ValueError("Query cannot be empty")
        if max_results <= 0 or max_results > 10000:
            raise ValueError("max_results must be between 1 and 10000")

        self._rate_limit()

        try:
            logger.info(f"Searching PubMed for: '{query}' (max: {max_results})")

            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="pub+date",  # Sort by publication date (newest first)
                retmode="xml",
            )

            records = Entrez.read(handle)
            handle.close()

            pmids = records.get("IdList", [])
            logger.info(f"Found {len(pmids)} papers for query: '{query}'")

            return pmids

        except Exception as e:
            logger.error(f"Error searching PubMed for query '{query}': {e}")
            raise APIError(f"PubMed search failed: {e}")

    def get_paper_details(self, pmids: List[str]) -> List[Paper]:
        """Retrieve detailed information for a list of PMIDs.

        Args:
            pmids: List of PubMed IDs

        Returns:
            List of Paper objects with detailed metadata

        Raises:
            APIError: If PubMed API call fails
        """
        if not pmids:
            return []

        papers = []
        batch_size = 100  # Process in batches to avoid API limits

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i : i + batch_size]
            batch_papers = self._fetch_paper_batch(batch)
            papers.extend(batch_papers)

            # Rate limiting between batches
            if i + batch_size < len(pmids):
                time.sleep(1.0 / self.rate_limit)

        logger.info(f"Retrieved details for {len(papers)} papers")
        return papers

    def fetch_abstracts(self, pmids: List[str]) -> Dict[str, str]:
        """Fetch abstracts for a list of PMIDs.

        Args:
            pmids: List of PubMed IDs

        Returns:
            Dictionary mapping PMID to abstract text
        """
        if not pmids:
            return {}

        abstracts = {}
        batch_size = 50  # Smaller batch size for abstracts

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i : i + batch_size]
            batch_abstracts = self._fetch_abstracts_batch(batch)
            abstracts.update(batch_abstracts)

            # Rate limiting
            if i + batch_size < len(pmids):
                time.sleep(1.0 / self.rate_limit)

        logger.info(f"Retrieved abstracts for {len(abstracts)} papers")
        return abstracts

    def _fetch_paper_batch(self, pmids: List[str]) -> List[Paper]:
        """Fetch a batch of paper details from PubMed.

        Args:
            pmids: List of PubMed IDs (max 100)

        Returns:
            List of Paper objects
        """
        self._rate_limit()

        try:
            handle = Entrez.efetch(
                db="pubmed", id=",".join(pmids), rettype="medline", retmode="xml"
            )

            records = Entrez.read(handle)
            handle.close()

            papers = []
            for record in records["PubmedArticle"]:
                try:
                    paper = self._parse_paper_record(record)
                    if paper:
                        papers.append(paper)
                except Exception as e:
                    pmid = record.get("MedlineCitation", {}).get("PMID", "Unknown")
                    logger.warning(f"Error parsing paper {pmid}: {e}")
                    continue

            return papers

        except Exception as e:
            logger.error(f"Error fetching paper batch {pmids[:3]}...: {e}")
            raise APIError(f"Failed to fetch paper details: {e}")

    def _fetch_abstracts_batch(self, pmids: List[str]) -> Dict[str, str]:
        """Fetch abstracts for a batch of PMIDs.

        Args:
            pmids: List of PubMed IDs

        Returns:
            Dictionary mapping PMID to abstract
        """
        self._rate_limit()

        try:
            handle = Entrez.efetch(
                db="pubmed", id=",".join(pmids), rettype="abstract", retmode="xml"
            )

            records = Entrez.read(handle)
            handle.close()

            abstracts = {}
            for record in records["PubmedArticle"]:
                try:
                    pmid = str(record["MedlineCitation"]["PMID"])
                    abstract = self._extract_abstract(record)
                    if abstract:
                        abstracts[pmid] = abstract
                except Exception as e:
                    logger.warning(f"Error extracting abstract: {e}")
                    continue

            return abstracts

        except Exception as e:
            logger.error(f"Error fetching abstracts for batch: {e}")
            return {}

    def _parse_paper_record(self, record: Dict) -> Optional[Paper]:
        """Parse a PubMed record into a Paper object.

        Args:
            record: PubMed XML record as dictionary

        Returns:
            Paper object or None if parsing fails
        """
        try:
            citation = record["MedlineCitation"]
            article = citation["Article"]

            # Extract PMID
            pmid = str(citation["PMID"])

            # Extract title
            title = article.get("ArticleTitle", "").strip()
            if not title:
                logger.warning(f"No title found for PMID {pmid}")
                return None

            # Extract authors
            authors = self._extract_authors(article)
            if not authors:
                logger.warning(f"No authors found for PMID {pmid}")
                authors = ["Unknown"]

            # Extract journal
            journal = self._extract_journal(article)
            if not journal:
                logger.warning(f"No journal found for PMID {pmid}")
                journal = "Unknown Journal"

            # Extract publication date
            pub_date = self._extract_publication_date(article)

            # Extract DOI
            doi = self._extract_doi(article)

            # Extract abstract
            abstract = self._extract_abstract(record)

            return Paper(
                pmid=pmid,
                title=title,
                authors=authors,
                journal=journal,
                publication_date=pub_date,
                abstract=abstract,
                doi=doi,
            )

        except Exception as e:
            logger.error(f"Error parsing paper record: {e}")
            return None

    def _extract_authors(self, article: Dict) -> List[str]:
        """Extract author names from article record."""
        authors = []

        author_list = article.get("AuthorList", [])
        for author in author_list:
            if "LastName" in author and "ForeName" in author:
                name = f"{author['ForeName']} {author['LastName']}"
                authors.append(name)
            elif "CollectiveName" in author:
                authors.append(author["CollectiveName"])

        return authors

    def _extract_journal(self, article: Dict) -> str:
        """Extract journal name from article record."""
        journal_info = article.get("Journal", {})

        # Try full journal title first
        title = journal_info.get("Title", "")
        if title:
            return title

        # Fall back to ISO abbreviation
        iso_abbrev = journal_info.get("ISOAbbreviation", "")
        if iso_abbrev:
            return iso_abbrev

        return "Unknown Journal"

    def _extract_publication_date(self, article: Dict) -> datetime:
        """Extract publication date from article record."""
        try:
            # Try ArticleDate first
            if "ArticleDate" in article and article["ArticleDate"]:
                date_info = article["ArticleDate"][0]
                year = int(date_info.get("Year", datetime.now().year))
                month = int(date_info.get("Month", 1))
                day = int(date_info.get("Day", 1))
                return datetime(year, month, day)

            # Fall back to Journal publication date
            journal = article.get("Journal", {})
            journal_issue = journal.get("JournalIssue", {})
            pub_date = journal_issue.get("PubDate", {})

            year = pub_date.get("Year")
            if year:
                month = pub_date.get("Month", "1")
                day = pub_date.get("Day", "1")

                # Handle month names
                month_map = {
                    "Jan": 1,
                    "Feb": 2,
                    "Mar": 3,
                    "Apr": 4,
                    "May": 5,
                    "Jun": 6,
                    "Jul": 7,
                    "Aug": 8,
                    "Sep": 9,
                    "Oct": 10,
                    "Nov": 11,
                    "Dec": 12,
                }

                if isinstance(month, str) and month in month_map:
                    month = month_map[month]
                else:
                    month = int(month) if str(month).isdigit() else 1

                day = int(day) if str(day).isdigit() else 1
                year = int(year)

                return datetime(year, month, day)

        except Exception as e:
            logger.warning(f"Error parsing publication date: {e}")

        # Default to current date if parsing fails
        return datetime.now()

    def _extract_doi(self, article: Dict) -> Optional[str]:
        """Extract DOI from article record."""
        try:
            elocation_id = article.get("ELocationID", [])
            for location in elocation_id:
                if location.attributes.get("EIdType") == "doi":
                    return str(location)
        except Exception:
            pass

        return None

    def _extract_abstract(self, record: Dict) -> Optional[str]:
        """Extract abstract text from record."""
        try:
            citation = record["MedlineCitation"]
            article = citation["Article"]
            abstract_info = article.get("Abstract", {})

            if "AbstractText" in abstract_info:
                abstract_parts = abstract_info["AbstractText"]
                if isinstance(abstract_parts, list):
                    # Join multiple abstract parts
                    return " ".join(str(part) for part in abstract_parts)
                else:
                    return str(abstract_parts)

        except Exception as e:
            logger.warning(f"Error extracting abstract: {e}")

        return None

    def _rate_limit(self) -> None:
        """Implement rate limiting for API calls."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        min_interval = 1.0 / self.rate_limit

        if time_since_last < min_interval:
            sleep_time = min_interval - time_since_last
            time.sleep(sleep_time)

        self.last_request_time = time.time()
