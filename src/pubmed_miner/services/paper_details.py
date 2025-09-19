"""
Enhanced paper details collection service.
"""

import re
import logging
from typing import Dict, List, Optional
from Bio import Entrez

from ..models import Paper
from ..utils.error_handler import retry_api_calls

logger = logging.getLogger(__name__)


class PaperDetailsService:
    """Enhanced service for collecting detailed paper information."""

    def __init__(self, email: str = "pubmed.miner@example.com"):
        """Initialize the paper details service.

        Args:
            email: Email address for Entrez API
        """
        self.email = email
        Entrez.email = self.email

    @retry_api_calls(max_attempts=3, delay=1.0)
    def enrich_paper_data(self, papers: List[Paper]) -> List[Paper]:
        """Enrich paper data with additional metadata.

        Args:
            papers: List of Paper objects to enrich

        Returns:
            List of enriched Paper objects
        """
        if not papers:
            return papers

        enriched_papers = []

        for paper in papers:
            try:
                enriched_paper = self._enrich_single_paper(paper)
                enriched_papers.append(enriched_paper)
            except Exception as e:
                logger.warning(f"Failed to enrich paper {paper.pmid}: {e}")
                # Keep original paper if enrichment fails
                enriched_papers.append(paper)

        logger.info(f"Enriched {len(enriched_papers)} papers")
        return enriched_papers

    def _enrich_single_paper(self, paper: Paper) -> Paper:
        """Enrich a single paper with additional data.

        Args:
            paper: Paper object to enrich

        Returns:
            Enriched Paper object
        """
        # Get additional metadata from PubMed
        additional_data = self._fetch_additional_metadata(paper.pmid)

        # Create enriched paper with updated data
        enriched_data = {
            "pmid": paper.pmid,
            "title": self._clean_title(paper.title),
            "authors": self._normalize_authors(paper.authors),
            "journal": self._normalize_journal_name(paper.journal),
            "publication_date": paper.publication_date,
            "abstract": self._clean_abstract(
                paper.abstract or additional_data.get("abstract")
            ),
            "doi": paper.doi or additional_data.get("doi"),
        }

        return Paper(**enriched_data)

    @retry_api_calls(max_attempts=2, delay=0.5)
    def _fetch_additional_metadata(self, pmid: str) -> Dict[str, Optional[str]]:
        """Fetch additional metadata for a paper.

        Args:
            pmid: PubMed ID

        Returns:
            Dictionary with additional metadata
        """
        try:
            handle = Entrez.efetch(
                db="pubmed", id=pmid, rettype="medline", retmode="xml"
            )

            records = Entrez.read(handle)
            handle.close()

            if not records["PubmedArticle"]:
                return {}

            record = records["PubmedArticle"][0]

            return {
                "abstract": self._extract_full_abstract(record),
                "doi": self._extract_doi_comprehensive(record),
                "keywords": self._extract_keywords(record),
                "mesh_terms": self._extract_mesh_terms(record),
            }

        except Exception as e:
            logger.warning(f"Failed to fetch additional metadata for {pmid}: {e}")
            return {}

    def _extract_full_abstract(self, record: Dict) -> Optional[str]:
        """Extract complete abstract including structured abstracts.

        Args:
            record: PubMed XML record

        Returns:
            Complete abstract text or None
        """
        try:
            article = record["MedlineCitation"]["Article"]
            abstract_info = article.get("Abstract", {})

            if "AbstractText" not in abstract_info:
                return None

            abstract_parts = abstract_info["AbstractText"]

            if isinstance(abstract_parts, list):
                # Handle structured abstracts
                structured_parts = []
                for part in abstract_parts:
                    if hasattr(part, "attributes") and "Label" in part.attributes:
                        label = part.attributes["Label"]
                        text = str(part)
                        structured_parts.append(f"{label}: {text}")
                    else:
                        structured_parts.append(str(part))
                return " ".join(structured_parts)
            else:
                return str(abstract_parts)

        except Exception as e:
            logger.warning(f"Error extracting full abstract: {e}")
            return None

    def _extract_doi_comprehensive(self, record: Dict) -> Optional[str]:
        """Extract DOI using multiple methods.

        Args:
            record: PubMed XML record

        Returns:
            DOI string or None
        """
        try:
            article = record["MedlineCitation"]["Article"]

            # Method 1: ELocationID
            elocation_ids = article.get("ELocationID", [])
            for elocation in elocation_ids:
                if hasattr(elocation, "attributes"):
                    if elocation.attributes.get("EIdType") == "doi":
                        return str(elocation)

            # Method 2: ArticleIdList in PubmedData
            if "PubmedData" in record:
                article_ids = record["PubmedData"].get("ArticleIdList", [])
                for article_id in article_ids:
                    if hasattr(article_id, "attributes"):
                        if article_id.attributes.get("IdType") == "doi":
                            return str(article_id)

        except Exception as e:
            logger.warning(f"Error extracting DOI: {e}")

        return None

    def _extract_keywords(self, record: Dict) -> List[str]:
        """Extract keywords from the record.

        Args:
            record: PubMed XML record

        Returns:
            List of keywords
        """
        keywords = []

        try:
            citation = record["MedlineCitation"]

            # Extract from KeywordList
            if "KeywordList" in citation:
                for keyword_list in citation["KeywordList"]:
                    for keyword in keyword_list:
                        keywords.append(str(keyword))

        except Exception as e:
            logger.warning(f"Error extracting keywords: {e}")

        return keywords

    def _extract_mesh_terms(self, record: Dict) -> List[str]:
        """Extract MeSH terms from the record.

        Args:
            record: PubMed XML record

        Returns:
            List of MeSH terms
        """
        mesh_terms = []

        try:
            citation = record["MedlineCitation"]

            # Extract from MeshHeadingList
            if "MeshHeadingList" in citation:
                for mesh_heading in citation["MeshHeadingList"]:
                    descriptor_name = mesh_heading.get("DescriptorName")
                    if descriptor_name:
                        mesh_terms.append(str(descriptor_name))

        except Exception as e:
            logger.warning(f"Error extracting MeSH terms: {e}")

        return mesh_terms

    def _clean_title(self, title: str) -> str:
        """Clean and normalize paper title.

        Args:
            title: Raw title string

        Returns:
            Cleaned title
        """
        if not title:
            return ""

        # Remove extra whitespace
        title = re.sub(r"\s+", " ", title.strip())

        # Remove trailing periods if they're not part of abbreviations
        if title.endswith(".") and not re.search(r"\b[A-Z]{2,}\.$", title):
            title = title[:-1]

        return title

    def _normalize_authors(self, authors: List[str]) -> List[str]:
        """Normalize author names.

        Args:
            authors: List of author names

        Returns:
            List of normalized author names
        """
        if not authors:
            return []

        normalized = []
        seen_authors = set()

        for author in authors:
            if not author or not author.strip():
                continue

            # Clean author name
            clean_author = re.sub(r"\s+", " ", author.strip())

            # Remove duplicates (case-insensitive)
            author_key = clean_author.lower()
            if author_key not in seen_authors:
                normalized.append(clean_author)
                seen_authors.add(author_key)

        return normalized

    def _normalize_journal_name(self, journal: str) -> str:
        """Normalize journal name.

        Args:
            journal: Raw journal name

        Returns:
            Normalized journal name
        """
        if not journal:
            return "Unknown Journal"

        # Clean whitespace
        journal = re.sub(r"\s+", " ", journal.strip())

        # Remove common suffixes that might vary
        suffixes_to_remove = [
            r"\s*\(Online\)$",
            r"\s*\(Print\)$",
            r"\s*:\s*official\s+.*$",
            r"\s*:\s*the\s+official\s+.*$",
        ]

        for suffix_pattern in suffixes_to_remove:
            journal = re.sub(suffix_pattern, "", journal, flags=re.IGNORECASE)

        return journal.strip()

    def _clean_abstract(self, abstract: Optional[str]) -> Optional[str]:
        """Clean and normalize abstract text.

        Args:
            abstract: Raw abstract text

        Returns:
            Cleaned abstract or None
        """
        if not abstract:
            return None

        # Remove extra whitespace
        abstract = re.sub(r"\s+", " ", abstract.strip())

        # Remove common prefixes
        prefixes_to_remove = [r"^Abstract:?\s*", r"^Background:?\s*", r"^Summary:?\s*"]

        for prefix_pattern in prefixes_to_remove:
            abstract = re.sub(prefix_pattern, "", abstract, flags=re.IGNORECASE)

        # Ensure minimum length
        if len(abstract) < 50:
            return None

        return abstract.strip()

    def calculate_relevance_score(self, paper: Paper, query: str) -> float:
        """Calculate relevance score between paper and search query.

        Args:
            paper: Paper object
            query: Search query string

        Returns:
            Relevance score between 0 and 1
        """
        if not query:
            return 0.0

        # Extract query terms
        query_terms = self._extract_query_terms(query)
        if not query_terms:
            return 0.0

        # Create searchable text from paper
        searchable_text = self._create_searchable_text(paper)

        # Calculate term matches
        matches = 0
        total_terms = len(query_terms)

        for term in query_terms:
            if term.lower() in searchable_text.lower():
                matches += 1

        # Calculate weighted score
        title_score = self._calculate_text_relevance(paper.title, query_terms)
        abstract_score = self._calculate_text_relevance(
            paper.abstract or "", query_terms
        )

        # Weighted combination
        relevance_score = (
            (matches / total_terms) * 0.4  # Term presence
            + title_score * 0.4  # Title relevance
            + abstract_score * 0.2  # Abstract relevance
        )

        return min(1.0, relevance_score)

    def _extract_query_terms(self, query: str) -> List[str]:
        """Extract meaningful terms from search query.

        Args:
            query: Search query string

        Returns:
            List of query terms
        """
        # Remove PubMed operators and common words
        operators = ["AND", "OR", "NOT", "(", ")", "[", "]", '"']
        stop_words = {"the", "a", "an", "in", "on", "at", "to", "for", "of", "with"}

        # Clean query
        clean_query = query
        for op in operators:
            clean_query = clean_query.replace(op, " ")

        # Extract terms
        terms = []
        for term in clean_query.split():
            term = term.strip().lower()
            if len(term) > 2 and term not in stop_words:
                terms.append(term)

        return terms

    def _create_searchable_text(self, paper: Paper) -> str:
        """Create searchable text from paper data.

        Args:
            paper: Paper object

        Returns:
            Combined searchable text
        """
        parts = [
            paper.title,
            " ".join(paper.authors),
            paper.journal,
            paper.abstract or "",
        ]

        return " ".join(part for part in parts if part)

    def _calculate_text_relevance(self, text: str, query_terms: List[str]) -> float:
        """Calculate relevance score for a text field.

        Args:
            text: Text to analyze
            query_terms: List of query terms

        Returns:
            Relevance score between 0 and 1
        """
        if not text or not query_terms:
            return 0.0

        text_lower = text.lower()
        matches = sum(1 for term in query_terms if term in text_lower)

        return matches / len(query_terms)
