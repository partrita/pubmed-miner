"""
Citation information collection service.
"""

import requests
import time
import logging
from typing import Dict, List, Optional, Union
from datetime import datetime

from ..models.cache import CitationCache
from ..utils.error_handler import retry_api_calls

logger = logging.getLogger(__name__)


class CitationService:
    """Service for collecting citation information from multiple sources."""

    def __init__(self, cache_manager=None):
        """Initialize citation service.

        Args:
            cache_manager: Optional cache manager for storing citation data
        """
        self.cache_manager = cache_manager
        self.last_request_times = {"crossref": 0.0, "semantic_scholar": 0.0, "pmc": 0.0}

        # Rate limits (requests per second)
        self.rate_limits = {"crossref": 50, "semantic_scholar": 100, "pmc": 3}

        logger.info("Initialized CitationService")

    def get_citation_count(self, pmid: str) -> Optional[int]:
        """Get citation count for a single PMID.

        Args:
            pmid: PubMed ID

        Returns:
            Citation count or None if not found
        """
        # Check cache first
        if self.cache_manager:
            cached_count = self._get_cached_citation(pmid)
            if cached_count is not None:
                return cached_count

        # Try multiple sources
        citation_count = None

        # Try Crossref first (most reliable)
        try:
            citation_count = self._get_crossref_citations(pmid)
            if citation_count is not None:
                self._cache_citation(pmid, citation_count, "crossref")
                return citation_count
        except Exception as e:
            logger.warning(f"Crossref failed for {pmid}: {e}")

        # Try Semantic Scholar as fallback
        try:
            citation_count = self._get_semantic_scholar_citations(pmid)
            if citation_count is not None:
                self._cache_citation(pmid, citation_count, "semantic_scholar")
                return citation_count
        except Exception as e:
            logger.warning(f"Semantic Scholar failed for {pmid}: {e}")

        # Try PMC as last resort
        try:
            citation_count = self._get_pmc_citations(pmid)
            if citation_count is not None:
                self._cache_citation(pmid, citation_count, "pmc")
                return citation_count
        except Exception as e:
            logger.warning(f"PMC failed for {pmid}: {e}")

        logger.warning(f"No citation data found for PMID {pmid}")
        return None

    def batch_get_citations(
        self, pmids: List[str], skip_errors: bool = True
    ) -> Dict[str, int]:
        """Get citation counts for multiple PMIDs.

        Args:
            pmids: List of PubMed IDs
            skip_errors: If True, skip PMIDs that fail and continue processing

        Returns:
            Dictionary mapping PMID to citation count (None values excluded)
        """
        citations = {}

        # Check cache for all PMIDs first
        uncached_pmids = []
        if self.cache_manager:
            for pmid in pmids:
                cached_count = self._get_cached_citation(pmid)
                if cached_count is not None:
                    citations[pmid] = cached_count
                else:
                    uncached_pmids.append(pmid)
        else:
            uncached_pmids = pmids

        if not uncached_pmids:
            return citations

        logger.info(f"Fetching citations for {len(uncached_pmids)} uncached PMIDs")

        # Process in batches to respect rate limits
        batch_size = 10
        for i in range(0, len(uncached_pmids), batch_size):
            batch = uncached_pmids[i : i + batch_size]
            batch_citations = self._fetch_citation_batch(batch, skip_errors=skip_errors)
            citations.update(batch_citations)

            # Rate limiting between batches
            if i + batch_size < len(uncached_pmids):
                time.sleep(0.5)

        logger.info(f"Retrieved citations for {len(citations)} PMIDs")
        return citations

    def _fetch_citation_batch(
        self, pmids: List[str], skip_errors: bool = True
    ) -> Dict[str, int]:
        """Fetch citations for a batch of PMIDs.

        Args:
            pmids: List of PubMed IDs
            skip_errors: If True, skip PMIDs that fail and continue processing

        Returns:
            Dictionary mapping PMID to citation count (None values excluded)
        """
        citations = {}

        for pmid in pmids:
            try:
                count = self.get_citation_count(pmid)
                if count is not None:
                    # Ensure count is a valid integer
                    safe_count = self._safe_citation_count(count)
                    citations[pmid] = safe_count
                elif not skip_errors:
                    # If not skipping errors, set to 0 as fallback
                    citations[pmid] = 0
            except Exception as e:
                logger.warning(f"Failed to get citations for {pmid}: {e}")
                if not skip_errors:
                    citations[pmid] = 0
                continue

        return citations

    @retry_api_calls(max_attempts=3, delay=1.0)
    def _get_crossref_citations(self, pmid: str) -> Optional[int]:
        """Get citation count from Crossref API.

        Args:
            pmid: PubMed ID

        Returns:
            Citation count or None
        """
        self._rate_limit("crossref")

        try:
            # First, get DOI from PubMed
            doi = self._get_doi_for_pmid(pmid)
            if not doi:
                return None

            # Query Crossref for citation count
            url = f"https://api.crossref.org/works/{doi}"
            headers = {
                "User-Agent": "PubMedMiner/1.0 (mailto:pubmed.miner@example.com)"
            }

            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()

            data = response.json()
            citation_count = data.get("message", {}).get("is-referenced-by-count", 0)

            logger.debug(f"Crossref citations for {pmid}: {citation_count}")
            return citation_count

        except requests.exceptions.RequestException as e:
            logger.warning(f"Crossref API error for {pmid}: {e}")
            return None
        except Exception as e:
            logger.warning(f"Error getting Crossref citations for {pmid}: {e}")
            return None

    @retry_api_calls(max_attempts=3, delay=1.0)
    def _get_semantic_scholar_citations(self, pmid: str) -> Optional[int]:
        """Get citation count from Semantic Scholar API.

        Args:
            pmid: PubMed ID

        Returns:
            Citation count or None
        """
        self._rate_limit("semantic_scholar")

        try:
            url = f"https://api.semanticscholar.org/graph/v1/paper/PMID:{pmid}"
            params = {"fields": "citationCount"}

            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()

            data = response.json()
            citation_count = data.get("citationCount", 0)

            logger.debug(f"Semantic Scholar citations for {pmid}: {citation_count}")
            return citation_count

        except requests.exceptions.RequestException as e:
            if response.status_code == 404:
                logger.debug(f"Paper {pmid} not found in Semantic Scholar")
                return None
            logger.warning(f"Semantic Scholar API error for {pmid}: {e}")
            return None
        except Exception as e:
            logger.warning(f"Error getting Semantic Scholar citations for {pmid}: {e}")
            return None

    def _get_pmc_citations(self, pmid: str) -> Optional[int]:
        """Get citation count from PMC (placeholder implementation).

        Args:
            pmid: PubMed ID

        Returns:
            Citation count or None
        """
        # PMC doesn't have a direct citation API, so this is a placeholder
        # In a real implementation, you might scrape PMC or use other methods
        logger.debug(f"PMC citation lookup not implemented for {pmid}")
        return None

    def _get_doi_for_pmid(self, pmid: str) -> Optional[str]:
        """Get DOI for a PMID using PubMed API.

        Args:
            pmid: PubMed ID

        Returns:
            DOI string or None
        """
        try:
            from Bio import Entrez

            handle = Entrez.efetch(
                db="pubmed", id=pmid, rettype="medline", retmode="xml"
            )

            records = Entrez.read(handle)
            handle.close()

            if not records["PubmedArticle"]:
                return None

            record = records["PubmedArticle"][0]

            # Try multiple methods to extract DOI
            # Method 1: ELocationID
            article = record["MedlineCitation"]["Article"]
            elocation_ids = article.get("ELocationID", [])

            for elocation in elocation_ids:
                if hasattr(elocation, "attributes"):
                    if elocation.attributes.get("EIdType") == "doi":
                        return str(elocation)

            # Method 2: ArticleIdList
            if "PubmedData" in record:
                article_ids = record["PubmedData"].get("ArticleIdList", [])
                for article_id in article_ids:
                    if hasattr(article_id, "attributes"):
                        if article_id.attributes.get("IdType") == "doi":
                            return str(article_id)

            return None

        except Exception as e:
            logger.warning(f"Error getting DOI for {pmid}: {e}")
            return None

    def _get_cached_citation(self, pmid: str) -> Optional[int]:
        """Get citation count from cache.

        Args:
            pmid: PubMed ID

        Returns:
            Cached citation count or None
        """
        if not self.cache_manager:
            return None

        try:
            cache_entry = self.cache_manager.get_citation(pmid)
            if cache_entry and not cache_entry.is_expired():
                logger.debug(
                    f"Using cached citation for {pmid}: {cache_entry.citation_count}"
                )
                return cache_entry.citation_count
        except Exception as e:
            logger.warning(f"Error reading citation cache for {pmid}: {e}")

        return None

    def _cache_citation(self, pmid: str, count: int, source: str) -> None:
        """Cache citation count.

        Args:
            pmid: PubMed ID
            count: Citation count
            source: Data source name
        """
        if not self.cache_manager:
            return

        try:
            cache_entry = CitationCache(
                pmid=pmid,
                citation_count=count,
                last_updated=datetime.now(),
                source=source,
            )
            self.cache_manager.save_citation(cache_entry)
            logger.debug(f"Cached citation for {pmid}: {count} from {source}")
        except Exception as e:
            logger.warning(f"Error caching citation for {pmid}: {e}")

    def _rate_limit(self, service: str) -> None:
        """Apply rate limiting for a service.

        Args:
            service: Service name ('crossref', 'semantic_scholar', 'pmc')
        """
        if service not in self.rate_limits:
            return

        current_time = time.time()
        last_request = self.last_request_times.get(service, 0)
        min_interval = 1.0 / self.rate_limits[service]

        time_since_last = current_time - last_request
        if time_since_last < min_interval:
            sleep_time = min_interval - time_since_last
            time.sleep(sleep_time)

        self.last_request_times[service] = time.time()

    def update_citation_cache(self, pmid: str, count: int) -> None:
        """Update citation count in cache.

        Args:
            pmid: PubMed ID
            count: New citation count
        """
        self._cache_citation(pmid, count, "manual_update")

    def clear_expired_cache(self, max_age_days: int = 7) -> int:
        """Clear expired cache entries.

        Args:
            max_age_days: Maximum age for cache entries in days

        Returns:
            Number of entries cleared
        """
        if not self.cache_manager:
            return 0

        try:
            cleared_count = self.cache_manager.clear_expired_citations(max_age_days)
            logger.info(f"Cleared {cleared_count} expired citation cache entries")
            return cleared_count
        except Exception as e:
            logger.error(f"Error clearing expired cache: {e}")
            return 0

    def get_cache_stats(self) -> Dict[str, int]:
        """Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        if not self.cache_manager:
            return {"total": 0, "expired": 0}

        try:
            return self.cache_manager.get_citation_cache_stats()
        except Exception as e:
            logger.error(f"Error getting cache stats: {e}")
            return {"total": 0, "expired": 0}

    def _safe_citation_count(self, count: Union[int, str, None]) -> int:
        """Safely convert citation count to integer, handling None and invalid values.

        Args:
            count: Citation count value to convert

        Returns:
            Safe integer citation count (>= 0)
        """
        if count is None:
            return 0

        try:
            int_count = int(count)
            return max(0, int_count)  # Ensure non-negative
        except (ValueError, TypeError):
            logger.debug(f"Failed to convert citation count {count} to int, using 0")
            return 0
