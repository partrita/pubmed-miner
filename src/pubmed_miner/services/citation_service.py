"""
Citation information collection service.
"""

import logging
import time
from datetime import datetime, timedelta
from typing import List, Optional, Tuple, Dict

import requests
from urllib3.util.retry import Retry

from ..models.cache import CitationCache

logger = logging.getLogger(__name__)

try:
    from Bio import Entrez as BioEntrez
    HAS_BIOENTRZ = True
except ImportError:
    BioEntrez = None
    HAS_BIOENTRZ = False


class CacheManager:
    """Manages citation count caching."""

    def __init__(self, cache_expiry_days: int = 7, max_cache_size: int = 5000):
        self.cache_expiry_days = cache_expiry_days
        self.max_cache_size = max_cache_size
        self.cache: Dict[str, CitationCache] = {}

    def get_citation(
        self, pmid: str, doi: Optional[str] = None
    ) -> Optional[CitationCache]:
        key = pmid
        if doi:
            key = f"{pmid}:{doi}"
        return self.cache.get(key)

    def save_citation(
        self, pmid: str, count: int, source: str, doi: Optional[str] = None
    ) -> None:
        key = pmid
        if doi:
            key = f"{pmid}:{doi}"
        self.cache[key] = CitationCache(
            pmid=pmid,
            citation_count=count,
            last_updated=datetime.now(),
            source=source,
        )

    def is_expired(self, citation: CitationCache) -> bool:
        return citation.is_expired(self.cache_expiry_days)

    def clear_expired(self) -> int:
        expired_keys = [
            k for k, v in self.cache.items() if self.is_expired(v)
        ]
        for k in expired_keys:
            del self.cache[k]
        return len(expired_keys)

    clear_expired_citations = clear_expired

    clear_expired_citations = clear_expired

    def get_all_papers(self) -> List[CitationCache]:
        return list(self.cache.values())

    def cleanup_unused(self, max_age_days: int = 30) -> int:
        cutoff = datetime.now() - timedelta(days=max_age_days)
        old_keys = [k for k, v in self.cache.items() if v.last_updated < cutoff]
        for k in old_keys:
            del self.cache[k]
        return len(old_keys)


class CitationService:
    """Service for collecting citation information for papers."""

    def __init__(
        self,
        cache_expiry_days: int = 7,
        max_cache_size: int = 5000,
        retry_attempts: int = 3,
        retry_delay: float = 1.0,
        crossref_email: Optional[str] = None,
        cache_manager: Optional[CacheManager] = None,
    ):
        self.logger = logging.getLogger(__name__)
        self.cache_manager = cache_manager or CacheManager(
            cache_expiry_days=cache_expiry_days,
            max_cache_size=max_cache_size,
        )
        self.crossref_base_url = "https://api.crossref.org/works"
        self.semantic_scholar_base_url = (
            "https://api.semanticscholar.org/graph/v1/paper"
        )
        self.retry_attempts = retry_attempts
        self.retry_delay = retry_delay
        self.session = None

        self._total_requests = 0
        self._cache_hits = 0
        self._cache_misses = 0
        self._api_calls = 0
        self._crossref_calls = 0
        self._semantic_scholar_calls = 0
        self._error_count = 0
        self._total_response_time = 0.0

        self.logger.info("Initialized CitationService")

    def get_citation_count(
        self, pmid: str, doi: Optional[str] = None
    ) -> int:
        """Get citation count for a paper using PMID or DOI."""
        if not pmid:
            raise ValueError("PMID cannot be empty")
        if not pmid.isdigit():
            raise ValueError("PMID must be numeric")

        self._total_requests += 1
        self._api_calls += 1

        cached = self.cache_manager.get_citation(pmid, doi)
        if cached and not self._is_cache_expired(cached):
            self._cache_hits += 1
            return cached.citation_count

        self._cache_misses += 1

        for attempt in range(self.retry_attempts):
            try:
                count, source = self._fetch_from_apis(pmid, doi)
                if count is not None:
                    self.cache_manager.save_citation(
                        pmid, count, source, doi
                    )
                    return count
            except requests.exceptions.HTTPError as e:
                if e.response is not None and e.response.status_code == 429:
                    retry_after = e.response.headers.get("Retry-After", "60")
                    try:
                        sleep_time = int(retry_after)
                    except (ValueError, TypeError):
                        sleep_time = 60
                    time.sleep(sleep_time)
                    continue
                self._error_count += 1
            except Exception as e:
                self._error_count += 1
                self.logger.error(
                    f"Failed to fetch citation count for PMID {pmid}: {e}"
                )

        return 0

    def _fetch_from_apis(
        self, pmid: str, doi: Optional[str] = None
    ) -> Tuple[Optional[int], str]:
        if doi:
            count = self._fetch_from_crossref(doi)
            if count is not None:
                return count, "crossref"

        count = self._fetch_from_semantic_scholar(pmid, doi)
        if count is not None:
            return count, "semantic_scholar"

        return None, None

    def _fetch_from_crossref(self, doi: str) -> Optional[int]:
        try:
            response = requests.get(
                f"{self.crossref_base_url}/{doi}",
                timeout=30,
            )
            response.raise_for_status()
            data = response.json()
            if isinstance(data, dict) and "message" in data:
                count = data["message"].get("is-referenced-by-count", 0)
            elif isinstance(data, dict):
                count = data.get("is-referenced-by-count", 0)
            else:
                count = 0
            self._crossref_calls += 1
            return count
        except requests.exceptions.RequestException:
            return None
        except Exception:
            return None

    def _fetch_from_semantic_scholar(
        self, pmid: str, doi: Optional[str] = None
    ) -> int:
        try:
            identifier = doi if doi else f"pubmed:{pmid}"
            url = (
                f"{self.semantic_scholar_base_url}/"
                f"{identifier}?fields=citationCount"
            )
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            count = data.get("citationCount", 0)
            self._semantic_scholar_calls += 1
            return count
        except requests.exceptions.RequestException:
            return 0
        except Exception:
            return 0

    def _is_cache_expired(self, citation: CitationCache) -> bool:
        return self.cache_manager.is_expired(citation)

    def batch_get_citation_counts(
        self, pmids: List[str], skip_errors: bool = False
    ) -> Dict[str, int]:
        results = {}
        for pmid in pmids:
            try:
                count = self.get_citation_count(pmid)
                results[pmid] = count
            except (ValueError, Exception):
                if skip_errors:
                    continue
                else:
                    raise
        return results

    def update_cache_settings(self, settings: Dict) -> None:
        if "cache_expiry_days" in settings:
            self.cache_manager.cache_expiry_days = settings[
                "cache_expiry_days"
            ]
        if "max_cache_size" in settings:
            self.cache_manager.max_cache_size = settings[
                "max_cache_size"
            ]

    def get_cached_papers(self) -> List[CitationCache]:
        return self.cache_manager.get_all_papers()

    def clear_expired_cache(self) -> int:
        return self.cache_manager.clear_expired()

    def get_statistics(self) -> Dict:
        return {
            "total_requests": self._total_requests,
            "cache_hits": self._cache_hits,
            "cache_misses": self._cache_misses,
            "api_calls": self._api_calls,
            "crossref_calls": self._crossref_calls,
            "semantic_scholar_calls": self._semantic_scholar_calls,
            "error_count": self._error_count,
            "average_response_time": self._total_response_time,
        }

    def _validate_doi(self, doi: str) -> bool:
        if not doi or not isinstance(doi, str):
            return False
        return doi.startswith("10.") and "/" in doi

    def _get_doi_from_pmid(self, pmid: str) -> Optional[str]:
        try:
            if not HAS_BIOENTRZ:
                self.logger.warning(
                    f"Bio.Entrez not available, cannot get DOI for PMID {pmid}"
                )
                return None

            Entrez = BioEntrez
            Entrez.email = "pubmed.miner@example.com"
            handle = Entrez.esearch(
                db="pubmed", term=pmid, retmode="xml"
            )
            record = Entrez.read(handle)
            handle.close()

            pmids = record.get("IdList", [])
            if not pmids:
                return None

            fetch_handle = Entrez.efetch(
                db="pubmed",
                id=pmid,
                rettype="medline",
                retmode="xml",
            )
            article_data = Entrez.read(fetch_handle)
            fetch_handle.close()

            article_ids = article_data.get("PubmedArticleSet", [])
            if isinstance(article_ids, dict):
                article_ids = [article_ids]
            for article in article_ids:
                id_list = article.get("PubmedData", {}).get(
                    "ArticleIdList", []
                )
                if isinstance(id_list, dict):
                    id_list = [id_list]
                for id_elem in id_list:
                    if (
                        isinstance(id_elem, dict)
                        and id_elem.get("IdType") == "doi"
                    ):
                        doi = id_elem.get("Id") or id_elem.get("#text", "")
                        if doi:
                            return doi

            return None
        except Exception as e:
            self.logger.warning(
                f"Failed to get DOI for PMID {pmid}: {e}"
            )
            return None

    def _clean_html(self, text: str) -> str:
        import re

        clean = re.sub(r"<[^>]+>", "", text)
        clean = clean.replace("&nbsp;", " ")
        clean = clean.replace("&amp;", "&")
        clean = clean.replace("&lt;", "<")
        clean = clean.replace("&gt;", ">")
        clean = clean.replace("&quot;", '"')
        clean = clean.strip()
        return clean

    def cleanup_unused_cache_entries(self, max_age_days: int = 30) -> int:
        return self.cache_manager.cleanup_unused(max_age_days)
