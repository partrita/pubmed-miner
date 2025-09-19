"""
Journal impact factor collection and estimation service.
"""

import re
import logging
from typing import Dict, Optional
from datetime import datetime
from difflib import SequenceMatcher

from ..models.cache import ImpactFactorCache
from ..utils.error_handler import retry_api_calls

logger = logging.getLogger(__name__)


class ImpactFactorService:
    """Service for collecting and estimating journal impact factors."""

    def __init__(self, cache_manager=None):
        """Initialize impact factor service.

        Args:
            cache_manager: Optional cache manager for storing impact factor data
        """
        self.cache_manager = cache_manager

        # Built-in impact factor data for common journals (2023 data)
        self.known_impact_factors = {
            "nature": 64.8,
            "science": 56.9,
            "cell": 64.5,
            "new england journal of medicine": 176.1,
            "lancet": 168.9,
            "nature medicine": 87.2,
            "nature genetics": 41.3,
            "nature biotechnology": 68.2,
            "plos one": 3.7,
            "scientific reports": 4.6,
            "proceedings of the national academy of sciences": 12.8,
            "journal of biological chemistry": 4.8,
            "nucleic acids research": 14.9,
            "bioinformatics": 6.9,
            "bmc bioinformatics": 3.3,
            "genome biology": 17.9,
            "nature communications": 17.7,
            "plos biology": 9.8,
            "molecular biology and evolution": 16.2,
            "genome research": 7.0,
        }

        logger.info("Initialized ImpactFactorService")

    def get_impact_factor(self, journal_name: str) -> Optional[float]:
        """Get impact factor for a journal.

        Args:
            journal_name: Journal name

        Returns:
            Impact factor or None if not found
        """
        if not journal_name or not journal_name.strip():
            return None

        # Normalize journal name
        normalized_name = self._normalize_journal_name(journal_name)

        # Check cache first
        if self.cache_manager:
            cached_if = self._get_cached_impact_factor(normalized_name)
            if cached_if is not None:
                return cached_if

        # Try built-in data
        impact_factor = self._get_builtin_impact_factor(normalized_name)
        if impact_factor is not None:
            self._cache_impact_factor(normalized_name, impact_factor, "builtin")
            return impact_factor

        # Try external APIs
        try:
            impact_factor = self._get_scimago_impact_factor(normalized_name)
            if impact_factor is not None:
                self._cache_impact_factor(normalized_name, impact_factor, "scimago")
                return impact_factor
        except Exception as e:
            logger.warning(f"SCImago lookup failed for {journal_name}: {e}")

        # Try estimation based on similar journals
        estimated_if = self.estimate_impact_factor(journal_name)
        if estimated_if is not None:
            self._cache_impact_factor(normalized_name, estimated_if, "estimated")
            return estimated_if

        logger.warning(f"No impact factor found for journal: {journal_name}")
        return None

    def estimate_impact_factor(self, journal_name: str) -> float:
        """Estimate impact factor based on journal name patterns and similar journals.

        Args:
            journal_name: Journal name

        Returns:
            Estimated impact factor
        """
        normalized_name = self._normalize_journal_name(journal_name)

        # Pattern-based estimation
        pattern_if = self._estimate_by_patterns(normalized_name)
        if pattern_if is not None:
            logger.info(f"Pattern-based estimate for {journal_name}: {pattern_if}")
            return pattern_if

        # Similarity-based estimation
        similar_if = self._estimate_by_similarity(normalized_name)
        if similar_if is not None:
            logger.info(f"Similarity-based estimate for {journal_name}: {similar_if}")
            return similar_if

        # Default estimation based on journal type
        default_if = self._estimate_by_type(normalized_name)
        logger.info(f"Default estimate for {journal_name}: {default_if}")
        return default_if

    def match_journal_name(self, partial_name: str) -> Optional[str]:
        """Find the best matching journal name from known journals.

        Args:
            partial_name: Partial or abbreviated journal name

        Returns:
            Best matching full journal name or None
        """
        if not partial_name:
            return None

        normalized_partial = self._normalize_journal_name(partial_name)
        best_match = None
        best_score = 0.0

        # Check against known journals
        for known_journal in self.known_impact_factors.keys():
            score = SequenceMatcher(None, normalized_partial, known_journal).ratio()
            if score > best_score and score > 0.6:  # Minimum similarity threshold
                best_score = score
                best_match = known_journal

        # Check cache for additional matches
        if self.cache_manager and best_score < 0.8:
            cached_matches = self.cache_manager.search_similar_journals(
                normalized_partial
            )
            for cached_journal in cached_matches:
                score = SequenceMatcher(
                    None,
                    normalized_partial,
                    self._normalize_journal_name(cached_journal.journal_name),
                ).ratio()
                if score > best_score:
                    best_score = score
                    best_match = cached_journal.journal_name

        return best_match if best_score > 0.6 else None

    def cache_impact_factor(self, journal: str, factor: float) -> None:
        """Manually cache an impact factor.

        Args:
            journal: Journal name
            factor: Impact factor value
        """
        self._cache_impact_factor(journal, factor, "manual")

    def _normalize_journal_name(self, journal_name: str) -> str:
        """Normalize journal name for consistent matching.

        Args:
            journal_name: Raw journal name

        Returns:
            Normalized journal name
        """
        if not journal_name:
            return ""

        # Convert to lowercase
        normalized = journal_name.lower().strip()

        # Remove common prefixes and suffixes
        patterns_to_remove = [
            r"^the\s+",
            r"\s*\(online\)$",
            r"\s*\(print\)$",
            r"\s*:\s*official.*$",
            r"\s*:\s*the\s+official.*$",
            r"\s*-\s*official.*$",
        ]

        for pattern in patterns_to_remove:
            normalized = re.sub(pattern, "", normalized, flags=re.IGNORECASE)

        # Standardize common abbreviations
        abbreviations = {
            r"\bj\b": "journal",
            r"\bproc\b": "proceedings",
            r"\bnat\b": "nature",
            r"\bsci\b": "science",
            r"\bmed\b": "medicine",
            r"\bbiol\b": "biology",
            r"\bchem\b": "chemistry",
            r"\bphys\b": "physics",
            r"\bint\b": "international",
            r"\bam\b": "american",
            r"\beur\b": "european",
        }

        for abbrev, full in abbreviations.items():
            normalized = re.sub(abbrev, full, normalized)

        # Clean up extra spaces
        normalized = re.sub(r"\s+", " ", normalized).strip()

        return normalized

    def _get_builtin_impact_factor(self, normalized_name: str) -> Optional[float]:
        """Get impact factor from built-in data.

        Args:
            normalized_name: Normalized journal name

        Returns:
            Impact factor or None
        """
        # Direct match
        if normalized_name in self.known_impact_factors:
            return self.known_impact_factors[normalized_name]

        # Fuzzy match
        for known_journal, impact_factor in self.known_impact_factors.items():
            similarity = SequenceMatcher(None, normalized_name, known_journal).ratio()
            if similarity > 0.9:  # Very high similarity threshold for built-in data
                logger.debug(
                    f"Fuzzy match: {normalized_name} -> {known_journal} (similarity: {similarity:.2f})"
                )
                return impact_factor

        return None

    @retry_api_calls(max_attempts=2, delay=1.0)
    def _get_scimago_impact_factor(self, journal_name: str) -> Optional[float]:
        """Get impact factor from SCImago API (placeholder implementation).

        Args:
            journal_name: Journal name

        Returns:
            Impact factor or None
        """
        # Note: This is a placeholder implementation
        # SCImago doesn't have a public API, so this would need to be implemented
        # using web scraping or a different data source

        logger.debug(f"SCImago lookup not implemented for {journal_name}")
        return None

    def _estimate_by_patterns(self, normalized_name: str) -> Optional[float]:
        """Estimate impact factor based on journal name patterns.

        Args:
            normalized_name: Normalized journal name

        Returns:
            Estimated impact factor or None
        """
        # High-impact patterns
        high_impact_patterns = [
            r"nature\s+",
            r"science\s+",
            r"cell\s+",
            r"new\s+england\s+journal",
            r"lancet",
            r"journal\s+of\s+clinical\s+investigation",
        ]

        for pattern in high_impact_patterns:
            if re.search(pattern, normalized_name):
                return 25.0  # High estimate for prestigious journals

        # Medium-impact patterns
        medium_impact_patterns = [
            r"proceedings\s+of\s+the\s+national\s+academy",
            r"journal\s+of\s+biological\s+chemistry",
            r"molecular\s+biology",
            r"genome\s+research",
            r"nucleic\s+acids\s+research",
        ]

        for pattern in medium_impact_patterns:
            if re.search(pattern, normalized_name):
                return 8.0  # Medium estimate

        # Open access patterns (typically lower IF)
        open_access_patterns = [
            r"plos\s+",
            r"bmc\s+",
            r"frontiers\s+in",
            r"scientific\s+reports",
        ]

        for pattern in open_access_patterns:
            if re.search(pattern, normalized_name):
                return 4.0  # Lower estimate for open access

        return None

    def _estimate_by_similarity(self, normalized_name: str) -> Optional[float]:
        """Estimate impact factor based on similarity to known journals.

        Args:
            normalized_name: Normalized journal name

        Returns:
            Estimated impact factor or None
        """
        best_similarity = 0.0
        best_impact_factor = None

        for known_journal, impact_factor in self.known_impact_factors.items():
            similarity = SequenceMatcher(None, normalized_name, known_journal).ratio()
            if similarity > best_similarity and similarity > 0.5:
                best_similarity = similarity
                best_impact_factor = impact_factor

        if best_impact_factor is not None:
            # Adjust estimate based on similarity
            adjustment_factor = 0.5 + (best_similarity * 0.5)  # 0.5 to 1.0
            return best_impact_factor * adjustment_factor

        return None

    def _estimate_by_type(self, normalized_name: str) -> float:
        """Provide default estimate based on journal type.

        Args:
            normalized_name: Normalized journal name

        Returns:
            Default impact factor estimate
        """
        # Medical journals tend to have higher impact factors
        if any(
            term in normalized_name
            for term in ["medicine", "medical", "clinical", "therapy"]
        ):
            return 6.0

        # Biology and life sciences
        if any(
            term in normalized_name
            for term in ["biology", "biological", "life", "molecular"]
        ):
            return 4.5

        # Computer science and bioinformatics
        if any(
            term in normalized_name
            for term in ["bioinformatics", "computational", "computer"]
        ):
            return 3.5

        # General science journals
        if any(term in normalized_name for term in ["science", "research", "journal"]):
            return 3.0

        # Default for unknown types
        return 2.5

    def _get_cached_impact_factor(self, journal_name: str) -> Optional[float]:
        """Get impact factor from cache.

        Args:
            journal_name: Journal name

        Returns:
            Cached impact factor or None
        """
        if not self.cache_manager:
            return None

        try:
            cache_entry = self.cache_manager.get_impact_factor(journal_name)
            if cache_entry and not cache_entry.is_expired():
                logger.debug(
                    f"Using cached impact factor for {journal_name}: {cache_entry.impact_factor}"
                )
                return cache_entry.impact_factor
        except Exception as e:
            logger.warning(f"Error reading impact factor cache for {journal_name}: {e}")

        return None

    def _cache_impact_factor(
        self, journal_name: str, impact_factor: float, source: str
    ) -> None:
        """Cache impact factor.

        Args:
            journal_name: Journal name
            impact_factor: Impact factor value
            source: Data source name
        """
        if not self.cache_manager:
            return

        try:
            cache_entry = ImpactFactorCache(
                journal_name=journal_name,
                impact_factor=impact_factor,
                year=datetime.now().year,
                last_updated=datetime.now(),
                source=source,
            )
            self.cache_manager.save_impact_factor(cache_entry)
            logger.debug(
                f"Cached impact factor for {journal_name}: {impact_factor} from {source}"
            )
        except Exception as e:
            logger.warning(f"Error caching impact factor for {journal_name}: {e}")

    def get_journal_statistics(self) -> Dict[str, any]:
        """Get statistics about cached journal data.

        Returns:
            Dictionary with journal statistics
        """
        stats = {
            "builtin_journals": len(self.known_impact_factors),
            "cached_journals": 0,
            "average_impact_factor": 0.0,
            "top_journals": [],
        }

        # Calculate average of built-in impact factors
        if self.known_impact_factors:
            stats["average_impact_factor"] = sum(
                self.known_impact_factors.values()
            ) / len(self.known_impact_factors)

            # Get top 5 journals by impact factor
            sorted_journals = sorted(
                self.known_impact_factors.items(), key=lambda x: x[1], reverse=True
            )
            stats["top_journals"] = sorted_journals[:5]

        # Add cache statistics if available
        if self.cache_manager:
            try:
                cache_stats = self.cache_manager.get_impact_factor_cache_stats()
                stats["cached_journals"] = cache_stats.get("total", 0)
            except Exception as e:
                logger.warning(f"Error getting cache statistics: {e}")

        return stats
