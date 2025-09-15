"""
Cache data models.
"""
from dataclasses import dataclass
from datetime import datetime
from typing import Optional


@dataclass
class CitationCache:
    """Cache entry for citation information."""
    pmid: str
    citation_count: int
    last_updated: datetime
    source: str

    def __post_init__(self):
        """Validate citation cache data."""
        if not self.pmid:
            raise ValueError("PMID cannot be empty")
        if self.citation_count < 0:
            raise ValueError("Citation count cannot be negative")
        if not self.source:
            raise ValueError("Source cannot be empty")

    def is_expired(self, max_age_days: int = 7) -> bool:
        """Check if cache entry is expired."""
        age = datetime.now() - self.last_updated
        return age.days > max_age_days


@dataclass
class ImpactFactorCache:
    """Cache entry for journal impact factor."""
    journal_name: str
    impact_factor: float
    year: int
    last_updated: datetime
    source: str

    def __post_init__(self):
        """Validate impact factor cache data."""
        if not self.journal_name:
            raise ValueError("Journal name cannot be empty")
        if self.impact_factor < 0:
            raise ValueError("Impact factor cannot be negative")
        if self.year < 1900 or self.year > datetime.now().year:
            raise ValueError("Year must be reasonable")
        if not self.source:
            raise ValueError("Source cannot be empty")

    def is_expired(self, max_age_days: int = 365) -> bool:
        """Check if cache entry is expired."""
        age = datetime.now() - self.last_updated
        return age.days > max_age_days


@dataclass
class PaperMetadataCache:
    """Cache entry for paper metadata."""
    pmid: str
    title: str
    authors_json: str  # JSON serialized list of authors
    journal: str
    publication_date: datetime
    abstract: Optional[str]
    doi: Optional[str]
    last_updated: datetime

    def __post_init__(self):
        """Validate paper metadata cache."""
        if not self.pmid:
            raise ValueError("PMID cannot be empty")
        if not self.title:
            raise ValueError("Title cannot be empty")
        if not self.journal:
            raise ValueError("Journal cannot be empty")

    def is_expired(self, max_age_days: int = 30) -> bool:
        """Check if cache entry is expired."""
        age = datetime.now() - self.last_updated
        return age.days > max_age_days