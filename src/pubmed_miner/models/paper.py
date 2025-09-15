"""
Data models for paper-related structures.
"""
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional


@dataclass
class Paper:
    """Core paper data model."""
    pmid: str
    title: str
    authors: List[str]
    journal: str
    publication_date: datetime
    abstract: Optional[str] = None
    doi: Optional[str] = None

    def __post_init__(self):
        """Validate paper data after initialization."""
        if not self.pmid:
            raise ValueError("PMID cannot be empty")
        if not self.title:
            raise ValueError("Title cannot be empty")
        if not self.authors:
            raise ValueError("Authors list cannot be empty")
        if not self.journal:
            raise ValueError("Journal cannot be empty")


@dataclass
class ScoredPaper(Paper):
    """Paper with scoring information."""
    citation_count: int = 0
    impact_factor: float = 0.0
    score: float = 0.0
    rank: int = 0

    def __post_init__(self):
        """Validate scored paper data."""
        super().__post_init__()
        if self.citation_count < 0:
            raise ValueError("Citation count cannot be negative")
        if self.impact_factor < 0:
            raise ValueError("Impact factor cannot be negative")
        if self.score < 0:
            raise ValueError("Score cannot be negative")
        if self.rank < 0:
            raise ValueError("Rank cannot be negative")