"""
Configuration data models.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any


@dataclass
class TopicConfig:
    """Configuration for a research topic."""

    name: str
    query: str
    max_papers: int = 1000
    essential_count: int = 15
    enabled: bool = True

    def __post_init__(self):
        """Validate topic configuration."""
        if not self.name:
            raise ValueError("Topic name cannot be empty")
        if not self.query:
            raise ValueError("Query cannot be empty")
        if self.max_papers <= 0:
            raise ValueError("Max papers must be positive")
        if self.essential_count <= 0:
            raise ValueError("Essential count must be positive")
        if self.essential_count > self.max_papers:
            raise ValueError("Essential count cannot exceed max papers")


@dataclass
class GitHubConfig:
    """GitHub integration configuration."""

    token: str
    repository: str
    issue_labels: List[str] = field(
        default_factory=lambda: ["essential-papers", "automated"]
    )

    def __post_init__(self):
        """Validate GitHub configuration."""
        # Allow empty token for local testing (will use mock mode)
        if self.token is None:
            self.token = "mock_token_for_local_testing"
        if not self.repository:
            raise ValueError("Repository cannot be empty")
        if "/" not in self.repository:
            raise ValueError("Repository must be in format 'owner/repo'")


@dataclass
class ScoringWeights:
    """Weights for paper scoring algorithm."""

    citation_weight: float = 0.4
    impact_factor_weight: float = 0.3
    recency_weight: float = 0.2
    relevance_weight: float = 0.1

    def __post_init__(self):
        """Validate scoring weights."""
        total = (
            self.citation_weight
            + self.impact_factor_weight
            + self.recency_weight
            + self.relevance_weight
        )
        if abs(total - 1.0) > 0.001:  # Allow small floating point errors
            raise ValueError(f"Weights must sum to 1.0, got {total}")

        weights = [
            self.citation_weight,
            self.impact_factor_weight,
            self.recency_weight,
            self.relevance_weight,
        ]
        if any(w < 0 for w in weights):
            raise ValueError("All weights must be non-negative")


@dataclass
class SystemConfig:
    """Overall system configuration."""

    topics: List[TopicConfig]
    github: GitHubConfig
    scoring_weights: ScoringWeights
    cache_settings: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate system configuration."""
        if not self.topics:
            raise ValueError("At least one topic must be configured")

        # Check for duplicate topic names
        topic_names = [topic.name for topic in self.topics]
        if len(topic_names) != len(set(topic_names)):
            raise ValueError("Topic names must be unique")
