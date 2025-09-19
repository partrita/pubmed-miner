"""
Unit tests for data models.
"""

import pytest
from datetime import datetime

from src.pubmed_miner.models import (
    Paper,
    ScoredPaper,
    TopicConfig,
    GitHubConfig,
    ScoringWeights,
    SystemConfig,
)


class TestPaper:
    """Test cases for Paper model."""

    def test_paper_creation_valid(self):
        """Test creating a valid paper."""
        paper = Paper(
            pmid="12345",
            title="Test Paper",
            authors=["John Doe", "Jane Smith"],
            journal="Nature",
            publication_date=datetime(2023, 1, 1),
        )

        assert paper.pmid == "12345"
        assert paper.title == "Test Paper"
        assert len(paper.authors) == 2
        assert paper.journal == "Nature"
        assert paper.abstract is None
        assert paper.doi is None

    def test_paper_creation_with_optional_fields(self):
        """Test creating a paper with optional fields."""
        paper = Paper(
            pmid="12345",
            title="Test Paper",
            authors=["John Doe"],
            journal="Nature",
            publication_date=datetime(2023, 1, 1),
            abstract="This is an abstract",
            doi="10.1038/nature12345",
        )

        assert paper.abstract == "This is an abstract"
        assert paper.doi == "10.1038/nature12345"

    def test_paper_validation_empty_pmid(self):
        """Test paper validation with empty PMID."""
        with pytest.raises(ValueError, match="PMID cannot be empty"):
            Paper(
                pmid="",
                title="Test Paper",
                authors=["John Doe"],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
            )

    def test_paper_validation_empty_title(self):
        """Test paper validation with empty title."""
        with pytest.raises(ValueError, match="Title cannot be empty"):
            Paper(
                pmid="12345",
                title="",
                authors=["John Doe"],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
            )

    def test_paper_validation_empty_authors(self):
        """Test paper validation with empty authors list."""
        with pytest.raises(ValueError, match="Authors list cannot be empty"):
            Paper(
                pmid="12345",
                title="Test Paper",
                authors=[],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
            )

    def test_paper_validation_empty_journal(self):
        """Test paper validation with empty journal."""
        with pytest.raises(ValueError, match="Journal cannot be empty"):
            Paper(
                pmid="12345",
                title="Test Paper",
                authors=["John Doe"],
                journal="",
                publication_date=datetime(2023, 1, 1),
            )


class TestScoredPaper:
    """Test cases for ScoredPaper model."""

    def test_scored_paper_creation_valid(self):
        """Test creating a valid scored paper."""
        scored_paper = ScoredPaper(
            pmid="12345",
            title="Test Paper",
            authors=["John Doe"],
            journal="Nature",
            publication_date=datetime(2023, 1, 1),
            citation_count=150,
            impact_factor=42.8,
            score=85.5,
            rank=1,
        )

        assert scored_paper.citation_count == 150
        assert scored_paper.impact_factor == 42.8
        assert scored_paper.score == 85.5
        assert scored_paper.rank == 1

    def test_scored_paper_validation_negative_citations(self):
        """Test scored paper validation with negative citations."""
        with pytest.raises(ValueError, match="Citation count cannot be negative"):
            ScoredPaper(
                pmid="12345",
                title="Test Paper",
                authors=["John Doe"],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
                citation_count=-10,
            )

    def test_scored_paper_validation_negative_impact_factor(self):
        """Test scored paper validation with negative impact factor."""
        with pytest.raises(ValueError, match="Impact factor cannot be negative"):
            ScoredPaper(
                pmid="12345",
                title="Test Paper",
                authors=["John Doe"],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
                impact_factor=-5.0,
            )

    def test_scored_paper_validation_negative_score(self):
        """Test scored paper validation with negative score."""
        with pytest.raises(ValueError, match="Score cannot be negative"):
            ScoredPaper(
                pmid="12345",
                title="Test Paper",
                authors=["John Doe"],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
                score=-10.0,
            )

    def test_scored_paper_validation_negative_rank(self):
        """Test scored paper validation with negative rank."""
        with pytest.raises(ValueError, match="Rank cannot be negative"):
            ScoredPaper(
                pmid="12345",
                title="Test Paper",
                authors=["John Doe"],
                journal="Nature",
                publication_date=datetime(2023, 1, 1),
                rank=-1,
            )


class TestTopicConfig:
    """Test cases for TopicConfig model."""

    def test_topic_config_creation_valid(self):
        """Test creating a valid topic configuration."""
        topic = TopicConfig(
            name="machine-learning",
            query="machine learning AND healthcare",
            max_papers=1000,
            essential_count=15,
            enabled=True,
        )

        assert topic.name == "machine-learning"
        assert topic.query == "machine learning AND healthcare"
        assert topic.max_papers == 1000
        assert topic.essential_count == 15
        assert topic.enabled is True

    def test_topic_config_defaults(self):
        """Test topic configuration with default values."""
        topic = TopicConfig(name="test-topic", query="test query")

        assert topic.max_papers == 1000
        assert topic.essential_count == 15
        assert topic.enabled is True

    def test_topic_config_validation_empty_name(self):
        """Test topic config validation with empty name."""
        with pytest.raises(ValueError, match="Topic name cannot be empty"):
            TopicConfig(name="", query="test query")

    def test_topic_config_validation_empty_query(self):
        """Test topic config validation with empty query."""
        with pytest.raises(ValueError, match="Query cannot be empty"):
            TopicConfig(name="test", query="")

    def test_topic_config_validation_invalid_max_papers(self):
        """Test topic config validation with invalid max_papers."""
        with pytest.raises(ValueError, match="Max papers must be positive"):
            TopicConfig(name="test", query="test", max_papers=0)

    def test_topic_config_validation_invalid_essential_count(self):
        """Test topic config validation with invalid essential_count."""
        with pytest.raises(ValueError, match="Essential count must be positive"):
            TopicConfig(name="test", query="test", essential_count=0)

    def test_topic_config_validation_essential_exceeds_max(self):
        """Test topic config validation when essential_count exceeds max_papers."""
        with pytest.raises(
            ValueError, match="Essential count cannot exceed max papers"
        ):
            TopicConfig(name="test", query="test", max_papers=10, essential_count=20)


class TestGitHubConfig:
    """Test cases for GitHubConfig model."""

    def test_github_config_creation_valid(self):
        """Test creating a valid GitHub configuration."""
        config = GitHubConfig(token="ghp_test_token", repository="owner/repo")

        assert config.token == "ghp_test_token"
        assert config.repository == "owner/repo"
        assert "essential-papers" in config.issue_labels
        assert "automated" in config.issue_labels

    def test_github_config_custom_labels(self):
        """Test GitHub config with custom labels."""
        config = GitHubConfig(
            token="ghp_test_token",
            repository="owner/repo",
            issue_labels=["research", "papers"],
        )

        assert config.issue_labels == ["research", "papers"]

    def test_github_config_validation_empty_token(self):
        """Test GitHub config validation with empty token."""
        with pytest.raises(ValueError, match="GitHub token cannot be empty"):
            GitHubConfig(token="", repository="owner/repo")

    def test_github_config_validation_empty_repository(self):
        """Test GitHub config validation with empty repository."""
        with pytest.raises(ValueError, match="Repository cannot be empty"):
            GitHubConfig(token="token", repository="")

    def test_github_config_validation_invalid_repository_format(self):
        """Test GitHub config validation with invalid repository format."""
        with pytest.raises(
            ValueError, match="Repository must be in format 'owner/repo'"
        ):
            GitHubConfig(token="token", repository="invalid-repo")


class TestScoringWeights:
    """Test cases for ScoringWeights model."""

    def test_scoring_weights_creation_valid(self):
        """Test creating valid scoring weights."""
        weights = ScoringWeights(
            citation_weight=0.4,
            impact_factor_weight=0.3,
            recency_weight=0.2,
            relevance_weight=0.1,
        )

        assert weights.citation_weight == 0.4
        assert weights.impact_factor_weight == 0.3
        assert weights.recency_weight == 0.2
        assert weights.relevance_weight == 0.1

    def test_scoring_weights_defaults(self):
        """Test scoring weights with default values."""
        weights = ScoringWeights()

        assert weights.citation_weight == 0.4
        assert weights.impact_factor_weight == 0.3
        assert weights.recency_weight == 0.2
        assert weights.relevance_weight == 0.1

    def test_scoring_weights_validation_sum_not_one(self):
        """Test scoring weights validation when weights don't sum to 1."""
        with pytest.raises(ValueError, match="Weights must sum to 1.0"):
            ScoringWeights(
                citation_weight=0.5,
                impact_factor_weight=0.5,
                recency_weight=0.2,
                relevance_weight=0.1,
            )

    def test_scoring_weights_validation_negative_weight(self):
        """Test scoring weights validation with negative weight."""
        with pytest.raises(ValueError, match="All weights must be non-negative"):
            ScoringWeights(
                citation_weight=-0.1,
                impact_factor_weight=0.4,
                recency_weight=0.4,
                relevance_weight=0.3,
            )

    def test_scoring_weights_validation_floating_point_tolerance(self):
        """Test scoring weights validation with floating point tolerance."""
        # This should not raise an error due to floating point tolerance
        weights = ScoringWeights(
            citation_weight=0.4,
            impact_factor_weight=0.3,
            recency_weight=0.2,
            relevance_weight=0.1000001,  # Slightly over due to floating point
        )

        assert weights is not None


class TestSystemConfig:
    """Test cases for SystemConfig model."""

    def test_system_config_creation_valid(self):
        """Test creating a valid system configuration."""
        topics = [
            TopicConfig(name="topic1", query="query1"),
            TopicConfig(name="topic2", query="query2"),
        ]
        github = GitHubConfig(token="token", repository="owner/repo")
        weights = ScoringWeights()

        config = SystemConfig(topics=topics, github=github, scoring_weights=weights)

        assert len(config.topics) == 2
        assert config.github.repository == "owner/repo"
        assert config.scoring_weights.citation_weight == 0.4

    def test_system_config_validation_empty_topics(self):
        """Test system config validation with empty topics list."""
        github = GitHubConfig(token="token", repository="owner/repo")
        weights = ScoringWeights()

        with pytest.raises(ValueError, match="At least one topic must be configured"):
            SystemConfig(topics=[], github=github, scoring_weights=weights)

    def test_system_config_validation_duplicate_topic_names(self):
        """Test system config validation with duplicate topic names."""
        topics = [
            TopicConfig(name="duplicate", query="query1"),
            TopicConfig(name="duplicate", query="query2"),
        ]
        github = GitHubConfig(token="token", repository="owner/repo")
        weights = ScoringWeights()

        with pytest.raises(ValueError, match="Topic names must be unique"):
            SystemConfig(topics=topics, github=github, scoring_weights=weights)
