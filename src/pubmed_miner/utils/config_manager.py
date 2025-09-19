"""
Configuration management utilities.
"""

import yaml
from pathlib import Path
from typing import List, Dict, Any
import logging

from ..models import TopicConfig, GitHubConfig, ScoringWeights, SystemConfig

logger = logging.getLogger(__name__)


class ConfigurationManager:
    """Manages system configuration loading and validation."""

    def __init__(self, config_dir: str = "config"):
        """Initialize configuration manager.

        Args:
            config_dir: Directory containing configuration files
        """
        self.config_dir = Path(config_dir)
        self.topics_file = self.config_dir / "topics.yaml"
        self.settings_file = self.config_dir / "settings.yaml"

    def load_topics(self) -> List[TopicConfig]:
        """Load topic configurations from YAML file.

        Returns:
            List of TopicConfig objects

        Raises:
            FileNotFoundError: If topics file doesn't exist
            ValueError: If configuration is invalid
        """
        if not self.topics_file.exists():
            raise FileNotFoundError(f"Topics file not found: {self.topics_file}")

        try:
            with open(self.topics_file, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)

            if not data or "topics" not in data:
                raise ValueError("Topics file must contain 'topics' key")

            topics = []
            for topic_data in data["topics"]:
                topic = TopicConfig(**topic_data)
                topics.append(topic)

            logger.info(f"Loaded {len(topics)} topics from {self.topics_file}")
            return topics

        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in topics file: {e}")
        except TypeError as e:
            raise ValueError(f"Invalid topic configuration: {e}")

    def load_system_config(self) -> SystemConfig:
        """Load complete system configuration.

        Returns:
            SystemConfig object with all settings
        """
        topics = self.load_topics()
        github_config = self.get_github_settings()
        scoring_weights = self.get_scoring_weights()
        cache_settings = self._load_cache_settings()

        return SystemConfig(
            topics=topics,
            github=github_config,
            scoring_weights=scoring_weights,
            cache_settings=cache_settings,
        )

    def get_github_settings(self) -> GitHubConfig:
        """Load GitHub configuration from settings file.

        Returns:
            GitHubConfig object
        """
        settings = self._load_settings()
        github_data = settings.get("github", {})

        # Get token from environment or settings
        import os

        token = os.getenv("GITHUB_TOKEN") or github_data.get("token", "")

        # For local testing, allow empty token (will use mock mode)
        if not token:
            logger.warning(
                "GitHub token not found - GitHub features will be disabled for local testing"
            )
            token = "mock_token_for_local_testing"

        return GitHubConfig(
            token=token,
            repository=github_data.get("repository", "local/test-repo"),
            issue_labels=github_data.get(
                "issue_labels", ["essential-papers", "automated"]
            ),
        )

    def get_scoring_weights(self) -> ScoringWeights:
        """Load scoring weights from settings file.

        Returns:
            ScoringWeights object with configured or default weights
        """
        settings = self._load_settings()
        weights_data = settings.get("scoring_weights", {})

        return ScoringWeights(
            citation_weight=weights_data.get("citation_weight", 0.4),
            impact_factor_weight=weights_data.get("impact_factor_weight", 0.3),
            recency_weight=weights_data.get("recency_weight", 0.2),
            relevance_weight=weights_data.get("relevance_weight", 0.1),
        )

    def get_pubmed_email(self) -> str:
        """Load PubMed email from settings file or environment.

        Returns:
            Email address for PubMed API
        """
        import os

        # First check environment variable
        env_email = os.getenv("PUBMED_EMAIL")
        if env_email:
            return env_email

        # Then check settings file
        settings = self._load_settings()
        pubmed_data = settings.get("pubmed", {})
        config_email = pubmed_data.get("email")

        if config_email and config_email != "pubmed.miner@example.com":
            return config_email

        # Default fallback
        logger.warning("PubMed email not configured - using default")
        return "pubmed.miner@example.com"

    def validate_config(self) -> bool:
        """Validate all configuration files.

        Returns:
            True if all configurations are valid

        Raises:
            ValueError: If any configuration is invalid
        """
        try:
            # Validate topics
            topics = self.load_topics()
            if not topics:
                raise ValueError("No topics configured")

            # Validate GitHub settings
            github_config = self.get_github_settings()

            # Validate scoring weights
            scoring_weights = self.get_scoring_weights()

            # Validate system config as a whole
            SystemConfig(
                topics=topics, github=github_config, scoring_weights=scoring_weights
            )

            logger.info("Configuration validation successful")
            return True

        except Exception as e:
            logger.error(f"Configuration validation failed: {e}")
            raise

    def create_default_configs(self) -> None:
        """Create default configuration files if they don't exist."""
        self.config_dir.mkdir(exist_ok=True)

        # Create default topics.yaml
        if not self.topics_file.exists():
            default_topics = {
                "topics": [
                    {
                        "name": "machine-learning-healthcare",
                        "query": "machine learning AND healthcare",
                        "max_papers": 1000,
                        "essential_count": 15,
                        "enabled": True,
                    },
                    {
                        "name": "covid-19-treatment",
                        "query": "COVID-19 AND treatment",
                        "max_papers": 500,
                        "essential_count": 10,
                        "enabled": True,
                    },
                ]
            }

            with open(self.topics_file, "w", encoding="utf-8") as f:
                yaml.dump(
                    default_topics, f, default_flow_style=False, allow_unicode=True
                )
            logger.info(f"Created default topics file: {self.topics_file}")

        # Create default settings.yaml
        if not self.settings_file.exists():
            default_settings = {
                "github": {
                    "repository": "username/repo-name",
                    "issue_labels": ["essential-papers", "automated"],
                },
                "scoring_weights": {
                    "citation_weight": 0.4,
                    "impact_factor_weight": 0.3,
                    "recency_weight": 0.2,
                    "relevance_weight": 0.1,
                },
                "cache_settings": {
                    "citation_cache_days": 7,
                    "impact_factor_cache_days": 365,
                    "paper_metadata_cache_days": 30,
                },
                "api_settings": {
                    "pubmed_rate_limit": 3,  # requests per second
                    "github_rate_limit": 5000,  # requests per hour
                    "retry_attempts": 3,
                    "retry_delay": 1.0,
                },
            }

            with open(self.settings_file, "w", encoding="utf-8") as f:
                yaml.dump(
                    default_settings, f, default_flow_style=False, allow_unicode=True
                )
            logger.info(f"Created default settings file: {self.settings_file}")

    def _load_settings(self) -> Dict[str, Any]:
        """Load settings from YAML file.

        Returns:
            Dictionary with settings data
        """
        if not self.settings_file.exists():
            logger.warning(
                f"Settings file not found: {self.settings_file}, using defaults"
            )
            return {}

        try:
            with open(self.settings_file, "r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except yaml.YAMLError as e:
            logger.error(f"Error loading settings file: {e}")
            return {}

    def _load_cache_settings(self) -> Dict[str, Any]:
        """Load cache-related settings.

        Returns:
            Dictionary with cache settings
        """
        settings = self._load_settings()
        return settings.get(
            "cache_settings",
            {
                "citation_cache_days": 7,
                "impact_factor_cache_days": 365,
                "paper_metadata_cache_days": 30,
            },
        )

    def update_topic(self, topic_name: str, **kwargs) -> None:
        """Update a specific topic configuration.

        Args:
            topic_name: Name of the topic to update
            **kwargs: Fields to update
        """
        topics = self.load_topics()

        # Find and update the topic
        topic_found = False
        for i, topic in enumerate(topics):
            if topic.name == topic_name:
                # Create updated topic
                topic_dict = {
                    "name": topic.name,
                    "query": topic.query,
                    "max_papers": topic.max_papers,
                    "essential_count": topic.essential_count,
                    "enabled": topic.enabled,
                }
                topic_dict.update(kwargs)
                topics[i] = TopicConfig(**topic_dict)
                topic_found = True
                break

        if not topic_found:
            raise ValueError(f"Topic '{topic_name}' not found")

        # Save updated topics
        topics_data = {
            "topics": [
                {
                    "name": t.name,
                    "query": t.query,
                    "max_papers": t.max_papers,
                    "essential_count": t.essential_count,
                    "enabled": t.enabled,
                }
                for t in topics
            ]
        }

        with open(self.topics_file, "w", encoding="utf-8") as f:
            yaml.dump(topics_data, f, default_flow_style=False, allow_unicode=True)

        logger.info(f"Updated topic '{topic_name}'")
