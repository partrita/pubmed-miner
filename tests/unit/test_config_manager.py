"""
Unit tests for ConfigurationManager.
"""

import pytest
import tempfile
import yaml
from pathlib import Path
from unittest.mock import patch

from src.pubmed_miner.utils.config_manager import ConfigurationManager


class TestConfigurationManager:
    """Test cases for ConfigurationManager class."""

    def setup_method(self):
        """Set up test fixtures."""
        # Create temporary directory for test configs
        self.temp_dir = tempfile.mkdtemp()
        self.config_manager = ConfigurationManager(self.temp_dir)

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_initialization(self):
        """Test configuration manager initialization."""
        assert self.config_manager.config_dir == Path(self.temp_dir)
        assert self.config_manager.topics_file == Path(self.temp_dir) / "topics.yaml"
        assert (
            self.config_manager.settings_file == Path(self.temp_dir) / "settings.yaml"
        )

    def test_load_topics_success(self):
        """Test successful topic loading."""
        # Create test topics file
        topics_data = {
            "topics": [
                {
                    "name": "machine-learning",
                    "query": "machine learning AND healthcare",
                    "max_papers": 1000,
                    "essential_count": 15,
                    "enabled": True,
                },
                {
                    "name": "covid-research",
                    "query": "COVID-19 AND treatment",
                    "max_papers": 500,
                    "essential_count": 10,
                    "enabled": False,
                },
            ]
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        topics = self.config_manager.load_topics()

        assert len(topics) == 2
        assert topics[0].name == "machine-learning"
        assert topics[0].query == "machine learning AND healthcare"
        assert topics[0].enabled is True
        assert topics[1].name == "covid-research"
        assert topics[1].enabled is False

    def test_load_topics_file_not_found(self):
        """Test topic loading when file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Topics file not found"):
            self.config_manager.load_topics()

    def test_load_topics_invalid_yaml(self):
        """Test topic loading with invalid YAML."""
        # Create invalid YAML file
        with open(self.config_manager.topics_file, "w") as f:
            f.write("invalid: yaml: content: [")

        with pytest.raises(ValueError, match="Invalid YAML in topics file"):
            self.config_manager.load_topics()

    def test_load_topics_missing_topics_key(self):
        """Test topic loading with missing 'topics' key."""
        # Create YAML without 'topics' key
        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump({"other_key": "value"}, f)

        with pytest.raises(ValueError, match="Topics file must contain 'topics' key"):
            self.config_manager.load_topics()

    def test_load_topics_invalid_topic_data(self):
        """Test topic loading with invalid topic data."""
        topics_data = {
            "topics": [
                {
                    "name": "test-topic",
                    # Missing required 'query' field
                    "max_papers": 1000,
                }
            ]
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        with pytest.raises(ValueError, match="Invalid topic configuration"):
            self.config_manager.load_topics()

    @patch.dict("os.environ", {"GITHUB_TOKEN": "test_token_123"})
    def test_get_github_settings_success(self):
        """Test successful GitHub settings loading."""
        settings_data = {
            "github": {
                "repository": "testuser/testrepo",
                "issue_labels": ["papers", "automated"],
            }
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        github_config = self.config_manager.get_github_settings()

        assert github_config.token == "test_token_123"
        assert github_config.repository == "testuser/testrepo"
        assert github_config.issue_labels == ["papers", "automated"]

    def test_get_github_settings_no_token(self):
        """Test GitHub settings when no token is available (should use mock token)."""
        with patch.dict("os.environ", {}, clear=True):
            github_config = self.config_manager.get_github_settings()
            assert github_config.token == "mock_token_for_local_testing"
            assert github_config.repository == "local/test-repo"

    def test_get_scoring_weights_success(self):
        """Test successful scoring weights loading."""
        settings_data = {
            "scoring_weights": {
                "citation_weight": 0.5,
                "impact_factor_weight": 0.3,
                "recency_weight": 0.1,
                "relevance_weight": 0.1,
            }
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        weights = self.config_manager.get_scoring_weights()

        assert weights.citation_weight == 0.5
        assert weights.impact_factor_weight == 0.3
        assert weights.recency_weight == 0.1
        assert weights.relevance_weight == 0.1

    def test_get_scoring_weights_defaults(self):
        """Test scoring weights with default values."""
        # Create empty settings file
        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump({}, f)

        weights = self.config_manager.get_scoring_weights()

        assert weights.citation_weight == 0.4
        assert weights.impact_factor_weight == 0.3
        assert weights.recency_weight == 0.2
        assert weights.relevance_weight == 0.1

    @patch.dict("os.environ", {"GITHUB_TOKEN": "test_token"})
    def test_load_system_config_success(self):
        """Test successful system configuration loading."""
        # Create topics file
        topics_data = {
            "topics": [
                {
                    "name": "test-topic",
                    "query": "test query",
                    "max_papers": 100,
                    "essential_count": 5,
                    "enabled": True,
                }
            ]
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        # Create settings file
        settings_data = {
            "github": {"repository": "test/repo", "issue_labels": ["test"]},
            "scoring_weights": {
                "citation_weight": 0.4,
                "impact_factor_weight": 0.3,
                "recency_weight": 0.2,
                "relevance_weight": 0.1,
            },
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        system_config = self.config_manager.load_system_config()

        assert len(system_config.topics) == 1
        assert system_config.github.repository == "test/repo"
        assert system_config.scoring_weights.citation_weight == 0.4

    @patch.dict("os.environ", {"GITHUB_TOKEN": "test_token"})
    def test_validate_config_success(self):
        """Test successful configuration validation."""
        # Create valid configuration files
        topics_data = {
            "topics": [
                {
                    "name": "test-topic",
                    "query": "test query",
                    "max_papers": 100,
                    "essential_count": 5,
                    "enabled": True,
                }
            ]
        }

        settings_data = {
            "github": {"repository": "test/repo", "issue_labels": ["test"]}
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        # Should not raise any exceptions
        assert self.config_manager.validate_config() is True

    def test_validate_config_no_topics(self):
        """Test configuration validation with no topics."""
        # Create empty topics file
        topics_data = {"topics": []}

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        with pytest.raises(ValueError, match="No topics configured"):
            self.config_manager.validate_config()

    def test_create_default_configs(self):
        """Test creation of default configuration files."""
        self.config_manager.create_default_configs()

        # Check that files were created
        assert self.config_manager.topics_file.exists()
        assert self.config_manager.settings_file.exists()

        # Check topics file content
        with open(self.config_manager.topics_file, "r") as f:
            topics_data = yaml.safe_load(f)

        assert "topics" in topics_data
        assert len(topics_data["topics"]) >= 2
        assert topics_data["topics"][0]["name"] == "machine-learning-healthcare"

        # Check settings file content
        with open(self.config_manager.settings_file, "r") as f:
            settings_data = yaml.safe_load(f)

        assert "github" in settings_data
        assert "scoring_weights" in settings_data

    def test_create_default_configs_existing_files(self):
        """Test that existing files are not overwritten."""
        # Create existing files
        existing_topics = {"topics": [{"name": "existing", "query": "existing"}]}
        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(existing_topics, f)

        existing_settings = {"github": {"repository": "existing/repo"}}
        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(existing_settings, f)

        # Call create_default_configs
        self.config_manager.create_default_configs()

        # Check that existing content is preserved
        with open(self.config_manager.topics_file, "r") as f:
            topics_data = yaml.safe_load(f)

        assert topics_data["topics"][0]["name"] == "existing"

    def test_update_topic_success(self):
        """Test successful topic update."""
        # Create initial topics
        topics_data = {
            "topics": [
                {
                    "name": "test-topic",
                    "query": "original query",
                    "max_papers": 100,
                    "essential_count": 5,
                    "enabled": True,
                }
            ]
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        # Update topic
        self.config_manager.update_topic(
            "test-topic", query="updated query", max_papers=200, enabled=False
        )

        # Verify update
        topics = self.config_manager.load_topics()
        assert len(topics) == 1
        assert topics[0].query == "updated query"
        assert topics[0].max_papers == 200
        assert topics[0].enabled is False
        assert topics[0].essential_count == 5  # Unchanged

    def test_update_topic_not_found(self):
        """Test updating a non-existent topic."""
        topics_data = {
            "topics": [
                {
                    "name": "existing-topic",
                    "query": "test query",
                    "max_papers": 100,
                    "essential_count": 5,
                    "enabled": True,
                }
            ]
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        with pytest.raises(ValueError, match="Topic 'nonexistent' not found"):
            self.config_manager.update_topic("nonexistent", query="new query")

    def test_load_cache_settings(self):
        """Test loading cache settings."""
        settings_data = {
            "cache_settings": {
                "citation_cache_days": 14,
                "impact_factor_cache_days": 180,
                "paper_metadata_cache_days": 60,
            }
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        cache_settings = self.config_manager._load_cache_settings()

        assert cache_settings["citation_cache_days"] == 14
        assert cache_settings["impact_factor_cache_days"] == 180
        assert cache_settings["paper_metadata_cache_days"] == 60

    def test_load_cache_settings_defaults(self):
        """Test loading cache settings with defaults."""
        # Create empty settings file
        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump({}, f)

        cache_settings = self.config_manager._load_cache_settings()

        assert cache_settings["citation_cache_days"] == 7
        assert cache_settings["impact_factor_cache_days"] == 365
        assert cache_settings["paper_metadata_cache_days"] == 30

    def test_load_settings_file_not_found(self):
        """Test loading settings when file doesn't exist."""
        settings = self.config_manager._load_settings()
        assert settings == {}

    def test_load_settings_invalid_yaml(self):
        """Test loading settings with invalid YAML."""
        with open(self.config_manager.settings_file, "w") as f:
            f.write("invalid: yaml: [")

        settings = self.config_manager._load_settings()
        assert settings == {}

    def test_github_token_from_settings(self):
        """Test GitHub token loading from settings file."""
        settings_data = {
            "github": {"token": "settings_token_123", "repository": "test/repo"}
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        # Ensure no environment token
        with patch.dict("os.environ", {}, clear=True):
            github_config = self.config_manager.get_github_settings()
            assert github_config.token == "settings_token_123"

    def test_environment_token_priority(self):
        """Test that environment token takes priority over settings."""
        settings_data = {
            "github": {"token": "settings_token", "repository": "test/repo"}
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        with patch.dict("os.environ", {"GITHUB_TOKEN": "env_token"}):
            github_config = self.config_manager.get_github_settings()
            assert github_config.token == "env_token"

    def test_concurrent_access(self):
        """Test concurrent access to configuration files."""
        import threading

        # Create initial config
        topics_data = {
            "topics": [
                {
                    "name": "test-topic",
                    "query": "test query",
                    "max_papers": 100,
                    "essential_count": 5,
                    "enabled": True,
                }
            ]
        }

        with open(self.config_manager.topics_file, "w") as f:
            yaml.dump(topics_data, f)

        results = []
        errors = []

        def read_config():
            try:
                topics = self.config_manager.load_topics()
                results.append(len(topics))
            except Exception as e:
                errors.append(e)

        # Create multiple threads
        threads = []
        for i in range(10):
            thread = threading.Thread(target=read_config)
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # All reads should succeed
        assert len(errors) == 0
        assert all(result == 1 for result in results)

    @patch.dict("os.environ", {"PUBMED_EMAIL": "env@example.com"})
    def test_get_pubmed_email_from_environment(self):
        """Test PubMed email loading from environment variable."""
        email = self.config_manager.get_pubmed_email()
        assert email == "env@example.com"

    def test_get_pubmed_email_from_settings(self):
        """Test PubMed email loading from settings file."""
        settings_data = {"pubmed": {"email": "config@example.com"}}

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        # Ensure no environment email
        with patch.dict("os.environ", {}, clear=True):
            email = self.config_manager.get_pubmed_email()
            assert email == "config@example.com"

    def test_get_pubmed_email_environment_priority(self):
        """Test that environment email takes priority over settings."""
        settings_data = {"pubmed": {"email": "config@example.com"}}

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        with patch.dict("os.environ", {"PUBMED_EMAIL": "env@example.com"}):
            email = self.config_manager.get_pubmed_email()
            assert email == "env@example.com"

    def test_get_pubmed_email_default_fallback(self):
        """Test PubMed email fallback to default."""
        # Create empty settings file
        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump({}, f)

        with patch.dict("os.environ", {}, clear=True):
            email = self.config_manager.get_pubmed_email()
            assert email == "pubmed.miner@example.com"

    def test_get_pubmed_email_ignore_default_in_config(self):
        """Test that default email in config is ignored."""
        settings_data = {
            "pubmed": {
                "email": "pubmed.miner@example.com"  # This is the default, should be ignored
            }
        }

        with open(self.config_manager.settings_file, "w") as f:
            yaml.dump(settings_data, f)

        with patch.dict("os.environ", {}, clear=True):
            email = self.config_manager.get_pubmed_email()
            assert email == "pubmed.miner@example.com"  # Falls back to default
