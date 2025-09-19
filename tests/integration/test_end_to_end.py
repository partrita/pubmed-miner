"""
End-to-end integration tests for the complete system.
"""

import pytest
import tempfile
import yaml
import os
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime

# Import the main automation script
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from automated_collection import AutomatedCollectionOrchestrator
from src.pubmed_miner.models import Paper


class TestEndToEndWorkflow:
    """End-to-end tests for the complete automated collection workflow."""

    def setup_method(self):
        """Set up test fixtures."""
        # Create temporary directory for test configs
        self.temp_dir = tempfile.mkdtemp()

        # Create test configuration files
        self.create_test_configs()

        # Set environment variables for testing
        self.original_env = os.environ.copy()
        os.environ["GITHUB_TOKEN"] = "test_token_12345"
        os.environ["PUBMED_EMAIL"] = "test@example.com"

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

        # Restore original environment
        os.environ.clear()
        os.environ.update(self.original_env)

    def create_test_configs(self):
        """Create test configuration files."""
        # Create config directory
        config_dir = Path(self.temp_dir) / "config"
        config_dir.mkdir(exist_ok=True)

        # Create topics configuration
        topics_data = {
            "topics": [
                {
                    "name": "machine-learning-healthcare",
                    "query": "machine learning AND healthcare",
                    "max_papers": 10,
                    "essential_count": 3,
                    "enabled": True,
                },
                {
                    "name": "covid-research",
                    "query": "COVID-19 AND treatment",
                    "max_papers": 5,
                    "essential_count": 2,
                    "enabled": True,
                },
                {
                    "name": "disabled-topic",
                    "query": "disabled topic",
                    "max_papers": 10,
                    "essential_count": 5,
                    "enabled": False,
                },
            ]
        }

        topics_file = config_dir / "topics.yaml"
        with open(topics_file, "w") as f:
            yaml.dump(topics_data, f)

        # Create settings configuration
        settings_data = {
            "github": {
                "repository": "testuser/testrepo",
                "issue_labels": ["essential-papers", "automated", "test"],
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
                "pubmed_rate_limit": 3,
                "github_rate_limit": 5000,
                "retry_attempts": 3,
                "retry_delay": 1.0,
            },
        }

        settings_file = config_dir / "settings.yaml"
        with open(settings_file, "w") as f:
            yaml.dump(settings_data, f)

        # Update working directory for the test
        os.chdir(self.temp_dir)

    def create_sample_papers(self, topic_name: str, count: int = 5):
        """Create sample papers for testing."""
        papers = []
        for i in range(count):
            paper = Paper(
                pmid=f"{topic_name}_{i + 1}",
                title=f"Sample Paper {i + 1} for {topic_name}",
                authors=[f"Author {i + 1}A", f"Author {i + 1}B"],
                journal=f"Journal {(i % 3) + 1}",
                publication_date=datetime(2023, (i % 12) + 1, 1),
                abstract=f"This is the abstract for paper {i + 1} about {topic_name}.",
                doi=f"10.1000/{topic_name}.{i + 1}",
            )
            papers.append(paper)
        return papers

    def test_complete_end_to_end_workflow(self):
        """Test the complete end-to-end workflow."""
        # Mock all external dependencies
        with patch("automated_collection.PaperCollectionService") as mock_paper_service:
            with patch("automated_collection.CitationService") as mock_citation_service:
                with patch(
                    "automated_collection.ImpactFactorService"
                ) as mock_if_service:
                    with patch(
                        "automated_collection.GitHubIssuesManager"
                    ) as mock_github_manager:
                        # Set up mock paper collection service
                        mock_paper_instance = Mock()
                        mock_paper_service.return_value = mock_paper_instance

                        # Mock paper search and details for each topic
                        def mock_search_papers(query, max_papers):
                            if "machine learning" in query:
                                return ["ml_1", "ml_2", "ml_3"]
                            elif "COVID-19" in query:
                                return ["covid_1", "covid_2"]
                            return []

                        def mock_get_paper_details(pmids):
                            papers = []
                            for pmid in pmids:
                                if pmid.startswith("ml_"):
                                    papers.extend(
                                        self.create_sample_papers("machine-learning", 1)
                                    )
                                elif pmid.startswith("covid_"):
                                    papers.extend(self.create_sample_papers("covid", 1))
                            return papers[: len(pmids)]

                        mock_paper_instance.search_papers.side_effect = (
                            mock_search_papers
                        )
                        mock_paper_instance.get_papers_details.side_effect = (
                            mock_get_paper_details
                        )

                        # Set up mock citation service
                        mock_citation_instance = Mock()
                        mock_citation_service.return_value = mock_citation_instance
                        mock_citation_instance.get_citation_count.side_effect = (
                            lambda pmid: hash(pmid) % 200
                        )

                        # Set up mock impact factor service
                        mock_if_instance = Mock()
                        mock_if_service.return_value = mock_if_instance
                        mock_if_instance.get_impact_factor.side_effect = (
                            lambda journal: hash(journal) % 50 + 1
                        )

                        # Set up mock GitHub manager
                        mock_github_instance = Mock()
                        mock_github_manager.return_value = mock_github_instance

                        def mock_create_or_update_issue(topic, papers):
                            return {
                                "number": hash(topic) % 100 + 1,
                                "created": True,
                                "html_url": f"https://github.com/testuser/testrepo/issues/{hash(topic) % 100 + 1}",
                            }

                        mock_github_instance.create_or_update_issue.side_effect = (
                            mock_create_or_update_issue
                        )

                        # Initialize and run orchestrator
                        orchestrator = AutomatedCollectionOrchestrator()
                        results = orchestrator.run_complete_workflow()

                        # Verify results
                        assert results["success"] is True
                        assert results["topics_processed"] == 2  # Only enabled topics
                        assert results["papers_collected"] > 0
                        assert results["issues_created"] == 2
                        assert len(results["errors"]) == 0

                        # Verify that disabled topics were not processed
                        search_calls = mock_paper_instance.search_papers.call_args_list
                        search_queries = [call[0][0] for call in search_calls]
                        assert not any(
                            "disabled topic" in query for query in search_queries
                        )

    def test_workflow_with_partial_failures(self):
        """Test workflow resilience with partial failures."""
        with patch("automated_collection.PaperCollectionService") as mock_paper_service:
            with patch("automated_collection.CitationService") as mock_citation_service:
                with patch(
                    "automated_collection.ImpactFactorService"
                ) as mock_if_service:
                    with patch(
                        "automated_collection.GitHubIssuesManager"
                    ) as mock_github_manager:
                        # Set up services with some failures
                        mock_paper_instance = Mock()
                        mock_paper_service.return_value = mock_paper_instance

                        def mock_search_with_failure(query, max_papers):
                            if "machine learning" in query:
                                return ["ml_1", "ml_2"]
                            elif "COVID-19" in query:
                                raise Exception("PubMed API error")
                            return []

                        mock_paper_instance.search_papers.side_effect = (
                            mock_search_with_failure
                        )
                        mock_paper_instance.get_papers_details.return_value = (
                            self.create_sample_papers("ml", 2)
                        )

                        # Citation service with intermittent failures
                        mock_citation_instance = Mock()
                        mock_citation_service.return_value = mock_citation_instance

                        def mock_citation_with_failure(pmid):
                            if "1" in pmid:
                                raise Exception("Citation API error")
                            return 100

                        mock_citation_instance.get_citation_count.side_effect = (
                            mock_citation_with_failure
                        )

                        # Impact factor service working
                        mock_if_instance = Mock()
                        mock_if_service.return_value = mock_if_instance
                        mock_if_instance.get_impact_factor.return_value = 10.0

                        # GitHub manager working
                        mock_github_instance = Mock()
                        mock_github_manager.return_value = mock_github_instance
                        mock_github_instance.create_or_update_issue.return_value = {
                            "number": 42,
                            "created": True,
                            "html_url": "https://github.com/testuser/testrepo/issues/42",
                        }

                        # Run orchestrator
                        orchestrator = AutomatedCollectionOrchestrator()
                        results = orchestrator.run_complete_workflow()

                        # Should still succeed partially
                        assert (
                            results["topics_processed"] == 1
                        )  # Only one topic succeeded
                        assert len(results["errors"]) == 1  # One topic failed
                        assert (
                            results["success"] is True
                        )  # Overall success due to partial completion

    def test_workflow_with_no_papers_found(self):
        """Test workflow when no papers are found for topics."""
        with patch("automated_collection.PaperCollectionService") as mock_paper_service:
            with patch(
                "automated_collection.GitHubIssuesManager"
            ) as mock_github_manager:
                # Set up paper service to return no results
                mock_paper_instance = Mock()
                mock_paper_service.return_value = mock_paper_instance
                mock_paper_instance.search_papers.return_value = []
                mock_paper_instance.get_papers_details.return_value = []

                # Set up GitHub manager
                mock_github_instance = Mock()
                mock_github_manager.return_value = mock_github_instance
                mock_github_instance.create_or_update_issue.return_value = {
                    "number": 43,
                    "created": True,
                    "html_url": "https://github.com/testuser/testrepo/issues/43",
                }

                # Run orchestrator
                orchestrator = AutomatedCollectionOrchestrator()
                results = orchestrator.run_complete_workflow()

                # Should still create issues with "no papers found" message
                assert results["success"] is True
                assert results["topics_processed"] == 2
                assert results["papers_collected"] == 0
                assert results["issues_created"] == 2

    def test_workflow_configuration_errors(self):
        """Test workflow behavior with configuration errors."""
        # Create invalid configuration
        config_dir = Path(self.temp_dir) / "config"

        # Overwrite with invalid topics file
        invalid_topics = {
            "topics": [
                {
                    "name": "",  # Invalid empty name
                    "query": "test query",
                    "max_papers": -10,  # Invalid negative value
                    "essential_count": 0,  # Invalid zero value
                    "enabled": True,
                }
            ]
        }

        topics_file = config_dir / "topics.yaml"
        with open(topics_file, "w") as f:
            yaml.dump(invalid_topics, f)

        # Should fail during initialization
        with pytest.raises(Exception):
            AutomatedCollectionOrchestrator()

    def test_workflow_missing_environment_variables(self):
        """Test workflow behavior with missing environment variables."""
        # Remove required environment variable
        del os.environ["GITHUB_TOKEN"]

        # Should fail during initialization
        with pytest.raises(Exception):
            AutomatedCollectionOrchestrator()

    def test_workflow_statistics_collection(self):
        """Test that workflow collects proper statistics."""
        with patch("automated_collection.PaperCollectionService") as mock_paper_service:
            with patch("automated_collection.CitationService") as mock_citation_service:
                with patch(
                    "automated_collection.ImpactFactorService"
                ) as mock_if_service:
                    with patch(
                        "automated_collection.GitHubIssuesManager"
                    ) as mock_github_manager:
                        # Set up mocks to return predictable data
                        mock_paper_instance = Mock()
                        mock_paper_service.return_value = mock_paper_instance
                        mock_paper_instance.search_papers.return_value = ["1", "2", "3"]
                        mock_paper_instance.get_papers_details.return_value = (
                            self.create_sample_papers("test", 3)
                        )

                        mock_citation_instance = Mock()
                        mock_citation_service.return_value = mock_citation_instance
                        mock_citation_instance.get_citation_count.return_value = 50

                        mock_if_instance = Mock()
                        mock_if_service.return_value = mock_if_instance
                        mock_if_instance.get_impact_factor.return_value = 5.0

                        mock_github_instance = Mock()
                        mock_github_manager.return_value = mock_github_instance
                        mock_github_instance.create_or_update_issue.return_value = {
                            "number": 1,
                            "created": True,
                            "html_url": "https://github.com/test/repo/issues/1",
                        }

                        # Run orchestrator
                        orchestrator = AutomatedCollectionOrchestrator()
                        results = orchestrator.run_complete_workflow()

                        # Verify statistics
                        assert "start_time" in results
                        assert "end_time" in results
                        assert "duration_seconds" in results
                        assert results["duration_seconds"] > 0
                        assert results["topics_processed"] == 2
                        assert (
                            results["papers_collected"] == 6
                        )  # 3 papers per topic, 2 topics
                        assert results["issues_created"] == 2

    def test_workflow_logging_integration(self):
        """Test that workflow properly integrates with logging system."""
        import logging
        from io import StringIO

        # Capture log output
        log_capture = StringIO()
        handler = logging.StreamHandler(log_capture)
        handler.setLevel(logging.INFO)

        logger = logging.getLogger()
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

        try:
            with patch(
                "automated_collection.PaperCollectionService"
            ) as mock_paper_service:
                with patch(
                    "automated_collection.GitHubIssuesManager"
                ) as mock_github_manager:
                    # Set up minimal mocks
                    mock_paper_instance = Mock()
                    mock_paper_service.return_value = mock_paper_instance
                    mock_paper_instance.search_papers.return_value = []
                    mock_paper_instance.get_papers_details.return_value = []

                    mock_github_instance = Mock()
                    mock_github_manager.return_value = mock_github_instance
                    mock_github_instance.create_or_update_issue.return_value = {
                        "number": 1,
                        "created": True,
                        "html_url": "https://github.com/test/repo/issues/1",
                    }

                    # Run orchestrator
                    orchestrator = AutomatedCollectionOrchestrator()
                    orchestrator.run_complete_workflow()

                    # Check that logs were generated
                    log_output = log_capture.getvalue()
                    assert (
                        "Starting automated essential papers collection workflow"
                        in log_output
                    )
                    assert "Workflow completed" in log_output

        finally:
            logger.removeHandler(handler)

    def test_workflow_error_recovery(self):
        """Test workflow error recovery and continuation."""
        with patch("automated_collection.PaperCollectionService") as mock_paper_service:
            with patch(
                "automated_collection.GitHubIssuesManager"
            ) as mock_github_manager:
                # Set up paper service with mixed success/failure
                mock_paper_instance = Mock()
                mock_paper_service.return_value = mock_paper_instance

                call_count = 0

                def mock_search_mixed_results(query, max_papers):
                    nonlocal call_count
                    call_count += 1
                    if call_count == 1:
                        # First topic succeeds
                        return ["1", "2"]
                    else:
                        # Second topic fails
                        raise Exception("API failure")

                mock_paper_instance.search_papers.side_effect = (
                    mock_search_mixed_results
                )
                mock_paper_instance.get_papers_details.return_value = (
                    self.create_sample_papers("test", 2)
                )

                # Set up GitHub manager
                mock_github_instance = Mock()
                mock_github_manager.return_value = mock_github_instance
                mock_github_instance.create_or_update_issue.return_value = {
                    "number": 1,
                    "created": True,
                    "html_url": "https://github.com/test/repo/issues/1",
                }

                # Run orchestrator
                orchestrator = AutomatedCollectionOrchestrator()
                results = orchestrator.run_complete_workflow()

                # Should continue processing despite one failure
                assert results["topics_processed"] == 1  # One succeeded
                assert len(results["errors"]) == 1  # One failed
                assert results["success"] is True  # Overall success
                assert (
                    results["issues_created"] == 1
                )  # Issue created for successful topic

    def test_workflow_performance_under_load(self):
        """Test workflow performance with larger datasets."""
        with patch("automated_collection.PaperCollectionService") as mock_paper_service:
            with patch("automated_collection.CitationService") as mock_citation_service:
                with patch(
                    "automated_collection.ImpactFactorService"
                ) as mock_if_service:
                    with patch(
                        "automated_collection.GitHubIssuesManager"
                    ) as mock_github_manager:
                        # Set up services to return larger datasets
                        mock_paper_instance = Mock()
                        mock_paper_service.return_value = mock_paper_instance

                        # Return many PMIDs
                        mock_paper_instance.search_papers.return_value = [
                            str(i) for i in range(100)
                        ]
                        mock_paper_instance.get_papers_details.return_value = (
                            self.create_sample_papers("perf", 100)
                        )

                        # Fast citation service
                        mock_citation_instance = Mock()
                        mock_citation_service.return_value = mock_citation_instance
                        mock_citation_instance.get_citation_count.return_value = 50

                        # Fast impact factor service
                        mock_if_instance = Mock()
                        mock_if_service.return_value = mock_if_instance
                        mock_if_instance.get_impact_factor.return_value = 5.0

                        # Fast GitHub manager
                        mock_github_instance = Mock()
                        mock_github_manager.return_value = mock_github_instance
                        mock_github_instance.create_or_update_issue.return_value = {
                            "number": 1,
                            "created": True,
                            "html_url": "https://github.com/test/repo/issues/1",
                        }

                        # Measure performance
                        import time

                        start_time = time.time()

                        orchestrator = AutomatedCollectionOrchestrator()
                        results = orchestrator.run_complete_workflow()

                        end_time = time.time()
                        execution_time = end_time - start_time

                        # Should complete within reasonable time
                        assert execution_time < 30.0  # 30 seconds max
                        assert results["success"] is True
                        assert (
                            results["papers_collected"] == 200
                        )  # 100 papers per topic, 2 topics
