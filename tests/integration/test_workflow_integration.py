"""
Integration tests for the complete workflow.
"""

import pytest
import tempfile
import yaml
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime

from src.pubmed_miner.utils.config_manager import ConfigurationManager
from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.services.citation_service import CitationService
from src.pubmed_miner.services.impact_factor_service import ImpactFactorService
from src.pubmed_miner.scoring.engine import ScoringEngine
from src.pubmed_miner.services.github_manager import GitHubIssuesManager
from src.pubmed_miner.models import Paper, ScoredPaper, TopicConfig, GitHubConfig


class TestWorkflowIntegration:
    """Integration tests for the complete essential papers workflow."""

    def setup_method(self):
        """Set up test fixtures."""
        # Create temporary directory for test configs
        self.temp_dir = tempfile.mkdtemp()

        # Create test configuration files
        self.create_test_configs()

        # Initialize services
        self.config_manager = ConfigurationManager(self.temp_dir)
        self.paper_collector = PaperCollectionService(email="test@example.com")
        self.citation_service = CitationService()
        self.impact_factor_service = ImpactFactorService()
        self.scoring_engine = ScoringEngine()

        # Mock GitHub config for testing
        github_config = GitHubConfig(
            token="test_token", repository="test/repo", issue_labels=["test"]
        )
        self.github_manager = GitHubIssuesManager(github_config)

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def create_test_configs(self):
        """Create test configuration files."""
        # Create topics configuration
        topics_data = {
            "topics": [
                {
                    "name": "machine-learning-test",
                    "query": "machine learning AND healthcare",
                    "max_papers": 10,
                    "essential_count": 3,
                    "enabled": True,
                }
            ]
        }

        topics_file = Path(self.temp_dir) / "topics.yaml"
        with open(topics_file, "w") as f:
            yaml.dump(topics_data, f)

        # Create settings configuration
        settings_data = {
            "github": {
                "repository": "test/repo",
                "issue_labels": ["test", "integration"],
            },
            "scoring_weights": {
                "citation_weight": 0.4,
                "impact_factor_weight": 0.3,
                "recency_weight": 0.2,
                "relevance_weight": 0.1,
            },
        }

        settings_file = Path(self.temp_dir) / "settings.yaml"
        with open(settings_file, "w") as f:
            yaml.dump(settings_data, f)

    @patch.dict("os.environ", {"GITHUB_TOKEN": "test_token_123"})
    def test_complete_workflow_success(self):
        """Test the complete workflow from configuration to GitHub issue creation."""

        # Step 1: Load configuration
        system_config = self.config_manager.load_system_config()
        assert len(system_config.topics) == 1
        topic = system_config.topics[0]

        # Step 2: Mock paper collection
        sample_papers = [
            Paper(
                pmid="12345",
                title="Machine Learning in Healthcare: A Review",
                authors=["John Doe", "Jane Smith"],
                journal="Nature Medicine",
                publication_date=datetime(2023, 1, 15),
                abstract="This paper reviews ML applications in healthcare.",
                doi="10.1038/nm.2023.12345",
            ),
            Paper(
                pmid="67890",
                title="AI for Medical Diagnosis",
                authors=["Bob Johnson"],
                journal="Science",
                publication_date=datetime(2023, 2, 10),
                abstract="AI applications in medical diagnosis.",
                doi="10.1126/science.2023.67890",
            ),
            Paper(
                pmid="11111",
                title="Deep Learning in Radiology",
                authors=["Alice Brown"],
                journal="Radiology",
                publication_date=datetime(2023, 3, 5),
                abstract="Deep learning applications in radiology.",
                doi="10.1148/radiol.2023.11111",
            ),
        ]

        with patch.object(self.paper_collector, "search_papers") as mock_search:
            mock_search.return_value = ["12345", "67890", "11111"]

            with patch.object(
                self.paper_collector, "get_paper_details"
            ) as mock_details:
                mock_details.return_value = sample_papers

                # Step 3: Mock citation and impact factor services
                citation_counts = {"12345": 150, "67890": 75, "11111": 200}
                impact_factors = {
                    "Nature Medicine": 87.2,
                    "Science": 56.9,
                    "Radiology": 29.1,
                }

                with patch.object(
                    self.citation_service, "batch_get_citation_counts"
                ) as mock_citations:
                    mock_citations.return_value = citation_counts

                    with patch.object(
                        self.impact_factor_service, "get_impact_factor"
                    ) as mock_if:
                        mock_if.side_effect = lambda journal: impact_factors.get(
                            journal, 0.0
                        )

                        # Step 4: Score and rank papers
                        scored_papers = []
                        for paper in sample_papers:
                            citations = citation_counts.get(paper.pmid, 0)
                            impact_factor = impact_factors.get(paper.journal, 0.0)

                            score = self.scoring_engine.calculate_paper_score(
                                paper, citations, impact_factor, topic.query
                            )

                            scored_paper = ScoredPaper(
                                pmid=paper.pmid,
                                title=paper.title,
                                authors=paper.authors,
                                journal=paper.journal,
                                publication_date=paper.publication_date,
                                abstract=paper.abstract,
                                doi=paper.doi,
                                citation_count=citations,
                                impact_factor=impact_factor,
                                score=score,
                                rank=0,
                            )
                            scored_papers.append(scored_paper)

                        # Rank papers
                        ranked_papers = self.scoring_engine.rank_papers(scored_papers)

                        # Select essential papers
                        essential_papers = self.scoring_engine.select_essential_papers(
                            ranked_papers, topic.essential_count
                        )

                        assert (
                            len(essential_papers) == 3
                        )  # All papers since we have only 3
                        assert essential_papers[0].rank == 1
                        assert essential_papers[1].rank == 2
                        assert essential_papers[2].rank == 3

                        # Step 5: Mock GitHub issue creation
                        with patch.object(
                            self.github_manager, "create_or_update_issue"
                        ) as mock_github:
                            mock_github.return_value = {
                                "number": 42,
                                "created": True,
                                "html_url": "https://github.com/test/repo/issues/42",
                            }

                            # Create GitHub issue
                            issue_result = self.github_manager.create_or_update_issue(
                                topic.name, essential_papers
                            )

                            assert issue_result["created"] is True
                            assert issue_result["number"] == 42

                            # Verify GitHub manager was called with correct parameters
                            mock_github.assert_called_once_with(
                                topic.name, essential_papers
                            )

    def test_workflow_with_api_failures(self):
        """Test workflow resilience with API failures."""
        topic = TopicConfig(
            name="test-topic",
            query="test query",
            max_papers=5,
            essential_count=2,
            enabled=True,
        )

        # Mock paper collection failure and retry
        with patch.object(self.paper_collector, "search_papers") as mock_search:
            # First call fails, second succeeds
            mock_search.side_effect = [
                Exception("PubMed API error"),
                ["12345", "67890"],
            ]

            sample_papers = [
                Paper(
                    pmid="12345",
                    title="Test Paper 1",
                    authors=["Author 1"],
                    journal="Test Journal",
                    publication_date=datetime(2023, 1, 1),
                    doi="10.1000/test.12345",
                ),
                Paper(
                    pmid="67890",
                    title="Test Paper 2",
                    authors=["Author 2"],
                    journal="Test Journal 2",
                    publication_date=datetime(2023, 2, 1),
                    doi="10.1000/test.67890",
                ),
            ]

            with patch.object(
                self.paper_collector, "get_paper_details"
            ) as mock_details:
                mock_details.return_value = sample_papers

                # Mock citation service with partial failures
                with patch.object(
                    self.citation_service, "get_citation_count"
                ) as mock_citations:

                    def citation_side_effect(pmid, **kwargs):
                        if pmid == "12345":
                            return 100
                        else:
                            raise Exception("Citation API error")

                    mock_citations.side_effect = citation_side_effect

                    # Mock impact factor service
                    with patch.object(
                        self.impact_factor_service, "get_impact_factor"
                    ) as mock_if:
                        mock_if.return_value = 5.0

                        # Process papers with error handling
                        scored_papers = []
                        for paper in sample_papers:
                            try:
                                citations = self.citation_service.get_citation_count(
                                    paper.pmid
                                )
                            except Exception:
                                citations = 0  # Fallback value

                            impact_factor = (
                                self.impact_factor_service.get_impact_factor(
                                    paper.journal
                                )
                            )

                            score = self.scoring_engine.calculate_paper_score(
                                paper, citations, impact_factor, topic.query
                            )

                            scored_paper = ScoredPaper(
                                pmid=paper.pmid,
                                title=paper.title,
                                authors=paper.authors,
                                journal=paper.journal,
                                publication_date=paper.publication_date,
                                citation_count=citations,
                                impact_factor=impact_factor,
                                score=score,
                                rank=0,
                            )
                            scored_papers.append(scored_paper)

                        # Should still have papers despite API failures
                        assert len(scored_papers) == 2
                        assert scored_papers[0].citation_count == 100
                        assert scored_papers[1].citation_count == 0  # Fallback

    def test_workflow_with_empty_results(self):
        """Test workflow behavior with empty search results."""
        topic = TopicConfig(
            name="empty-topic",
            query="nonexistent query",
            max_papers=10,
            essential_count=5,
            enabled=True,
        )

        # Mock empty search results
        with patch.object(self.paper_collector, "search_papers") as mock_search:
            mock_search.return_value = []

            with patch.object(
                self.paper_collector, "get_paper_details"
            ) as mock_details:
                mock_details.return_value = []

                # Process empty results
                papers = self.paper_collector.get_paper_details([])
                assert papers == []

                # Score empty list
                scored_papers = self.scoring_engine.rank_papers([])
                assert scored_papers == []

                # Select from empty list
                essential_papers = self.scoring_engine.select_essential_papers(
                    [], topic.essential_count
                )
                assert essential_papers == []

                # GitHub issue should still be created with "no papers" message
                with patch.object(
                    self.github_manager, "create_or_update_issue"
                ) as mock_github:
                    mock_github.return_value = {
                        "number": 43,
                        "created": True,
                        "html_url": "https://github.com/test/repo/issues/43",
                    }

                    issue_result = self.github_manager.create_or_update_issue(
                        topic.name, essential_papers
                    )

                    assert issue_result["created"] is True
                    mock_github.assert_called_once_with(topic.name, [])

    def test_workflow_performance_metrics(self):
        """Test workflow performance measurement."""
        import time

        topic = TopicConfig(
            name="performance-test",
            query="performance test",
            max_papers=100,
            essential_count=10,
            enabled=True,
        )

        # Create larger dataset for performance testing
        sample_papers = []
        for i in range(50):
            paper = Paper(
                pmid=str(10000 + i),
                title=f"Test Paper {i}",
                authors=[f"Author {i}"],
                journal=f"Journal {i % 5}",
                publication_date=datetime(2023, 1, 1),
                doi=f"10.1000/test.{10000 + i}",
            )
            sample_papers.append(paper)

        with patch.object(self.paper_collector, "search_papers") as mock_search:
            mock_search.return_value = [str(10000 + i) for i in range(50)]

            with patch.object(
                self.paper_collector, "get_paper_details"
            ) as mock_details:
                mock_details.return_value = sample_papers

                # Mock fast citation and impact factor services
                with patch.object(
                    self.citation_service, "batch_get_citation_counts"
                ) as mock_citations:
                    citation_counts = {str(10000 + i): i * 10 for i in range(50)}
                    mock_citations.return_value = citation_counts

                    with patch.object(
                        self.impact_factor_service, "get_impact_factor"
                    ) as mock_if:
                        mock_if.return_value = 5.0

                        # Measure workflow performance
                        start_time = time.time()

                        # Score all papers
                        scored_papers = []
                        for paper in sample_papers:
                            citations = citation_counts.get(paper.pmid, 0)
                            impact_factor = 5.0

                            score = self.scoring_engine.calculate_paper_score(
                                paper, citations, impact_factor, topic.query
                            )

                            scored_paper = ScoredPaper(
                                pmid=paper.pmid,
                                title=paper.title,
                                authors=paper.authors,
                                journal=paper.journal,
                                publication_date=paper.publication_date,
                                citation_count=citations,
                                impact_factor=impact_factor,
                                score=score,
                                rank=0,
                            )
                            scored_papers.append(scored_paper)

                        # Rank and select
                        ranked_papers = self.scoring_engine.rank_papers(scored_papers)
                        essential_papers = self.scoring_engine.select_essential_papers(
                            ranked_papers, topic.essential_count
                        )

                        end_time = time.time()
                        processing_time = end_time - start_time

                        # Verify results
                        assert len(essential_papers) == topic.essential_count
                        assert processing_time < 5.0  # Should complete within 5 seconds

                        # Verify ranking is correct
                        for i in range(len(essential_papers) - 1):
                            assert (
                                essential_papers[i].score
                                >= essential_papers[i + 1].score
                            )

    def test_configuration_validation_integration(self):
        """Test integration of configuration validation with services."""
        # Test with invalid configuration
        invalid_topics_data = {
            "topics": [
                {
                    "name": "invalid-topic",
                    "query": "",  # Empty query
                    "max_papers": -10,  # Invalid value
                    "essential_count": 0,  # Invalid value
                    "enabled": True,
                }
            ]
        }

        invalid_topics_file = Path(self.temp_dir) / "invalid_topics.yaml"
        with open(invalid_topics_file, "w") as f:
            yaml.dump(invalid_topics_data, f)

        invalid_config_manager = ConfigurationManager(self.temp_dir)
        invalid_config_manager.topics_file = invalid_topics_file

        # Should raise validation error
        with pytest.raises(ValueError):
            invalid_config_manager.load_topics()

    def test_cache_integration(self):
        """Test integration with caching systems."""
        # Test that cache is properly used across services
        with patch.object(
            self.citation_service.cache_manager, "get_citation"
        ) as mock_cache_get:
            with patch.object(
                self.citation_service.cache_manager, "save_citation"
            ) as mock_cache_save:
                # First call - cache miss
                mock_cache_get.return_value = None

                with patch(
                    "src.pubmed_miner.services.citation_service.requests.get"
                ) as mock_requests:
                    mock_response = Mock()
                    mock_response.status_code = 200
                    mock_response.json.return_value = {
                        "message": {"is-referenced-by-count": 100}
                    }
                    mock_requests.return_value = mock_response

                    # First call should fetch from API and cache
                    count1 = self.citation_service.get_citation_count(
                        "12345", doi="10.1000/test"
                    )
                    assert count1 == 100
                    mock_cache_save.assert_called_once()

                # Second call - cache hit
                from src.pubmed_miner.models.cache import CitationCache

                cached_citation = CitationCache(
                    pmid="12345",
                    citation_count=100,
                    last_updated=datetime.now(),
                    source="crossref",
                )
                mock_cache_get.return_value = cached_citation

                count2 = self.citation_service.get_citation_count("12345")
                assert count2 == 100
                # Should not make additional API calls

    def test_error_propagation_integration(self):
        """Test how errors propagate through the integrated workflow."""
        from src.pubmed_miner.utils.error_handler import APIError

        topic = TopicConfig(
            name="error-test",
            query="error test",
            max_papers=5,
            essential_count=2,
            enabled=True,
        )

        # Test API error propagation
        with patch.object(self.paper_collector, "search_papers") as mock_search:
            mock_search.side_effect = APIError("PubMed API unavailable")

            with pytest.raises(APIError, match="PubMed API unavailable"):
                self.paper_collector.search_papers(topic.query)

        # Test validation error propagation
        with pytest.raises(ValueError, match="Query cannot be empty"):
            self.paper_collector.search_papers("")

    def test_concurrent_workflow_execution(self):
        """Test concurrent execution of workflow components."""
        import threading

        results = []
        errors = []

        def run_scoring_workflow(topic_name):
            try:
                # Mock data for concurrent test
                sample_paper = Paper(
                    pmid=f"pmid_{topic_name}",
                    title=f"Paper for {topic_name}",
                    authors=["Test Author"],
                    journal="Test Journal",
                    publication_date=datetime(2023, 1, 1),
                )

                score = self.scoring_engine.calculate_paper_score(
                    sample_paper, 100, 5.0, f"query for {topic_name}"
                )

                scored_paper = ScoredPaper(
                    pmid=sample_paper.pmid,
                    title=sample_paper.title,
                    authors=sample_paper.authors,
                    journal=sample_paper.journal,
                    publication_date=sample_paper.publication_date,
                    citation_count=100,
                    impact_factor=5.0,
                    score=score,
                    rank=1,
                )

                results.append(scored_paper)

            except Exception as e:
                errors.append(e)

        # Create multiple threads
        threads = []
        for i in range(10):
            thread = threading.Thread(target=run_scoring_workflow, args=(f"topic_{i}",))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # All workflows should complete successfully
        assert len(errors) == 0
        assert len(results) == 10

        # All results should have valid scores
        for result in results:
            assert 0 <= result.score <= 100
            assert result.rank == 1
