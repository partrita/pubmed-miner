"""
Unit tests for GitHubIssuesManager.
"""

import pytest
from unittest.mock import Mock, patch
from datetime import datetime

from src.pubmed_miner.services.github_manager import GitHubIssuesManager
from src.pubmed_miner.models import ScoredPaper, GitHubConfig
from src.pubmed_miner.utils.error_handler import GitHubError


class TestGitHubIssuesManager:
    """Test cases for GitHubIssuesManager class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.config = GitHubConfig(
            token="mock_token_for_local_testing",
            repository="testuser/testrepo",
            issue_labels=["essential-papers", "automated"],
        )
        self.manager = GitHubIssuesManager(self.config)

        # Sample scored papers for testing
        self.sample_papers = [
            ScoredPaper(
                pmid="12345",
                title="Machine Learning in Healthcare",
                authors=["John Doe", "Jane Smith"],
                journal="Nature Medicine",
                publication_date=datetime(2023, 1, 15),
                citation_count=150,
                impact_factor=87.2,
                score=95.5,
                rank=1,
            ),
            ScoredPaper(
                pmid="67890",
                title="AI Applications in Drug Discovery",
                authors=["Bob Johnson"],
                journal="Science",
                publication_date=datetime(2023, 2, 10),
                citation_count=75,
                impact_factor=56.9,
                score=88.2,
                rank=2,
            ),
        ]

    def test_initialization(self):
        """Test manager initialization."""
        assert self.manager.config.repository == "testuser/testrepo"
        assert self.manager.base_url == "https://api.github.com"
        assert self.manager.mock_mode is True  # Should be in mock mode
        assert "Accept" in self.manager.headers

    def test_mock_mode_create_issue(self):
        """Test issue creation in mock mode."""
        result = self.manager.create_or_update_issue("test-topic", self.sample_papers)

        assert result["mock_mode"] is True
        assert result["created"] is True
        assert "number" in result
        assert result["title"] == "[Essential Papers] test-topic"
        assert "html_url" in result
        assert "body" in result

    def test_mock_mode_find_existing_issue(self):
        """Test finding existing issue in mock mode."""
        issue = self.manager.find_existing_issue("test-topic")
        assert issue is None  # Mock mode always returns None

    def test_mock_mode_validate_access(self):
        """Test access validation in mock mode."""
        assert self.manager.validate_access() is True

    @patch("src.pubmed_miner.services.github_manager.requests.get")
    def test_find_existing_issue_found(self, mock_get):
        """Test finding an existing issue."""
        # Mock GitHub API response with existing issue
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = [
            {
                "number": 42,
                "title": "[Essential Papers] machine-learning-healthcare",
                "state": "open",
                "body": "Existing issue body",
                "html_url": "https://github.com/testuser/testrepo/issues/42",
            }
        ]
        mock_get.return_value = mock_response

        issue = self.manager._find_existing_issue("machine-learning-healthcare")

        assert issue is not None
        assert issue["number"] == 42
        assert issue["title"] == "[Essential Papers] machine-learning-healthcare"

    @patch("src.pubmed_miner.services.github_manager.requests.get")
    def test_find_existing_issue_not_found(self, mock_get):
        """Test when no existing issue is found."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = []
        mock_get.return_value = mock_response

        issue = self.manager._find_existing_issue("nonexistent-topic")
        assert issue is None

    @patch("src.pubmed_miner.services.github_manager.requests.get")
    def test_find_existing_issue_api_error(self, mock_get):
        """Test API error when finding existing issue."""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_response.text = "Repository not found"
        mock_get.return_value = mock_response

        with pytest.raises(GitHubError, match="Failed to search for existing issues"):
            self.manager._find_existing_issue("test-topic")

    @patch("src.pubmed_miner.services.github_manager.requests.post")
    def test_create_issue_success(self, mock_post):
        """Test successful issue creation."""
        mock_response = Mock()
        mock_response.status_code = 201
        mock_response.json.return_value = {
            "number": 43,
            "title": "[Essential Papers] test-topic",
            "html_url": "https://github.com/testuser/testrepo/issues/43",
        }
        mock_post.return_value = mock_response

        issue = self.manager._create_issue("test-topic", "Issue body content")

        assert issue["number"] == 43
        assert issue["title"] == "[Essential Papers] test-topic"

        # Verify API call
        mock_post.assert_called_once()
        call_args = mock_post.call_args
        assert call_args[0][0].endswith("/repos/testuser/testrepo/issues")

        # Check request body
        request_data = call_args[1]["json"]
        assert request_data["title"] == "[Essential Papers] test-topic"
        assert request_data["body"] == "Issue body content"
        assert request_data["labels"] == ["essential-papers", "automated"]

    @patch("src.pubmed_miner.services.github_manager.requests.patch")
    def test_update_issue_success(self, mock_patch):
        """Test successful issue update."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "number": 42,
            "title": "[Essential Papers] test-topic",
            "html_url": "https://github.com/testuser/testrepo/issues/42",
        }
        mock_patch.return_value = mock_response

        issue = self.manager._update_issue(42, "Updated body content")

        assert issue["number"] == 42
        mock_patch.assert_called_once()

    def test_format_issue_body(self):
        """Test issue body formatting."""
        body = self.manager._format_issue_body("test-topic", self.sample_papers)

        # Check that body contains expected elements
        assert "test-topic" in body
        assert "Machine Learning in Healthcare" in body
        assert "AI Applications in Drug Discovery" in body
        assert "John Doe" in body
        assert "Nature Medicine" in body
        assert "PMID: 12345" in body
        assert "Score: 95.5" in body
        assert "Rank #1" in body

        # Check markdown formatting
        assert "## Essential Papers for test-topic" in body
        assert "### 1. Machine Learning in Healthcare" in body
        assert "**Authors:** John Doe, Jane Smith" in body
        assert "**Journal:** Nature Medicine" in body

    def test_format_issue_body_empty_papers(self):
        """Test issue body formatting with empty papers list."""
        body = self.manager._format_issue_body("test-topic", [])

        assert "No essential papers found" in body
        assert "test-topic" in body

    @patch.object(GitHubIssuesManager, "_find_existing_issue")
    @patch.object(GitHubIssuesManager, "_create_issue")
    def test_create_or_update_issue_new(self, mock_create, mock_find):
        """Test creating a new issue when none exists."""
        # Mock no existing issue
        mock_find.return_value = None

        # Mock successful creation
        mock_create.return_value = {
            "number": 43,
            "created": True,
            "html_url": "https://github.com/testuser/testrepo/issues/43",
        }

        result = self.manager.create_or_update_issue("test-topic", self.sample_papers)

        assert result["created"] is True
        assert result["number"] == 43
        mock_find.assert_called_once_with("test-topic")
        mock_create.assert_called_once()

    @patch.object(GitHubIssuesManager, "_find_existing_issue")
    @patch.object(GitHubIssuesManager, "_update_issue")
    def test_create_or_update_issue_existing(self, mock_update, mock_find):
        """Test updating an existing issue."""
        # Mock existing issue
        mock_find.return_value = {
            "number": 42,
            "title": "[Essential Papers] test-topic",
        }

        # Mock successful update
        mock_update.return_value = {
            "number": 42,
            "updated": True,
            "html_url": "https://github.com/testuser/testrepo/issues/42",
        }

        result = self.manager.create_or_update_issue("test-topic", self.sample_papers)

        assert result["updated"] is True
        assert result["number"] == 42
        mock_find.assert_called_once_with("test-topic")
        mock_update.assert_called_once_with(42, mock_update.call_args[0][1])

    def test_create_or_update_issue_empty_papers(self):
        """Test creating/updating issue with empty papers list."""
        with patch.object(self.manager, "_find_existing_issue") as mock_find:
            mock_find.return_value = None

            with patch.object(self.manager, "_create_issue") as mock_create:
                mock_create.return_value = {"number": 43, "created": True}

                result = self.manager.create_or_update_issue("test-topic", [])

                # Should still create issue with "no papers found" message
                assert result["created"] is True
                mock_create.assert_called_once()

    @patch("src.pubmed_miner.services.github_manager.requests.post")
    def test_add_comment_to_issue(self, mock_post):
        """Test adding a comment to an issue."""
        mock_response = Mock()
        mock_response.status_code = 201
        mock_response.json.return_value = {"id": 123456, "body": "Test comment"}
        mock_post.return_value = mock_response

        comment = self.manager.add_comment_to_issue(42, "Test comment")

        assert comment["id"] == 123456
        assert comment["body"] == "Test comment"
        mock_post.assert_called_once()

    def test_validate_topic_name(self):
        """Test topic name validation."""
        # Valid topic names
        assert self.manager._validate_topic_name("machine-learning") is True
        assert self.manager._validate_topic_name("covid_19_research") is True
        assert self.manager._validate_topic_name("ai123") is True

        # Invalid topic names
        assert self.manager._validate_topic_name("") is False
        assert self.manager._validate_topic_name("topic with spaces") is False
        assert self.manager._validate_topic_name("topic@special") is False
        assert self.manager._validate_topic_name(None) is False

    def test_sanitize_topic_name(self):
        """Test topic name sanitization."""
        assert (
            self.manager._sanitize_topic_name("Machine Learning") == "machine-learning"
        )
        assert (
            self.manager._sanitize_topic_name("COVID-19 Research")
            == "covid-19-research"
        )
        assert self.manager._sanitize_topic_name("AI & ML") == "ai-ml"

    @patch("src.pubmed_miner.services.github_manager.requests.get")
    def test_get_repository_info(self, mock_get):
        """Test getting repository information."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "name": "testrepo",
            "full_name": "testuser/testrepo",
            "open_issues_count": 5,
            "private": False,
        }
        mock_get.return_value = mock_response

        repo_info = self.manager.get_repository_info()

        assert repo_info["name"] == "testrepo"
        assert repo_info["open_issues_count"] == 5

    @patch("src.pubmed_miner.services.github_manager.requests.get")
    def test_get_repository_info_not_found(self, mock_get):
        """Test getting repository info when repo doesn't exist."""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_response.text = "Not Found"
        mock_get.return_value = mock_response

        with pytest.raises(GitHubError, match="Repository not found"):
            self.manager.get_repository_info()

    def test_get_statistics(self):
        """Test statistics collection."""
        stats = self.manager.get_statistics()

        required_keys = [
            "issues_created",
            "issues_updated",
            "comments_added",
            "api_calls",
            "error_count",
            "last_activity",
        ]

        for key in required_keys:
            assert key in stats

    @patch("src.pubmed_miner.services.github_manager.requests.get")
    def test_list_issues_for_topic(self, mock_get):
        """Test listing issues for a specific topic."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = [
            {
                "number": 42,
                "title": "[Essential Papers] machine-learning",
                "state": "open",
            },
            {
                "number": 43,
                "title": "[Essential Papers] covid-research",
                "state": "open",
            },
        ]
        mock_get.return_value = mock_response

        issues = self.manager.list_issues_for_topic("machine-learning")

        assert len(issues) == 1
        assert issues[0]["number"] == 42

    def test_format_paper_summary(self):
        """Test individual paper summary formatting."""
        paper = self.sample_papers[0]
        summary = self.manager._format_paper_summary(paper, 1)

        assert "### 1. Machine Learning in Healthcare" in summary
        assert "**Authors:** John Doe, Jane Smith" in summary
        assert "**Journal:** Nature Medicine (IF: 87.2)" in summary
        assert "**Citations:** 150" in summary
        assert "**Score:** 95.5" in summary
        assert "[PubMed](https://pubmed.ncbi.nlm.nih.gov/12345/)" in summary

    def test_rate_limiting_handling(self):
        """Test GitHub API rate limiting handling."""
        with patch("src.pubmed_miner.services.github_manager.requests.get") as mock_get:
            # Mock rate limit response
            mock_response = Mock()
            mock_response.status_code = 403
            mock_response.headers = {
                "X-RateLimit-Remaining": "0",
                "X-RateLimit-Reset": str(int(datetime.now().timestamp()) + 3600),
            }
            mock_response.json.return_value = {"message": "API rate limit exceeded"}
            mock_get.return_value = mock_response

            with pytest.raises(GitHubError, match="GitHub API rate limit exceeded"):
                self.manager._find_existing_issue("test-topic")

    def test_authentication_error(self):
        """Test handling of authentication errors."""
        with patch("src.pubmed_miner.services.github_manager.requests.get") as mock_get:
            mock_response = Mock()
            mock_response.status_code = 401
            mock_response.text = "Bad credentials"
            mock_get.return_value = mock_response

            with pytest.raises(GitHubError, match="GitHub authentication failed"):
                self.manager._find_existing_issue("test-topic")

    def test_close_issue(self):
        """Test closing an issue."""
        with patch(
            "src.pubmed_miner.services.github_manager.requests.patch"
        ) as mock_patch:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.return_value = {"number": 42, "state": "closed"}
            mock_patch.return_value = mock_response

            result = self.manager.close_issue(42)

            assert result["state"] == "closed"
            mock_patch.assert_called_once()

    def test_reopen_issue(self):
        """Test reopening a closed issue."""
        with patch(
            "src.pubmed_miner.services.github_manager.requests.patch"
        ) as mock_patch:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.return_value = {"number": 42, "state": "open"}
            mock_patch.return_value = mock_response

            result = self.manager.reopen_issue(42)

            assert result["state"] == "open"
            mock_patch.assert_called_once()

    def test_concurrent_operations(self):
        """Test thread safety of concurrent operations."""
        import threading

        results = []
        errors = []

        def create_issue(topic):
            try:
                with patch.object(self.manager, "_find_existing_issue") as mock_find:
                    mock_find.return_value = None
                    with patch.object(self.manager, "_create_issue") as mock_create:
                        mock_create.return_value = {"number": 42, "created": True}
                        result = self.manager.create_or_update_issue(
                            topic, self.sample_papers
                        )
                        results.append(result)
            except Exception as e:
                errors.append(e)

        # Create multiple threads
        threads = []
        for i in range(5):
            thread = threading.Thread(target=create_issue, args=(f"topic-{i}",))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # All operations should succeed
        assert len(errors) == 0
        assert len(results) == 5

    def test_mock_mode_initialization(self):
        """Test initialization in mock mode."""
        mock_config = GitHubConfig(
            token="mock_token_for_local_testing",
            repository="test/repo",
            issue_labels=["test"],
        )
        mock_manager = GitHubIssuesManager(mock_config)

        assert mock_manager.mock_mode is True
        assert "Authorization" not in mock_manager.headers

    def test_empty_token_uses_mock_mode(self):
        """Test that empty token automatically uses mock mode."""
        mock_config = GitHubConfig(
            token="", repository="test/repo", issue_labels=["test"]
        )
        mock_manager = GitHubIssuesManager(mock_config)

        assert mock_manager.mock_mode is True

    def test_none_token_uses_mock_mode(self):
        """Test that None token automatically uses mock mode."""
        mock_config = GitHubConfig(
            token=None, repository="test/repo", issue_labels=["test"]
        )
        mock_manager = GitHubIssuesManager(mock_config)

        assert mock_manager.mock_mode is True
