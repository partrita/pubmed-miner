"""
GitHub Issues management service.
"""

import re
import time
import logging
from datetime import datetime
from typing import List, Dict, Optional, Any, Tuple, Union

import requests
from requests.exceptions import HTTPError

from ..models import ScoredPaper, GitHubConfig

logger = logging.getLogger(__name__)


class GitHubError(Exception):
    """Exception for GitHub API errors."""

    pass


class GitHubIssuesManager:
    """Manages GitHub issues for tracking essential papers."""

    def __init__(
        self,
        config: Optional[GitHubConfig] = None,
        token: Optional[str] = None,
        mock_mode: bool = False,
    ):
        """Initialize the GitHub issues manager.

        Args:
            config: GitHub configuration object
            token: GitHub personal access token
            mock_mode: If True, operate in mock mode without API calls
        """
        self.logger = logging.getLogger(__name__)
        self.config = config
        self.token = token
        self.mock_mode = mock_mode
        self._rate_limit_remaining = 5000
        self._rate_limit_reset = 0
        self._last_request_time = 0
        self._min_request_interval = 1.0
        self.headers = {
            "Accept": "application/vnd.github.v3+json",
            "User-Agent": "PubMed-Miner",
        }
        self.base_url = "https://api.github.com"

        if token:
            self.headers["Authorization"] = f"token {token}"
            self.logger.info("GitHub token configured")

        if not config:
            config = GitHubConfig()
            self.config = config

        if not self.token and self.config:
            self.token = self.config.token

        # Auto-detect mock mode: when token is absent or is a placeholder/local-testing token
        if mock_mode is not True and (not self.token or str(self.token).startswith("mock_token")):
            self.mock_mode = True

        if self.mock_mode:
            self.logger.info("GitHub manager initialized in mock mode")
        else:
            self.logger.info(
                f"GitHub manager initialized for repository: {self.config.repository}"
            )

        self._stats_issues_created = 0
        self._stats_issues_updated = 0
        self._stats_comments_added = 0
        self._stats_api_calls = 0
        self._stats_rate_limit_hits = 0
        self._stats_errors = 0

    def create_issue(
        self, title: str, body: str, labels: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """Create a new GitHub issue.

        Args:
            title: Issue title
            body: Issue body content
            labels: List of label names

        Returns:
            Created issue data

        Raises:
            GitHubError: If issue creation fails
        """
        if self.mock_mode:
            return self._mock_create_issue(title, body, labels)
        return self._create_issue_api(title, body, labels)

    def _mock_create_issue(
        self, title: str, body: str, labels: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """Create a mock issue for testing.

        Args:
            title: Issue title
            body: Issue body content
            labels: List of label names

        Returns:
            Mock issue data
        """
        issue_number = 1000 + (self._stats_issues_created % 1000)
        result = {
            "number": issue_number,
            "title": title,
            "body": body,
            "state": "open",
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "created_at": datetime.now().isoformat(),
            "labels": labels if labels else [],
            "mock_mode": True,
        }
        self._stats_issues_created += 1
        self.logger.debug(f"Mock issue created: #{issue_number}")
        return result

    def _create_issue_api(
        self, title: str, body: str, labels: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """Create an issue via GitHub API.

        Args:
            title: Issue title
            body: Issue body content
            labels: List of label names

        Returns:
            Created issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues"
        data: Dict[str, Any] = {
            "title": title,
            "body": body,
            "labels": labels if labels else ["essential-papers"],
        }

        try:
            response = requests.post(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            issue = response.json()
            self._stats_issues_created += 1
            self._stats_api_calls += 1
            self.logger.info(f"Created issue: #{issue['number']}")
            return issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to create issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to create issue: {e}")

    def update_issue(
        self, issue_number: int, title: str, body: str
    ) -> Dict[str, Any]:
        """Update an existing GitHub issue.

        Args:
            issue_number: Issue number to update
            title: New title
            body: New body content

        Returns:
            Updated issue data

        Raises:
            GitHubError: If update fails
        """
        if self.mock_mode:
            return self._mock_update_issue(issue_number, title, body)
        return self._update_issue_api(issue_number, title, body)

    def _mock_update_issue(
        self, issue_number: int, title: str, body: str
    ) -> Dict[str, Any]:
        """Update a mock issue.

        Args:
            issue_number: Issue number
            title: New title
            body: New body content

        Returns:
            Mock issue data
        """
        result = {
            "number": issue_number,
            "title": title,
            "body": body,
            "state": "open",
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "updated_at": datetime.now().isoformat(),
            "mock_mode": True,
        }
        self._stats_issues_updated += 1
        self.logger.debug(f"Mock issue updated: #{issue_number}")
        return result

    def _update_issue_api(
        self, issue_number: int, title: str, body: str
    ) -> Dict[str, Any]:
        """Update an issue via GitHub API.

        Args:
            issue_number: Issue number
            title: New title
            body: New body content

        Returns:
            Updated issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}"
        data: Dict[str, Any] = {
            "title": title,
            "body": body,
        }

        try:
            response = requests.patch(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            issue = response.json()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Updated issue: #{issue_number}")
            return issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to update issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to update issue: {e}")

    def close_issue(self, issue_number: int) -> Dict[str, Any]:
        """Close a GitHub issue.

        Args:
            issue_number: Issue number to close

        Returns:
            Updated issue data

        Raises:
            GitHubError: If closing fails
        """
        if self.mock_mode:
            return self._mock_close_issue(issue_number)
        return self._close_issue_api(issue_number)

    def _mock_close_issue(self, issue_number: int) -> Dict[str, Any]:
        """Close a mock issue.

        Args:
            issue_number: Issue number

        Returns:
            Mock issue data
        """
        result = {
            "number": issue_number,
            "state": "closed",
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "mock_mode": True,
        }
        self._stats_issues_updated += 1
        self.logger.debug(f"Mock issue closed: #{issue_number}")
        return result

    def _close_issue_api(self, issue_number: int) -> Dict[str, Any]:
        """Close an issue via GitHub API.

        Args:
            issue_number: Issue number

        Returns:
            Updated issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}"
        data: Dict[str, Any] = {"state": "closed"}

        try:
            response = requests.patch(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            issue = response.json()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Closed issue: #{issue_number}")
            return issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to close issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to close issue: {e}")

    def reopen_issue(self, issue_number: int) -> Dict[str, Any]:
        """Reopen a closed GitHub issue.

        Args:
            issue_number: Issue number to reopen

        Returns:
            Updated issue data

        Raises:
            GitHubError: If reopening fails
        """
        if self.mock_mode:
            return self._mock_reopen_issue(issue_number)
        return self._reopen_issue_api(issue_number)

    def _mock_reopen_issue(self, issue_number: int) -> Dict[str, Any]:
        """Reopen a mock issue.

        Args:
            issue_number: Issue number

        Returns:
            Mock issue data
        """
        result = {
            "number": issue_number,
            "state": "open",
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "mock_mode": True,
        }
        self._stats_issues_updated += 1
        self.logger.debug(f"Mock issue reopened: #{issue_number}")
        return result

    def _reopen_issue_api(self, issue_number: int) -> Dict[str, Any]:
        """Reopen an issue via GitHub API.

        Args:
            issue_number: Issue number

        Returns:
            Updated issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}"
        data: Dict[str, Any] = {"state": "open"}

        try:
            response = requests.patch(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            issue = response.json()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Reopened issue: #{issue_number}")
            return issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to reopen issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to reopen issue: {e}")

    def add_comment_to_issue(
        self, issue_number: int, comment: str
    ) -> Dict[str, Any]:
        """Add a comment to a GitHub issue.

        Args:
            issue_number: Issue number
            comment: Comment text

        Returns:
            Created comment data

        Raises:
            GitHubError: If adding comment fails
        """
        if self.mock_mode:
            return self._mock_add_comment(issue_number, comment)
        return self._add_comment_api(issue_number, comment)

    def _mock_add_comment(
        self, issue_number: int, comment: str
    ) -> Dict[str, Any]:
        """Add a mock comment.

        Args:
            issue_number: Issue number
            comment: Comment text

        Returns:
            Mock comment data
        """
        result = {
            "id": 5000000 + (self._stats_comments_added % 1000000),
            "body": comment,
            "issue_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "created_at": datetime.now().isoformat(),
            "mock_mode": True,
        }
        self._stats_comments_added += 1
        self.logger.debug(f"Mock comment added to issue #{issue_number}")
        return result

    def _add_comment_api(
        self, issue_number: int, comment: str
    ) -> Dict[str, Any]:
        """Add a comment via GitHub API.

        Args:
            issue_number: Issue number
            comment: Comment text

        Returns:
            Created comment data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}/comments"
        data: Dict[str, Any] = {"body": comment}

        try:
            response = requests.post(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            comment_data = response.json()
            self._stats_comments_added += 1
            self._stats_api_calls += 1
            self.logger.info(f"Added comment to issue #{issue_number}")
            return comment_data
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to add comment: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to add comment: {e}")

    def add_labels_to_issue(
        self, issue_number: int, labels: List[str]
    ) -> Dict[str, Any]:
        """Add labels to a GitHub issue.

        Args:
            issue_number: Issue number
            labels: List of label names

        Returns:
            Updated issue data

        Raises:
            GitHubError: If adding labels fails
        """
        if self.mock_mode:
            return self._mock_add_labels(issue_number, labels)
        return self._add_labels_api(issue_number, labels)

    def _mock_add_labels(
        self, issue_number: int, labels: List[str]
    ) -> Dict[str, Any]:
        """Add mock labels.

        Args:
            issue_number: Issue number
            labels: List of label names

        Returns:
            Mock issue data
        """
        result = {
            "number": issue_number,
            "labels": labels,
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "mock_mode": True,
        }
        self._stats_issues_updated += 1
        self.logger.debug(f"Mock labels added to issue #{issue_number}")
        return result

    def _add_labels_api(
        self, issue_number: int, labels: List[str]
    ) -> Dict[str, Any]:
        """Add labels via GitHub API.

        Args:
            issue_number: Issue number
            labels: List of label names

        Returns:
            Updated issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}/labels"

        try:
            response = requests.post(
                url, headers=self.headers, json=labels, timeout=30
            )
            response.raise_for_status()
            label_data = response.json()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Added labels to issue #{issue_number}")
            return label_data
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to add labels: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to add labels: {e}")

    def remove_labels_from_issue(
        self, issue_number: int, labels: List[str]
    ) -> Dict[str, Any]:
        """Remove labels from a GitHub issue.

        Args:
            issue_number: Issue number
            labels: List of label names to remove

        Returns:
            Updated issue data

        Raises:
            GitHubError: If removing labels fails
        """
        if self.mock_mode:
            return self._mock_remove_labels(issue_number, labels)
        return self._remove_labels_api(issue_number, labels)

    def _mock_remove_labels(
        self, issue_number: int, labels: List[str]
    ) -> Dict[str, Any]:
        """Remove mock labels.

        Args:
            issue_number: Issue number
            labels: List of label names

        Returns:
            Mock issue data
        """
        result = {
            "number": issue_number,
            "labels": [],
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "mock_mode": True,
        }
        self._stats_issues_updated += 1
        self.logger.debug(f"Mock labels removed from issue #{issue_number}")
        return result

    def _remove_labels_api(
        self, issue_number: int, labels: List[str]
    ) -> Dict[str, Any]:
        """Remove labels via GitHub API.

        Args:
            issue_number: Issue number
            labels: List of label names to remove

        Returns:
            Updated issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}/labels"

        try:
            response = requests.delete(
                url, headers=self.headers, json=labels, timeout=30
            )
            response.raise_for_status()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Removed labels from issue #{issue_number}")
            return {}
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to remove labels: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to remove labels: {e}")

    def get_issue(self, issue_number: int) -> Dict[str, Any]:
        """Get a GitHub issue by number.

        Args:
            issue_number: Issue number

        Returns:
            Issue data

        Raises:
            GitHubError: If retrieval fails
        """
        if self.mock_mode:
            return self._mock_get_issue(issue_number)
        return self._get_issue_api(issue_number)

    def _mock_get_issue(self, issue_number: int) -> Dict[str, Any]:
        """Get a mock issue.

        Args:
            issue_number: Issue number

        Returns:
            Mock issue data
        """
        result = {
            "number": issue_number,
            "title": f"Mock issue #{issue_number}",
            "state": "open",
            "body": "Mock issue body",
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "mock_mode": True,
        }
        self._stats_api_calls += 1
        return result

    def _get_issue_api(self, issue_number: int) -> Dict[str, Any]:
        """Get an issue via GitHub API.

        Args:
            issue_number: Issue number

        Returns:
            Issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}"

        try:
            response = requests.get(
                url, headers=self.headers, timeout=30
            )
            response.raise_for_status()
            issue = response.json()
            self._stats_api_calls += 1
            return issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            status = getattr(getattr(e, "response", None), "status_code", None)
            if status is None:
                import re
                m = re.search(r"(\d{3})", str(e))
                if m:
                    status = int(m.group(1))
            if status == 403:
                raise GitHubError("GitHub API rate limit exceeded")
            if status == 401:
                raise GitHubError("GitHub authentication failed")
            raise GitHubError(f"Failed to get issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to get issue: {e}")

    def list_issues_for_topic(
        self, topic: str, state: str = "all"
    ) -> List[Dict[str, Any]]:
        """List GitHub issues filtered by topic label.

        Args:
            topic: Topic name to filter by
            state: Issue state filter ('open', 'closed', 'all')

        Returns:
            List of issues

        Raises:
            GitHubError: If retrieval fails
        """
        if self.mock_mode:
            return self._mock_list_issues(topic, state)
        return self._list_issues_api(topic, state)

    def _mock_list_issues(
        self, topic: str, state: str = "all"
    ) -> List[Dict[str, Any]]:
        """List mock issues.

        Args:
            topic: Topic name
            state: State filter

        Returns:
            List of mock issues
        """
        issues = []
        for i in range(1000, 1003):
            issue = {
                "number": i,
                "title": f"[{self.config.issue_prefix}] {topic}",
                "state": "open" if state == "all" else state,
                "labels": [f"topic-{topic}"],
                "html_url": f"https://github.com/{self.config.repository}/issues/{i}",
                "mock_mode": True,
            }
            issues.append(issue)
        self._stats_api_calls += 1
        return issues

    def _list_issues_api(
        self, topic: str, state: str = "all"
    ) -> List[Dict[str, Any]]:
        """List issues via GitHub API.

        Args:
            topic: Topic name to filter by
            state: Issue state filter

        Returns:
            List of issues

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues"
        params: Dict[str, Any] = {
            "state": state,
            "labels": f"topic-{topic}",
            "per_page": 100,
        }

        try:
            response = requests.get(
                url, headers=self.headers, params=params, timeout=30
            )
            response.raise_for_status()
            issues = response.json()
            self._stats_api_calls += 1
            # Defense: filter by the requested topic label even if server did
            topic_label = f"topic-{topic}"
            issues = [i for i in issues if topic_label in i.get("labels", [])]
            return issues
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            status = getattr(getattr(e, "response", None), "status_code", None)
            if status is None:
                import re
                m = re.search(r"(\d{3})", str(e))
                if m:
                    status = int(m.group(1))
            if status == 403:
                raise GitHubError("GitHub API rate limit exceeded")
            if status == 401:
                raise GitHubError("GitHub authentication failed")
            raise GitHubError(f"Failed to list issues: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to list issues: {e}")

    def get_repository_info(self) -> Dict[str, Any]:
        """Get repository information from GitHub.

        Returns:
            Repository data

        Raises:
            GitHubError: If retrieval fails
        """
        if self.mock_mode:
            return self._mock_get_repository_info()
        return self._get_repository_info_api()

    def _mock_get_repository_info(self) -> Dict[str, Any]:
        """Get mock repository info.

        Returns:
            Mock repository data
        """
        result = {
            "name": self.config.repository.split("/")[-1],
            "full_name": self.config.repository,
            "open_issues_count": 42,
            "private": False,
            "mock_mode": True,
        }
        self._stats_api_calls += 1
        return result

    def _get_repository_info_api(self) -> Dict[str, Any]:
        """Get repository info via GitHub API.

        Returns:
            Repository data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}"
        self.logger.info(f"Fetching repository info via API: {url}")

        try:
            response = requests.get(
                url, headers=self.headers, timeout=30
            )
            self._stats_api_calls += 1
            self.logger.info(
                f"Repository info API response: {response.status_code}"
            )
            response.raise_for_status()
            repo_data = response.json()
            self.logger.info(f"Successfully fetched repository info")
            return repo_data
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            status = getattr(getattr(e, "response", None), "status_code", None)
            if status is None:
                import re
                m = re.search(r"(\d{3})", str(e))
                if m:
                    status = int(m.group(1))
            if status == 403:
                raise GitHubError("GitHub API rate limit exceeded")
            if status == 401:
                raise GitHubError("GitHub authentication failed")
            self.logger.error(f"Failed to get repository info: {e}")
            raise GitHubError(f"Failed to get repository info: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            self.logger.error(f"Failed to get repository info: {e}")
            raise GitHubError(f"Failed to get repository info: {e}")

    def validate_access(self) -> bool:
        """Validate GitHub access by checking rate limit.

        Returns:
            True if access is valid, False otherwise
        """
        if self.mock_mode:
            return True
        return self._validate_access_api()

    def _validate_access_api(self) -> bool:
        """Validate access via GitHub API.

        Returns:
            True if access is valid, False otherwise
        """
        try:
            url = f"{self.base_url}/rate_limit"
            response = requests.get(url, headers=self.headers, timeout=10)
            if response.status_code == 200:
                data = response.json()
                self._rate_limit_remaining = data.get("resources", {}).get(
                    "core", {}
                ).get("remaining", 5000)
                self._rate_limit_reset = data.get("resources", {}).get(
                    "core", {}
                ).get("reset", 0)
                return True
            return False
        except Exception as e:
            self.logger.error(f"GitHub access validation failed: {e}")
            return False

    def validate_and_sanitize_topic(
        self, topic: str
    ) -> Tuple[str, str]:
        """Validate and sanitize a topic name.

        Args:
            topic: Topic name to validate and sanitize

        Returns:
            Tuple of (validated_topic, sanitized_topic)

        Raises:
            ValueError: If topic is invalid
        """
        if not self._validate_topic_name(topic):
            raise ValueError(
                f"Invalid topic name: '{topic}'. "
                "Topic names must be alphanumeric with hyphens or underscores."
            )
        sanitized = self._sanitize_topic_name(topic)
        return topic, sanitized

    def _validate_topic_name(self, topic: str) -> bool:
        """Validate a topic name.

        Args:
            topic: Topic name to validate

        Returns:
            True if valid, False otherwise
        """
        if not topic or not isinstance(topic, str):
            return False
        if len(topic) < 1 or len(topic) > 100:
            return False
        pattern = r"^[a-zA-Z0-9_-]+$"
        return bool(re.match(pattern, topic))

    def _sanitize_topic_name(self, topic: str) -> str:
        """Sanitize a topic name for use in GitHub issues.

        Args:
            topic: Topic name to sanitize

        Returns:
            Sanitized topic name
        """
        sanitized = topic.lower()
        sanitized = re.sub(r"\s+", "-", sanitized)
        sanitized = re.sub(r"[^a-z0-9\-_]", "", sanitized)
        sanitized = re.sub(r"-+", "-", sanitized)
        sanitized = sanitized.strip("-_")
        if not sanitized:
            sanitized = "untitled"
        return sanitized

    def _format_issue_body(
        self, topic: str, papers: List[ScoredPaper]
    ) -> str:
        """Format issue body with paper list.

        Args:
            topic: Topic name
            papers: List of ScoredPaper objects

        Returns:
            Formatted issue body string
        """
        date_str = datetime.now().strftime("%Y-%m-%d")
        body = f"## Essential Papers for {topic}\n\n"
        body += f"**Generated:** {date_str}\n\n"

        if not papers:
            body += "*No essential papers found for this topic.*\n"
            return body

        body += f"**Total Papers:** {len(papers)}\n\n"

        for rank, paper in enumerate(papers[:20], start=1):
            body += self._format_paper_summary(paper, rank)

        if len(papers) > 20:
            body += f"\n*...and {len(papers) - 20} more papers*\n"

        return body

    def _format_paper_summary(
        self, paper: ScoredPaper, rank: int
    ) -> str:
        """Format a paper summary for issue body.

        Args:
            paper: ScoredPaper object
            rank: Paper rank (1-based)

        Returns:
            Formatted paper summary string
        """
        summary = (
            f"### {rank}. {paper.title}\n\n"
            f"- **PMID:** {paper.pmid}\n"
            f"- **Authors:** {', '.join(paper.authors) if paper.authors else 'Unknown'}\n"
            f"- **Journal:** {paper.journal}"
            + (f" (IF: {paper.impact_factor:.1f})" if paper.impact_factor is not None else "")
            + "\n"
            + (f"- **Published:** {paper.publication_date.strftime('%Y-%m-%d')}\n" if paper.publication_date else "")
            + f"- **Score:** {paper.score:.1f}\n"
            + (f"- **Citations:** {paper.citation_count}\n" if paper.citation_count is not None else "")
            + f"- **Rank:** #{rank}\n"
            + f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/{paper.pmid}/)\n"
            + "\n"
        )
        return summary

    def create_or_update_issue(
        self, topic: str, papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Create or update a GitHub issue for a topic.

        Args:
            topic: Topic name for the issue
            papers: List of scored papers

        Returns:
            Created or updated issue data

        Raises:
            GitHubError: If operation fails
            ValueError: If topic is invalid
        """
        validated_topic, sanitized_topic = self.validate_and_sanitize_topic(topic)

        date_str = datetime.now().strftime("%Y-%m-%d")
        issue_title = f"{date_str}: {sanitized_topic} papers"

        existing_issue = self.find_existing_issue_for_date(
            sanitized_topic, date_str
        )

        if existing_issue:
            self.logger.info(
                f"Updating existing issue #{existing_issue['number']} for {topic}"
            )
            issue_body = self._format_issue_body(topic, papers)
            updated = self._update_issue(
                existing_issue, papers
            )
            self._stats_issues_updated += 1
            return updated

        self.logger.info(f"Creating new issue for topic: {topic}")
        issue_body = self._format_issue_body(topic, papers)
        created = self._create_issue(issue_title, issue_body)
        self._stats_issues_created += 1
        return created

    def find_existing_issue_for_date(
        self, topic: str, date_str: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """Find an existing issue for a topic on a specific date.

        Args:
            topic: Topic name
            date_str: Date string in YYYY-MM-DD format

        Returns:
            Existing issue data or None
        """
        if self.mock_mode:
            return self._mock_find_existing_issue(topic, date_str)
        return self._find_existing_issue_api(topic, date_str)

    def _mock_find_existing_issue(
        self, topic: str, date_str: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """Find mock existing issue.

        Args:
            topic: Topic name
            date_str: Date string

        Returns:
            Mock issue or None
        """
        self._stats_api_calls += 1
        return None

    def find_existing_issue(self, topic: str) -> Optional[Dict[str, Any]]:
        """Find an existing issue for a topic (mock-safe public entry)."""
        if self.mock_mode:
            return None
        return self._find_existing_issue(topic)

    def _find_existing_issue(self, topic: str) -> Optional[Dict[str, Any]]:
        """Find an existing issue by topic (API, with error translation)."""
        try:
            return self._find_existing_issue_api(topic, None)
        except GitHubError as e:
            raise GitHubError(str(e))

    def _find_existing_issue_api(
        self, topic: str, date_str: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """Find existing issue via GitHub API.

        Args:
            topic: Topic name
            date_str: Date string

        Returns:
            Existing issue or None
        """
        self._rate_limit_check()

        search_term = f"[{self.config.issue_prefix}] {topic}"
        url = f"{self.base_url}/repos/{self.config.repository}/issues"
        params: Dict[str, Any] = {
            "state": "all",
            "per_page": 100,
        }

        try:
            response = requests.get(
                url, headers=self.headers, params=params, timeout=30
            )
            response.raise_for_status()
            issues = response.json()
            self._stats_api_calls += 1

            for issue in issues:
                if issue.get("title", "").startswith(f"[{self.config.issue_prefix}]"):
                    return issue
            return None
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            status = getattr(getattr(e, "response", None), "status_code", None)
            if status is None:
                # Fallback: try to parse status from exception string (e.g. "403 Forbidden")
                import re as _re
                m = _re.search(r"(\d{3})\s", str(e))
                if m:
                    status = int(m.group(1))
            if status == 403:
                raise GitHubError("GitHub API rate limit exceeded")
            if status == 401:
                raise GitHubError("GitHub authentication failed")
            raise GitHubError(f"404 Client Error: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to search for existing issues: {e}")

    def _create_issue(
        self, title: str, body: str
    ) -> Dict[str, Any]:
        """Create a GitHub issue.

        Args:
            title: Issue title
            body: Issue body content

        Returns:
            Created issue data

        Raises:
            GitHubError: If creation fails
        """
        if self.mock_mode:
            return self._mock_create_issue_inner(title, body)
        return self._create_issue_api_inner(title, body)

    def _mock_create_issue_inner(
        self, title: str, body: str
    ) -> Dict[str, Any]:
        """Create a mock issue.

        Args:
            title: Issue title
            body: Issue body content

        Returns:
            Mock issue data
        """
        issue_data = {
            "number": 2000 + (self._stats_issues_created % 1000),
            "title": title,
            "body": body,
            "state": "open",
            "html_url": f"https://github.com/{self.config.repository}/issues/{2000 + (self._stats_issues_created % 1000)}",
            "created": True,
            "mock_mode": True,
        }
        self._stats_issues_created += 1
        self.logger.debug(f"Mock issue created: {title}")
        return issue_data

    def _create_issue_api_inner(
        self, title: str, body: str
    ) -> Dict[str, Any]:
        """Create an issue via GitHub API.

        Args:
            title: Issue title
            body: Issue body content

        Returns:
            Created issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()

        url = f"{self.base_url}/repos/{self.config.repository}/issues"
        data: Dict[str, Any] = {
            "title": title,
            "body": body,
            "labels": [self.config.issue_prefix.lower(), "automated"],
        }

        try:
            response = requests.post(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            issue = response.json()
            self._stats_issues_created += 1
            self._stats_api_calls += 1
            self.logger.info(f"Created issue: #{issue['number']}")
            return issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to create issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to create issue: {e}")

    def _update_issue(
        self,
        target: Union[int, Dict[str, Any]],
        papers: List[ScoredPaper],
    ) -> Dict[str, Any]:
        """Update an existing GitHub issue.

        Args:
            target: Issue number or existing issue data dict
            papers: List of ScoredPaper objects

        Returns:
            Updated issue data

        Raises:
            GitHubError: If update fails
        """
        if self.mock_mode:
            return self._mock_update_issue(target, papers)
        return self._update_issue_api(target, papers)

    def _mock_update_issue(
        self, target: Union[int, Dict[str, Any]], papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Update a mock issue.

        Args:
            target: Issue number or existing issue data
            papers: List of papers

        Returns:
            Mock issue data
        """
        if isinstance(target, dict):
            issue = target
        else:
            existing = self.find_existing_issue_for_date(
                "", None
            )
            issue = existing if existing else {"number": target, "title": ""}
        updated_issue = issue.copy()
        updated_issue["body"] = self._format_issue_body(
            updated_issue.get("title", "").replace(
                self.config.issue_prefix + "] ", ""
            ),
            papers,
        )
        updated_issue["updated_at"] = datetime.now().isoformat()
        updated_issue["updated"] = True
        updated_issue["mock_mode"] = True
        self._stats_issues_updated += 1
        self.logger.debug(f"Mock issue updated: #{updated_issue.get('number')}")
        return updated_issue

    def _update_issue_api(
        self, target: Union[int, Dict[str, Any]], papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Update an issue via GitHub API.

        Args:
            target: Issue number or existing issue data
            papers: List of papers

        Returns:
            Updated issue data

        Raises:
            GitHubError: If API call fails
        """
        self._rate_limit_check()
        if isinstance(target, dict):
            issue_number = target.get("number")
            title = target.get("title", "")
        else:
            issue_number = target
            title = ""
        if not issue_number:
            raise GitHubError("Issue must have a number")
        topic_name = title.replace(f"[{self.config.issue_prefix}] ", "")
        body = self._format_issue_body(topic_name, papers)
        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}"
        data: Dict[str, Any] = {
            "title": title,
            "body": body,
        }
        try:
            response = requests.patch(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            updated_issue = response.json()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Updated issue: #{issue_number}")
            return updated_issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to update issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to update issue: {e}")

    def _update_issue_api_inner(
        self, issue: Dict[str, Any], papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Update an issue via GitHub API (legacy inner)."""
        self._rate_limit_check()

        issue_number = issue.get("number")
        if not issue_number:
            raise GitHubError("Issue must have a number")

        topic_name = issue.get("title", "").replace(
            f"[{self.config.issue_prefix}] ", ""
        )
        body = self._format_issue_body(topic_name, papers)

        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue_number}"
        data: Dict[str, Any] = {
            "title": issue.get("title", ""),
            "body": body,
        }

        try:
            response = requests.patch(
                url, headers=self.headers, json=data, timeout=30
            )
            response.raise_for_status()
            updated_issue = response.json()
            self._stats_issues_updated += 1
            self._stats_api_calls += 1
            self.logger.info(f"Updated issue: #{issue_number}")
            return updated_issue
        except HTTPError as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to update issue: {e}")
        except Exception as e:
            self._stats_errors += 1
            self._stats_api_calls += 1
            raise GitHubError(f"Failed to update issue: {e}")

    def _rate_limit_check(self) -> None:
        """Check and enforce rate limiting."""
        current_time = time.time()
        time_since_last = current_time - self._last_request_time

        if time_since_last < self._min_request_interval:
            sleep_time = self._min_request_interval - time_since_last
            if not self.mock_mode:
                time.sleep(sleep_time)
                self._stats_rate_limit_hits += 1
            self.logger.debug(f"Rate limited: sleeping {sleep_time:.2f}s")

        self._last_request_time = time.time()

        if self._rate_limit_remaining < 100 and self._rate_limit_remaining > 0:
            if self._rate_limit_reset > current_time:
                sleep_time = self._rate_limit_reset - current_time + 1
                if not self.mock_mode:
                    time.sleep(sleep_time)
                self.logger.warning(
                    f"Rate limit low ({self._rate_limit_remaining}), "
                    f"sleeping {sleep_time:.1f}s"
                )

    def get_statistics(self) -> Dict[str, Any]:
        """Get statistics for GitHub operations.

        Returns:
            Dictionary with statistics
        """
        return {
            "issues_created": self._stats_issues_created,
            "issues_updated": self._stats_issues_updated,
            "comments_added": self._stats_comments_added,
            "api_calls": self._stats_api_calls,
            "rate_limit_hits": self._stats_rate_limit_hits,
            "errors": self._stats_errors,
            "error_count": self._stats_errors,
            "last_activity": datetime.now().isoformat(),
            "mock_mode": self.mock_mode,
        }

    def get_statistics_summary(self) -> Dict[str, Any]:
        """Get a summary of statistics.

        Returns:
            Dictionary with summary statistics
        """
        return {
            "total_issues": self._stats_issues_created + self._stats_issues_updated,
            "issues_created": self._stats_issues_created,
            "issues_updated": self._stats_issues_updated,
            "comments_added": self._stats_comments_added,
            "mock_mode": self.mock_mode,
            "repository": self.config.repository if self.config else None,
        }

    def close_issue_by_number(self, issue_number: int) -> Dict[str, Any]:
        """Close an issue by its number.

        Args:
            issue_number: The GitHub issue number

        Returns:
            Updated issue data

        Raises:
            GitHubError: If closing fails
        """
        return self.close_issue(issue_number)

    def reopen_issue_by_number(self, issue_number: int) -> Dict[str, Any]:
        """Reopen an issue by its number.

        Args:
            issue_number: The GitHub issue number

        Returns:
            Updated issue data

        Raises:
            GitHubError: If reopening fails
        """
        return self.reopen_issue(issue_number)

    def retry_with_delay(
        self, func, *args, max_retries: int = 3, **kwargs
    ) -> Any:
        """Retry a function with exponential backoff.

        Args:
            func: Function to retry
            *args: Function arguments
            max_retries: Maximum number of retries
            **kwargs: Function keyword arguments

        Returns:
            Function result

        Raises:
            Exception: If all retries fail
        """
        for attempt in range(max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if attempt == max_retries - 1:
                    raise
                delay = 2**attempt
                self.logger.warning(
                    f"Attempt {attempt + 1}/{max_retries} failed: {e}. "
                    f"Retrying in {delay}s..."
                )
                time.sleep(delay)
