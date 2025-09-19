"""
GitHub Issues management service.
"""

import requests
import json
import logging
from typing import List, Dict, Optional, Any
from datetime import datetime

from ..models import ScoredPaper, GitHubConfig
from ..utils.error_handler import APIError, retry_api_calls

logger = logging.getLogger(__name__)


class GitHubIssuesManager:
    """Manages GitHub Issues for essential papers tracking."""

    def __init__(self, config: GitHubConfig):
        """Initialize GitHub Issues manager.

        Args:
            config: GitHub configuration
        """
        self.config = config
        self.base_url = "https://api.github.com"

        # Check if we're in mock mode (for local testing without token)
        self.mock_mode = (
            config.token == "mock_token_for_local_testing"
            or not config.token
            or config.token.startswith("mock_")
        )

        if self.mock_mode:
            logger.info(
                "GitHubIssuesManager initialized in MOCK MODE for local testing"
            )
            self.headers = {
                "Accept": "application/vnd.github.v3+json",
                "User-Agent": "PubMedMiner/1.0",
            }
        else:
            self.headers = {
                "Authorization": f"token {config.token}",
                "Accept": "application/vnd.github.v3+json",
                "User-Agent": "PubMedMiner/1.0",
            }
            logger.info(f"Initialized GitHubIssuesManager for {config.repository}")

    def create_or_update_issue(
        self, topic: str, papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Create or update an issue for a topic's essential papers.

        Args:
            topic: Topic name
            papers: List of essential papers

        Returns:
            Dictionary with issue information

        Raises:
            APIError: If GitHub API call fails
        """
        if self.mock_mode:
            return self._mock_create_or_update_issue(topic, papers)

        issue_title = f"[Essential Papers] {topic}"

        # Check if issue already exists
        existing_issue = self.find_existing_issue(topic)

        if existing_issue:
            logger.info(
                f"Updating existing issue #{existing_issue['number']} for topic: {topic}"
            )
            return self._update_issue(existing_issue, papers)
        else:
            logger.info(f"Creating new issue for topic: {topic}")
            return self._create_issue(issue_title, papers)

    def find_existing_issue(self, topic: str) -> Optional[Dict[str, Any]]:
        """Find existing issue for a topic.

        Args:
            topic: Topic name

        Returns:
            Issue dictionary or None if not found
        """
        if self.mock_mode:
            return self._mock_find_existing_issue(topic)

        try:
            search_query = (
                f'repo:{self.config.repository} is:issue "[Essential Papers] {topic}"'
            )
            url = f"{self.base_url}/search/issues"
            params = {"q": search_query}

            response = requests.get(url, headers=self.headers, params=params)
            response.raise_for_status()

            data = response.json()
            issues = data.get("items", [])

            if issues:
                # Return the first matching issue
                issue = issues[0]
                logger.debug(
                    f"Found existing issue #{issue['number']} for topic: {topic}"
                )
                return issue

            return None

        except requests.exceptions.RequestException as e:
            logger.error(f"Error searching for existing issue: {e}")
            return None

    def format_paper_list(self, papers: List[ScoredPaper]) -> str:
        """Format papers list as markdown for issue body.

        Args:
            papers: List of scored papers

        Returns:
            Markdown formatted paper list
        """
        if not papers:
            return "No essential papers found for this topic."

        # Sort papers by rank
        sorted_papers = sorted(papers, key=lambda p: p.rank)

        # Create markdown content
        lines = [
            "# Essential Papers",
            "",
            f"**Last Updated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"**Total Papers:** {len(papers)}",
            "",
            "## Top Papers by Importance Score",
            "",
        ]

        for paper in sorted_papers:
            # Format paper entry
            paper_lines = [
                f"### {paper.rank}. {paper.title}",
                "",
                f"**Authors:** {', '.join(paper.authors)}",
                f"**Journal:** {paper.journal}",
                f"**Publication Date:** {paper.publication_date.strftime('%Y-%m-%d')}",
                f"**PMID:** [{paper.pmid}](https://pubmed.ncbi.nlm.nih.gov/{paper.pmid}/)",
                f"**Citations:** {paper.citation_count}",
                f"**Impact Factor:** {paper.impact_factor:.2f}",
                f"**Importance Score:** {paper.score:.1f}/100",
                "",
            ]

            # Add DOI if available
            if paper.doi:
                paper_lines.insert(
                    -1, f"**DOI:** [{paper.doi}](https://doi.org/{paper.doi})"
                )

            # Add abstract if available (truncated)
            if paper.abstract:
                abstract = (
                    paper.abstract[:300] + "..."
                    if len(paper.abstract) > 300
                    else paper.abstract
                )
                paper_lines.extend([f"**Abstract:** {abstract}", ""])

            lines.extend(paper_lines)
            lines.append("---")
            lines.append("")

        # Add footer
        lines.extend(
            [
                "## About This List",
                "",
                "This list is automatically generated based on:",
                "- **Citation Count** (40%): Number of times the paper has been cited",
                "- **Journal Impact Factor** (30%): Impact factor of the publishing journal",
                "- **Publication Recency** (20%): How recently the paper was published",
                "- **Query Relevance** (10%): How well the paper matches the search query",
                "",
                "Papers are ranked by their composite importance score and updated daily.",
                "",
                f"*Generated by PubMed Miner - {datetime.now().strftime('%Y-%m-%d')}*",
            ]
        )

        return "\n".join(lines)

    @retry_api_calls(max_attempts=3, delay=1.0)
    def _create_issue(self, title: str, papers: List[ScoredPaper]) -> Dict[str, Any]:
        """Create a new GitHub issue.

        Args:
            title: Issue title
            papers: List of papers for issue body

        Returns:
            Created issue data

        Raises:
            APIError: If issue creation fails
        """
        url = f"{self.base_url}/repos/{self.config.repository}/issues"

        body = self.format_paper_list(papers)

        issue_data = {"title": title, "body": body, "labels": self.config.issue_labels}

        try:
            response = requests.post(
                url, headers=self.headers, data=json.dumps(issue_data)
            )
            response.raise_for_status()

            issue = response.json()
            logger.info(f"Created issue #{issue['number']}: {title}")

            return issue

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to create issue: {e}")
            raise APIError(f"GitHub issue creation failed: {e}")

    @retry_api_calls(max_attempts=3, delay=1.0)
    def _update_issue(
        self, issue: Dict[str, Any], papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Update an existing GitHub issue.

        Args:
            issue: Existing issue data
            papers: Updated list of papers

        Returns:
            Updated issue data

        Raises:
            APIError: If issue update fails
        """
        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue['number']}"

        body = self.format_paper_list(papers)

        update_data = {"body": body}

        try:
            response = requests.patch(
                url, headers=self.headers, data=json.dumps(update_data)
            )
            response.raise_for_status()

            updated_issue = response.json()
            logger.info(f"Updated issue #{updated_issue['number']}")

            return updated_issue

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to update issue #{issue['number']}: {e}")
            raise APIError(f"GitHub issue update failed: {e}")

    def add_comment_for_changes(
        self, issue: Dict[str, Any], changes: List[str]
    ) -> None:
        """Add a comment to an issue describing changes.

        Args:
            issue: Issue data
            changes: List of change descriptions
        """
        if not changes:
            return

        comment_body = "## 📊 Paper List Updates\n\n"
        comment_body += (
            f"**Updated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n\n"
        )
        comment_body += "### Changes:\n"

        for change in changes:
            comment_body += f"- {change}\n"

        comment_body += "\n*Automated update by PubMed Miner*"

        self._add_comment(issue, comment_body)

    @retry_api_calls(max_attempts=3, delay=1.0)
    def _add_comment(self, issue: Dict[str, Any], body: str) -> Dict[str, Any]:
        """Add a comment to an issue.

        Args:
            issue: Issue data
            body: Comment body

        Returns:
            Created comment data

        Raises:
            APIError: If comment creation fails
        """
        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue['number']}/comments"

        comment_data = {"body": body}

        try:
            response = requests.post(
                url, headers=self.headers, data=json.dumps(comment_data)
            )
            response.raise_for_status()

            comment = response.json()
            logger.info(f"Added comment to issue #{issue['number']}")

            return comment

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to add comment to issue #{issue['number']}: {e}")
            raise APIError(f"GitHub comment creation failed: {e}")

    def close_issue(self, topic: str) -> bool:
        """Close an issue for a topic.

        Args:
            topic: Topic name

        Returns:
            True if issue was closed, False if not found
        """
        existing_issue = self.find_existing_issue(topic)

        if not existing_issue:
            logger.warning(f"No issue found to close for topic: {topic}")
            return False

        return self._close_issue(existing_issue)

    @retry_api_calls(max_attempts=3, delay=1.0)
    def _close_issue(self, issue: Dict[str, Any]) -> bool:
        """Close a GitHub issue.

        Args:
            issue: Issue data

        Returns:
            True if successful

        Raises:
            APIError: If issue closing fails
        """
        url = f"{self.base_url}/repos/{self.config.repository}/issues/{issue['number']}"

        update_data = {"state": "closed"}

        try:
            response = requests.patch(
                url, headers=self.headers, data=json.dumps(update_data)
            )
            response.raise_for_status()

            logger.info(f"Closed issue #{issue['number']}")
            return True

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to close issue #{issue['number']}: {e}")
            raise APIError(f"GitHub issue closing failed: {e}")

    def get_repository_info(self) -> Dict[str, Any]:
        """Get repository information.

        Returns:
            Repository information

        Raises:
            APIError: If repository access fails
        """
        url = f"{self.base_url}/repos/{self.config.repository}"

        try:
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()

            repo_info = response.json()
            logger.debug(f"Retrieved repository info for {self.config.repository}")

            return repo_info

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to get repository info: {e}")
            raise APIError(f"GitHub repository access failed: {e}")

    def list_essential_paper_issues(self) -> List[Dict[str, Any]]:
        """List all essential paper issues in the repository.

        Returns:
            List of issue data
        """
        try:
            search_query = (
                f'repo:{self.config.repository} is:issue "[Essential Papers]"'
            )
            url = f"{self.base_url}/search/issues"
            params = {"q": search_query, "sort": "updated", "order": "desc"}

            response = requests.get(url, headers=self.headers, params=params)
            response.raise_for_status()

            data = response.json()
            issues = data.get("items", [])

            logger.info(f"Found {len(issues)} essential paper issues")
            return issues

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to list essential paper issues: {e}")
            return []

    def validate_access(self) -> bool:
        """Validate GitHub API access and repository permissions.

        Returns:
            True if access is valid (always True in mock mode)
        """
        if self.mock_mode:
            logger.info("[MOCK MODE] GitHub access validation - returning True")
            return True

        try:
            # Test repository access
            repo_info = self.get_repository_info()

            # Check if we have issues permission
            permissions = repo_info.get("permissions", {})
            can_create_issues = permissions.get("push", False) or permissions.get(
                "admin", False
            )

            if not can_create_issues:
                logger.error("Insufficient permissions to create issues")
                return False

            logger.info("GitHub API access validated successfully")
            return True

        except Exception as e:
            logger.error(f"GitHub API access validation failed: {e}")
            return False

    def get_rate_limit_info(self) -> Dict[str, Any]:
        """Get GitHub API rate limit information.

        Returns:
            Rate limit information
        """
        try:
            url = f"{self.base_url}/rate_limit"
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()

            rate_limit_info = response.json()
            logger.debug("Retrieved GitHub API rate limit info")

            return rate_limit_info

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to get rate limit info: {e}")
            return {}

    def detect_paper_changes(
        self, old_papers: List[ScoredPaper], new_papers: List[ScoredPaper]
    ) -> List[str]:
        """Detect changes between old and new paper lists.

        Args:
            old_papers: Previous paper list
            new_papers: New paper list

        Returns:
            List of change descriptions
        """
        changes = []

        # Create lookup dictionaries
        old_dict = {p.pmid: p for p in old_papers}
        new_dict = {p.pmid: p for p in new_papers}

        # Find new papers
        new_pmids = set(new_dict.keys()) - set(old_dict.keys())
        if new_pmids:
            changes.append(f"Added {len(new_pmids)} new papers")
            for pmid in list(new_pmids)[:3]:  # Show first 3
                paper = new_dict[pmid]
                changes.append(f"  ➕ {paper.title} (Rank #{paper.rank})")
            if len(new_pmids) > 3:
                changes.append(f"  ... and {len(new_pmids) - 3} more")

        # Find removed papers
        removed_pmids = set(old_dict.keys()) - set(new_dict.keys())
        if removed_pmids:
            changes.append(f"Removed {len(removed_pmids)} papers")

        # Find rank changes
        rank_changes = []
        for pmid in set(old_dict.keys()) & set(new_dict.keys()):
            old_rank = old_dict[pmid].rank
            new_rank = new_dict[pmid].rank
            if old_rank != new_rank:
                paper = new_dict[pmid]
                direction = "↗️" if new_rank < old_rank else "↘️"
                rank_changes.append(
                    f"  {direction} {paper.title}: #{old_rank} → #{new_rank}"
                )

        if rank_changes:
            changes.append(f"Rank changes for {len(rank_changes)} papers")
            changes.extend(rank_changes[:5])  # Show first 5 changes
            if len(rank_changes) > 5:
                changes.append(f"  ... and {len(rank_changes) - 5} more rank changes")

        return changes

    # Mock methods for local testing without GitHub token
    def _mock_create_or_update_issue(
        self, topic: str, papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Mock implementation for creating/updating issues."""
        logger.info(f"[MOCK] Creating/updating issue for topic: {topic}")
        logger.info(f"[MOCK] Would process {len(papers)} papers")

        # Simulate issue data
        mock_issue = {
            "number": hash(topic) % 10000,  # Generate consistent fake issue number
            "title": f"[Essential Papers] {topic}",
            "body": self.format_paper_list(papers),
            "state": "open",
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat(),
            "html_url": f"https://github.com/{self.config.repository}/issues/{hash(topic) % 10000}",
            "user": {"login": "mock_user"},
        }

        # Log some paper details for verification
        if papers:
            logger.info(
                f"[MOCK] Top paper: {papers[0].title} (Score: {papers[0].score:.1f})"
            )

        return mock_issue

    def _mock_find_existing_issue(self, topic: str) -> Optional[Dict[str, Any]]:
        """Mock implementation for finding existing issues."""
        logger.info(f"[MOCK] Searching for existing issue for topic: {topic}")

        # Simulate that no existing issue is found (always create new)
        # In real implementation, this would search GitHub
        return None

    def _mock_add_comment(self, issue: Dict[str, Any], body: str) -> Dict[str, Any]:
        """Mock implementation for adding comments."""
        logger.info(f"[MOCK] Adding comment to issue #{issue['number']}")
        logger.debug(f"[MOCK] Comment body: {body[:100]}...")

        return {
            "id": hash(body) % 100000,
            "body": body,
            "created_at": datetime.now().isoformat(),
            "user": {"login": "mock_user"},
        }

    def _mock_create_or_update_issue(
        self, topic: str, papers: List[ScoredPaper]
    ) -> Dict[str, Any]:
        """Mock version of create_or_update_issue for local testing.

        Args:
            topic: Topic name
            papers: List of essential papers

        Returns:
            Mock issue information
        """
        import random

        issue_number = random.randint(1, 1000)

        logger.info(f"[MOCK MODE] Would create/update issue for topic: {topic}")
        logger.info(f"[MOCK MODE] Issue would contain {len(papers)} papers")

        # Log paper details for verification
        for i, paper in enumerate(papers[:3]):  # Show first 3 papers
            logger.info(
                f"[MOCK MODE] Paper {i + 1}: {paper.title} (Score: {paper.score:.1f})"
            )

        if len(papers) > 3:
            logger.info(f"[MOCK MODE] ... and {len(papers) - 3} more papers")

        return {
            "number": issue_number,
            "title": f"[Essential Papers] {topic}",
            "html_url": f"https://github.com/{self.config.repository}/issues/{issue_number}",
            "state": "open",
            "created": True,
            "mock_mode": True,
            "body": self.format_paper_list(papers),
        }

    def _mock_find_existing_issue(self, topic: str) -> Optional[Dict[str, Any]]:
        """Mock version of find_existing_issue for local testing.

        Args:
            topic: Topic name

        Returns:
            None (always assume no existing issue in mock mode)
        """
        logger.debug(f"[MOCK MODE] Searching for existing issue for topic: {topic}")
        return None  # Always assume no existing issue for simplicity
