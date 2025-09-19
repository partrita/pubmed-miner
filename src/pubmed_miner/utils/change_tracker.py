"""
Change tracking utilities for monitoring paper list updates.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

from ..models import ScoredPaper

logger = logging.getLogger(__name__)


class ChangeTracker:
    """Tracks changes in paper lists between runs."""

    def __init__(self, cache_dir: str = "cache"):
        """Initialize change tracker.

        Args:
            cache_dir: Directory to store tracking data
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.tracking_file = self.cache_dir / "paper_tracking.json"

        logger.info(f"Initialized ChangeTracker with cache: {self.tracking_file}")

    def save_paper_snapshot(self, topic: str, papers: List[ScoredPaper]) -> None:
        """Save a snapshot of papers for a topic.

        Args:
            topic: Topic name
            papers: List of scored papers
        """
        try:
            # Load existing tracking data
            tracking_data = self._load_tracking_data()

            # Create snapshot
            snapshot = {
                "timestamp": datetime.now().isoformat(),
                "paper_count": len(papers),
                "papers": {},
            }

            for paper in papers:
                snapshot["papers"][paper.pmid] = {
                    "title": paper.title,
                    "authors": paper.authors,
                    "journal": paper.journal,
                    "publication_date": paper.publication_date.isoformat(),
                    "citation_count": paper.citation_count,
                    "impact_factor": paper.impact_factor,
                    "score": paper.score,
                    "rank": paper.rank,
                    "doi": paper.doi,
                }

            # Update tracking data
            tracking_data[topic] = snapshot

            # Save to file
            self._save_tracking_data(tracking_data)

            logger.info(f"Saved snapshot for topic '{topic}' with {len(papers)} papers")

        except Exception as e:
            logger.error(f"Error saving paper snapshot for {topic}: {e}")

    def get_previous_snapshot(self, topic: str) -> Optional[Dict]:
        """Get the previous snapshot for a topic.

        Args:
            topic: Topic name

        Returns:
            Previous snapshot data or None
        """
        try:
            tracking_data = self._load_tracking_data()
            return tracking_data.get(topic)
        except Exception as e:
            logger.error(f"Error loading previous snapshot for {topic}: {e}")
            return None

    def detect_changes(
        self, topic: str, current_papers: List[ScoredPaper]
    ) -> Dict[str, any]:
        """Detect changes between current and previous paper lists.

        Args:
            topic: Topic name
            current_papers: Current list of papers

        Returns:
            Dictionary with change information
        """
        previous_snapshot = self.get_previous_snapshot(topic)

        if not previous_snapshot:
            logger.info(f"No previous snapshot found for topic '{topic}'")
            return {
                "has_changes": True,
                "is_first_run": True,
                "summary": f"First run - tracking {len(current_papers)} papers",
                "changes": [],
                "new_papers": len(current_papers),
                "removed_papers": 0,
                "rank_changes": 0,
                "score_changes": 0,
            }

        # Convert current papers to dict for comparison
        current_dict = {p.pmid: p for p in current_papers}
        previous_dict = previous_snapshot.get("papers", {})

        # Detect different types of changes
        changes = self._analyze_changes(previous_dict, current_dict)

        # Determine if there are significant changes
        has_changes = (
            changes["new_papers"] > 0
            or changes["removed_papers"] > 0
            or changes["rank_changes"] > 0
            or changes["significant_score_changes"] > 0
        )

        change_summary = self._create_change_summary(changes)

        result = {
            "has_changes": has_changes,
            "is_first_run": False,
            "summary": change_summary,
            "changes": changes["detailed_changes"],
            "new_papers": changes["new_papers"],
            "removed_papers": changes["removed_papers"],
            "rank_changes": changes["rank_changes"],
            "score_changes": changes["significant_score_changes"],
            "previous_count": len(previous_dict),
            "current_count": len(current_dict),
            "previous_timestamp": previous_snapshot.get("timestamp"),
        }

        logger.info(f"Change detection for '{topic}': {change_summary}")
        return result

    def _analyze_changes(
        self, previous_dict: Dict[str, Dict], current_dict: Dict[str, ScoredPaper]
    ) -> Dict[str, any]:
        """Analyze changes between previous and current paper dictionaries.

        Args:
            previous_dict: Previous papers data
            current_dict: Current papers data

        Returns:
            Dictionary with detailed change analysis
        """
        changes = {
            "new_papers": 0,
            "removed_papers": 0,
            "rank_changes": 0,
            "significant_score_changes": 0,
            "detailed_changes": [],
        }

        previous_pmids = set(previous_dict.keys())
        current_pmids = set(current_dict.keys())

        # New papers
        new_pmids = current_pmids - previous_pmids
        changes["new_papers"] = len(new_pmids)

        for pmid in new_pmids:
            paper = current_dict[pmid]
            changes["detailed_changes"].append(
                {
                    "type": "new_paper",
                    "pmid": pmid,
                    "title": paper.title,
                    "rank": paper.rank,
                    "score": paper.score,
                    "description": f"New paper added at rank #{paper.rank}: {paper.title}",
                }
            )

        # Removed papers
        removed_pmids = previous_pmids - current_pmids
        changes["removed_papers"] = len(removed_pmids)

        for pmid in removed_pmids:
            prev_paper = previous_dict[pmid]
            changes["detailed_changes"].append(
                {
                    "type": "removed_paper",
                    "pmid": pmid,
                    "title": prev_paper["title"],
                    "previous_rank": prev_paper["rank"],
                    "description": f"Paper removed (was rank #{prev_paper['rank']}): {prev_paper['title']}",
                }
            )

        # Changes in existing papers
        common_pmids = previous_pmids & current_pmids

        for pmid in common_pmids:
            prev_paper = previous_dict[pmid]
            curr_paper = current_dict[pmid]

            # Rank changes
            prev_rank = prev_paper["rank"]
            curr_rank = curr_paper.rank

            if prev_rank != curr_rank:
                changes["rank_changes"] += 1
                direction = "improved" if curr_rank < prev_rank else "declined"
                changes["detailed_changes"].append(
                    {
                        "type": "rank_change",
                        "pmid": pmid,
                        "title": curr_paper.title,
                        "previous_rank": prev_rank,
                        "current_rank": curr_rank,
                        "direction": direction,
                        "description": f"Rank change: {curr_paper.title} moved from #{prev_rank} to #{curr_rank}",
                    }
                )

            # Significant score changes (>5 points)
            prev_score = prev_paper["score"]
            curr_score = curr_paper.score
            score_diff = abs(curr_score - prev_score)

            if score_diff > 5.0:
                changes["significant_score_changes"] += 1
                direction = "increased" if curr_score > prev_score else "decreased"
                changes["detailed_changes"].append(
                    {
                        "type": "score_change",
                        "pmid": pmid,
                        "title": curr_paper.title,
                        "previous_score": prev_score,
                        "current_score": curr_score,
                        "score_diff": score_diff,
                        "direction": direction,
                        "description": f"Score change: {curr_paper.title} score {direction} by {score_diff:.1f} points",
                    }
                )

            # Citation count changes
            prev_citations = prev_paper.get("citation_count", 0)
            curr_citations = curr_paper.citation_count

            if curr_citations != prev_citations:
                citation_diff = curr_citations - prev_citations
                if abs(citation_diff) >= 5:  # Only track significant citation changes
                    changes["detailed_changes"].append(
                        {
                            "type": "citation_change",
                            "pmid": pmid,
                            "title": curr_paper.title,
                            "previous_citations": prev_citations,
                            "current_citations": curr_citations,
                            "citation_diff": citation_diff,
                            "description": f"Citation update: {curr_paper.title} citations changed by {citation_diff:+d}",
                        }
                    )

        return changes

    def _create_change_summary(self, changes: Dict[str, any]) -> str:
        """Create a human-readable summary of changes.

        Args:
            changes: Changes analysis data

        Returns:
            Summary string
        """
        summary_parts = []

        if changes["new_papers"] > 0:
            summary_parts.append(f"{changes['new_papers']} new papers")

        if changes["removed_papers"] > 0:
            summary_parts.append(f"{changes['removed_papers']} removed papers")

        if changes["rank_changes"] > 0:
            summary_parts.append(f"{changes['rank_changes']} rank changes")

        if changes["significant_score_changes"] > 0:
            summary_parts.append(
                f"{changes['significant_score_changes']} significant score changes"
            )

        if not summary_parts:
            return "No significant changes detected"

        return ", ".join(summary_parts)

    def get_topic_history(self, topic: str, limit: int = 10) -> List[Dict]:
        """Get historical snapshots for a topic.

        Args:
            topic: Topic name
            limit: Maximum number of historical entries

        Returns:
            List of historical snapshot data
        """
        try:
            history_file = self.cache_dir / f"history_{topic}.json"

            if not history_file.exists():
                return []

            with open(history_file, "r", encoding="utf-8") as f:
                history = json.load(f)

            # Return most recent entries first
            return history[-limit:] if len(history) > limit else history

        except Exception as e:
            logger.error(f"Error loading topic history for {topic}: {e}")
            return []

    def save_to_history(self, topic: str, papers: List[ScoredPaper]) -> None:
        """Save current papers to historical record.

        Args:
            topic: Topic name
            papers: Current papers list
        """
        try:
            history_file = self.cache_dir / f"history_{topic}.json"

            # Load existing history
            history = []
            if history_file.exists():
                with open(history_file, "r", encoding="utf-8") as f:
                    history = json.load(f)

            # Create new entry
            entry = {
                "timestamp": datetime.now().isoformat(),
                "paper_count": len(papers),
                "top_papers": [],
            }

            # Save top 5 papers for quick reference
            top_papers = sorted(papers, key=lambda p: p.rank)[:5]
            for paper in top_papers:
                entry["top_papers"].append(
                    {
                        "pmid": paper.pmid,
                        "title": paper.title,
                        "rank": paper.rank,
                        "score": paper.score,
                    }
                )

            history.append(entry)

            # Keep only last 50 entries
            if len(history) > 50:
                history = history[-50:]

            # Save updated history
            with open(history_file, "w", encoding="utf-8") as f:
                json.dump(history, f, indent=2, ensure_ascii=False)

            logger.debug(f"Saved history entry for topic '{topic}'")

        except Exception as e:
            logger.error(f"Error saving topic history for {topic}: {e}")

    def cleanup_old_data(self, days_to_keep: int = 30) -> int:
        """Clean up old tracking data.

        Args:
            days_to_keep: Number of days of data to keep

        Returns:
            Number of entries cleaned up
        """
        try:
            cutoff_date = datetime.now().timestamp() - (days_to_keep * 24 * 3600)
            cleaned_count = 0

            # Clean up main tracking file
            tracking_data = self._load_tracking_data()
            topics_to_remove = []

            for topic, snapshot in tracking_data.items():
                try:
                    snapshot_time = datetime.fromisoformat(
                        snapshot["timestamp"]
                    ).timestamp()
                    if snapshot_time < cutoff_date:
                        topics_to_remove.append(topic)
                except (ValueError, KeyError):
                    # Remove invalid entries
                    topics_to_remove.append(topic)

            for topic in topics_to_remove:
                del tracking_data[topic]
                cleaned_count += 1

            if topics_to_remove:
                self._save_tracking_data(tracking_data)

            # Clean up history files
            for history_file in self.cache_dir.glob("history_*.json"):
                try:
                    with open(history_file, "r", encoding="utf-8") as f:
                        history = json.load(f)

                    # Filter out old entries
                    filtered_history = []
                    for entry in history:
                        try:
                            entry_time = datetime.fromisoformat(
                                entry["timestamp"]
                            ).timestamp()
                            if entry_time >= cutoff_date:
                                filtered_history.append(entry)
                            else:
                                cleaned_count += 1
                        except (ValueError, KeyError):
                            cleaned_count += 1

                    # Save filtered history or remove file if empty
                    if filtered_history:
                        with open(history_file, "w", encoding="utf-8") as f:
                            json.dump(filtered_history, f, indent=2, ensure_ascii=False)
                    else:
                        history_file.unlink()

                except Exception as e:
                    logger.warning(f"Error cleaning history file {history_file}: {e}")

            logger.info(f"Cleaned up {cleaned_count} old tracking entries")
            return cleaned_count

        except Exception as e:
            logger.error(f"Error during cleanup: {e}")
            return 0

    def _load_tracking_data(self) -> Dict:
        """Load tracking data from file.

        Returns:
            Tracking data dictionary
        """
        if not self.tracking_file.exists():
            return {}

        try:
            with open(self.tracking_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Error loading tracking data: {e}")
            return {}

    def _save_tracking_data(self, data: Dict) -> None:
        """Save tracking data to file.

        Args:
            data: Tracking data to save
        """
        try:
            with open(self.tracking_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except IOError as e:
            logger.error(f"Error saving tracking data: {e}")

    def get_statistics(self) -> Dict[str, any]:
        """Get tracking statistics.

        Returns:
            Dictionary with tracking statistics
        """
        try:
            tracking_data = self._load_tracking_data()

            stats = {
                "tracked_topics": len(tracking_data),
                "total_papers_tracked": 0,
                "oldest_snapshot": None,
                "newest_snapshot": None,
                "topics": [],
            }

            timestamps = []

            for topic, snapshot in tracking_data.items():
                paper_count = snapshot.get("paper_count", 0)
                timestamp = snapshot.get("timestamp")

                stats["total_papers_tracked"] += paper_count
                stats["topics"].append(
                    {
                        "name": topic,
                        "paper_count": paper_count,
                        "last_updated": timestamp,
                    }
                )

                if timestamp:
                    timestamps.append(timestamp)

            if timestamps:
                timestamps.sort()
                stats["oldest_snapshot"] = timestamps[0]
                stats["newest_snapshot"] = timestamps[-1]

            return stats

        except Exception as e:
            logger.error(f"Error getting tracking statistics: {e}")
            return {"tracked_topics": 0, "total_papers_tracked": 0}
