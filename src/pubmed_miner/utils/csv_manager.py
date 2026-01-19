"""
CSV management utilities for saving and loading paper collection data.
"""

import csv
import os
import logging
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime

from ..models import Paper, ScoredPaper

logger = logging.getLogger(__name__)


class CSVManager:
    """Manager for handling CSV operations on paper collections."""

    # CSV column headers
    HEADERS = [
        "pmid",
        "title",
        "authors",
        "journal",
        "publication_date",
        "doi",
        "abstract",
        "topic",
    ]

    SCORED_HEADERS = HEADERS + [
        "citation_count",
        "impact_factor",
        "score",
        "rank",
    ]

    @staticmethod
    def save_papers(
        papers: List[Paper],
        filepath: str,
        include_scoring: bool = False,
        append: bool = False,
    ) -> None:
        """Save papers to a CSV file.

        Args:
            papers: List of Paper or ScoredPaper objects to save
            filepath: Path to the CSV file
            include_scoring: If True, include scoring information for ScoredPaper objects
            append: If True, append to existing file; otherwise, overwrite

        Raises:
            ValueError: If filepath is invalid or papers list is empty
            IOError: If file operation fails
        """
        if not papers:
            raise ValueError("Papers list cannot be empty")

        filepath = Path(filepath)

        # Create parent directories if needed
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # Determine headers and check if we should include scoring
        headers = CSVManager.SCORED_HEADERS if include_scoring else CSVManager.HEADERS
        mode = "a" if append and filepath.exists() else "w"
        write_header = mode == "w" or not filepath.exists()

        try:
            with open(filepath, mode=mode, newline="", encoding="utf-8") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=headers)

                # Write header only if creating new file or file is empty
                if write_header:
                    writer.writeheader()

                # Write paper rows
                for paper in papers:
                    row = CSVManager._paper_to_dict(paper, include_scoring)
                    writer.writerow(row)

            file_mode = "appended to" if append and mode == "a" else "created"
            logger.info(
                f"Successfully {file_mode} CSV file: {filepath} "
                f"({len(papers)} papers saved)"
            )

        except IOError as e:
            logger.error(f"Error writing to CSV file {filepath}: {e}")
            raise IOError(f"Failed to save papers to CSV: {e}")
        except Exception as e:
            logger.error(f"Unexpected error saving papers to CSV: {e}")
            raise

    @staticmethod
    def load_papers(filepath: str) -> List[Dict]:
        """Load papers from a CSV file.

        Args:
            filepath: Path to the CSV file

        Returns:
            List of dictionaries containing paper data

        Raises:
            FileNotFoundError: If file does not exist
            IOError: If file operation fails
        """
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"CSV file not found: {filepath}")

        papers = []

        try:
            with open(filepath, mode="r", newline="", encoding="utf-8") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    papers.append(row)

            logger.info(f"Successfully loaded {len(papers)} papers from {filepath}")
            return papers

        except IOError as e:
            logger.error(f"Error reading CSV file {filepath}: {e}")
            raise IOError(f"Failed to load papers from CSV: {e}")
        except Exception as e:
            logger.error(f"Unexpected error loading papers from CSV: {e}")
            raise

    @staticmethod
    def _paper_to_dict(paper: Paper, include_scoring: bool = False) -> Dict:
        """Convert a Paper object to a dictionary for CSV writing.

        Args:
            paper: Paper or ScoredPaper object
            include_scoring: If True, include scoring fields for ScoredPaper

        Returns:
            Dictionary representation of the paper
        """
        # Convert authors list to comma-separated string
        authors_str = "; ".join(paper.authors) if paper.authors else ""

        # Format publication date
        pub_date_str = (
            paper.publication_date.isoformat()
            if isinstance(paper.publication_date, datetime)
            else str(paper.publication_date)
        )

        row = {
            "pmid": paper.pmid,
            "title": paper.title,
            "authors": authors_str,
            "journal": paper.journal,
            "publication_date": pub_date_str,
            "doi": paper.doi or "",
            "abstract": paper.abstract or "",
            "topic": paper.topic or "",
        }

        # Add scoring information if available and requested
        if include_scoring and isinstance(paper, ScoredPaper):
            row.update(
                {
                    "citation_count": paper.citation_count,
                    "impact_factor": paper.impact_factor,
                    "score": paper.score,
                    "rank": paper.rank,
                }
            )

        return row

    @staticmethod
    def append_papers(papers: List[Paper], filepath: str) -> None:
        """Append papers to an existing CSV file.

        Args:
            papers: List of Paper objects to append
            filepath: Path to the CSV file

        Raises:
            ValueError: If papers list is empty
            IOError: If file operation fails
        """
        if not papers:
            raise ValueError("Papers list cannot be empty")

        CSVManager.save_papers(papers, filepath, append=True)

    @staticmethod
    def update_collection(
        papers: List[ScoredPaper], filepath: str
    ) -> None:
        """Update collection with scored papers.

        Args:
            papers: List of ScoredPaper objects
            filepath: Path to the CSV file
        """
        CSVManager.save_papers(papers, filepath, include_scoring=True)
