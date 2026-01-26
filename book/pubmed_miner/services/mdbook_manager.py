"""
MdBook management service.
"""

import logging
import os
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any

from ..models import ScoredPaper

logger = logging.getLogger(__name__)


class MdBookManager:
    """Manages mdbook content generation for essential papers."""

    def __init__(self, book_root: str = "."):
        """Initialize MdBook manager.

        Args:
            book_root: Root directory of the mdbook (containing src/)
        """
        self.book_root = Path(book_root)
        self.src_dir = self.book_root / "src"
        self.summary_path = self.src_dir / "SUMMARY.md"
        
        # Ensure src directory exists
        self.src_dir.mkdir(parents=True, exist_ok=True)

    def create_daily_page(self, topic: str, papers: List[ScoredPaper]) -> str:
        """Create a markdown page for the day's papers.

        Args:
            topic: Topic name
            papers: List of essential papers

        Returns:
            Path to the created file relative to src/
        """
        current_date = datetime.now()
        year = current_date.strftime("%Y")
        month = current_date.strftime("%m")
        day = current_date.strftime("%d")
        
        # Create directory structure: src/year/month/
        topic_slug = topic.lower().replace(" ", "_")
        daily_dir = self.src_dir / year / month
        daily_dir.mkdir(parents=True, exist_ok=True)
        
        # File name: day_topic.md
        filename = f"{day}_{topic_slug}.md"
        file_path = daily_dir / filename
        
        content = self._format_page_content(topic, papers, current_date)
        
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
            
        logger.info(f"Created daily page: {file_path}")
        
        # Return relative path for SUMMARY.md
        return f"{year}/{month}/{filename}"

    def update_summary(self, relative_path: str, topic: str):
        """Update SUMMARY.md to include the new page.
        
        Args:
            relative_path: Path to the new page relative to src/
            topic: Topic name
        """
        current_date = datetime.now()
        date_str = current_date.strftime("%Y-%m-%d")
        link_text = f"{date_str}: {topic}"
        new_entry = f"    - [{link_text}]({relative_path})\n"
        
        # Check if year section exists, if not add it
        year_section = f"- [{current_date.year}]()" # Dummy link or just text
        
        # Simple appending strategy for now. 
        # Ideally, we'd want to insert it in the correct chronological order.
        # Let's read existing content first.
        
        if not self.summary_path.exists():
            with open(self.summary_path, "w") as f:
                f.write("# Summary\n\n- [Introduction](README.md)\n\n# Daily Updates\n")

        with open(self.summary_path, "r") as f:
            lines = f.readlines()

        # We want to add the new entry under "Daily Updates"
        # Let's structure it by Year -> Month -> Day/Topic
        
        # However, to keep it simple and requested "daily page addition":
        # We will append to the end, or insert at top of a list.
        # Let's try to group by Year/Month if possible, or just a flat list under "Daily Updates".
        
        # Let's just append for now, but formatted nicely.
        
        # Check if the entry already exists
        if any(relative_path in line for line in lines):
            logger.info(f"Entry for {relative_path} already exists in SUMMARY.md")
            return

        # Prepare year/month headers
        year_header = f"- {current_date.year}"
        month_name = current_date.strftime("%B")
        month_header = f"  - {month_name}"
        
        # We need to construct the file content carefully.
        # This is a simple parser/updater.
        
        new_lines = []
        in_daily_updates = False
        year_found = False
        month_found = False
        
        # We will reconstruct the file.
        # This is complex to do robustly with just appending.
        # Let's just append to the end for simplicity, user can reorganize if needed.
        # Or better: simple text search.
        
        content = "".join(lines)
        
        # Ensure we have the headers
        if year_header not in content:
            # Append Year header
            with open(self.summary_path, "a") as f:
                f.write(f"\n{year_header}\n")
                f.write(f"{month_header}\n")
                f.write(f"    - [{link_text}]({relative_path})\n")
        elif month_header not in content:
            # Year exists, but month doesn't (simplified check)
             with open(self.summary_path, "a") as f:
                f.write(f"{month_header}\n")
                f.write(f"    - [{link_text}]({relative_path})\n")
        else:
            # Both exist, just append item
             with open(self.summary_path, "a") as f:
                f.write(f"    - [{link_text}]({relative_path})\n")

        logger.info("Updated SUMMARY.md")

    def _format_page_content(self, topic: str, papers: List[ScoredPaper], date: datetime) -> str:
        """Format papers list as markdown page.

        Args:
            topic: Topic name
            papers: List of scored papers
            date: Date object

        Returns:
            Markdown formatted content
        """
        date_str = date.strftime('%Y-%m-%d')
        
        lines = [
            f"# {date_str}: {topic}",
            "",
            f"**Total Papers:** {len(papers)}",
            "",
            "## Top Papers by Importance Score",
            "",
        ]
        
        if not papers:
             lines.append("No essential papers found for this topic today.")
             return "\n".join(lines)

        # Sort papers by rank
        sorted_papers = sorted(papers, key=lambda p: p.rank)

        for paper in sorted_papers:
            # Format paper entry using HTML for better control (optional, but MD is safer for mdbook)
            # We'll use the custom CSS classes we defined
            
            lines.append(f'<div class="paper-entry">')
            lines.append(f'<div class="paper-title">{paper.rank}. {paper.title}</div>')
            
            # Meta info
            authors = ', '.join(paper.authors)
            meta = [
                f"**Authors:** {authors}",
                f"**Journal:** {paper.journal}",
                f"**Date:** {paper.publication_date.strftime('%Y-%m-%d')}",
                f"**PMID:** [{paper.pmid}](https://pubmed.ncbi.nlm.nih.gov/{paper.pmid}/)",
            ]
            
            if paper.doi:
                meta.append(f"**DOI:** [{paper.doi}](https://doi.org/{paper.doi})")
                
            meta.append(f"**Score:** {paper.score:.1f} (Citations: {paper.citation_count}, IF: {paper.impact_factor:.1f})")
            
            lines.append(f'<div class="paper-meta">')
            lines.append(" | ".join(meta))
            lines.append('</div>')

            # Abstract
            if paper.abstract:
                lines.append(f'<div class="paper-abstract">')
                lines.append(f"<details><summary>Abstract</summary>")
                lines.append(f"<p>{paper.abstract}</p>")
                lines.append(f"</details>")
                lines.append('</div>')

            lines.append('</div>') # End paper-entry
            lines.append("")

        return "\n".join(lines)
