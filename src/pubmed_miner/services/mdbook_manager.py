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
            book_root: Root directory of the mdbook (containing book_src/)
        """
        self.book_root = Path(book_root)
        self.src_dir = self.book_root / "book_src"
        self.summary_path = self.src_dir / "SUMMARY.md"
        
        # Ensure book_src directory exists
        self.src_dir.mkdir(parents=True, exist_ok=True)

    def create_daily_page(self, topic: str, papers: List[ScoredPaper]) -> str:
        """Create a markdown page for the day's papers.

        Args:
            topic: Topic name
            papers: List of essential papers

        Returns:
            Path to the created file relative to book_src/
        """
        current_date = datetime.now()
        year = current_date.strftime("%Y")
        month = current_date.strftime("%m")
        day = current_date.strftime("%d")
        
        # Create directory structure: book_src/year/month/
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

    def update_monthly_page(self, topic: str, papers: List[ScoredPaper], date: datetime = None) -> str:
        """Update a monthly markdown page with new papers in a table format.

        Args:
            topic: Topic name
            papers: List of essential papers
            date: Optional date to use for year/month. Defaults to current date.

        Returns:
            Path to the monthly file relative to book_src/
        """
        current_date = date or datetime.now()
        year = current_date.strftime("%Y")
        month = current_date.strftime("%m")
        month_name = current_date.strftime("%B")
        
        # Monthly file: book_src/year/month.md (e.g., book_src/2026/02.md)
        monthly_dir = self.src_dir / year
        monthly_dir.mkdir(parents=True, exist_ok=True)
        filename = f"{month}.md"
        file_path = monthly_dir / filename
        
        # 1. Read existing PMIDs to avoid duplicates
        import re
        existing_pmids = set()
        if file_path.exists():
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
                # Simple PMID extraction from [PMID](https://pubmed.ncbi.nlm.nih.gov/12345/)
                existing_pmids = set(re.findall(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/', content))

        # 2. Filter out duplicates
        new_papers = [p for p in papers if str(p.pmid) not in existing_pmids]
        
        if not new_papers:
            logger.info(f"No new papers to add for {topic} in {year}-{month}")
            return f"{year}/{filename}"

        # 3. Create or Update table
        if not file_path.exists():
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(f"# {year} {month_name}\n\n")
                f.write("| 날짜 | 주제 | 제목 | 저널 | 점수 | 링크 |\n")
                f.write("| :--- | :--- | :--- | :--- | :--- | :--- |\n")

        with open(file_path, "a", encoding="utf-8") as f:
            # Sort by rank to ensure most important ones are added first if multiple new ones
            for paper in sorted(new_papers, key=lambda p: p.rank):
                date_str = current_date.strftime("%Y-%m-%d")
                pmid_link = f"[PMID](https://pubmed.ncbi.nlm.nih.gov/{paper.pmid}/)"
                doi_link = f", [DOI](https://doi.org/{paper.doi})" if paper.doi else ""
                
                # Title might contain | character, which breaks MD table
                safe_title = paper.title.replace("|", "\\|")
                
                f.write(f"| {date_str} | {topic} | {safe_title} | {paper.journal} | {paper.score:.1f} | {pmid_link}{doi_link} |\n")

        logger.info(f"Updated monthly page: {file_path}")
        return f"{year}/{filename}"

    def update_summary(self, relative_path: str, topic: str):
        """Update SUMMARY.md to include the page.
        
        Args:
            relative_path: Path to the page relative to book_src/
            topic: Topic name
        """
        # We want to use monthly pages now
        is_monthly = relative_path.count('/') == 1 and relative_path.endswith('.md')
        
        if is_monthly:
            # Extract month from relative_path (e.g., "2026/02.md")
            parts = relative_path.split('/')
            year = parts[0]
            month_num = parts[1].replace('.md', '')
            try:
                dt = datetime.strptime(f"{year}-{month_num}-01", "%Y-%m-%d")
                month_name = dt.strftime("%B")
                link_text = f"{month_name} {year}"
            except ValueError:
                link_text = f"Month {month_num} {year}"
        else:
            current_date = datetime.now()
            date_str = current_date.strftime("%Y-%m-%d")
            link_text = f"{date_str}: {topic}"
            
        new_entry = f"    - [{link_text}]({relative_path})\n"
        
        if not self.summary_path.exists():
            with open(self.summary_path, "w", encoding="utf-8") as f:
                f.write("# 요약\n\n- [소개](README.md)\n\n# 업데이트\n")

        with open(self.summary_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        # Check if the entry already exists
        if any(relative_path in line for line in lines):
            logger.info(f"Entry for {relative_path} already exists in SUMMARY.md")
            return

        # Prepare year/month headers
        year_header = f"- [{current_date.year}]()"
        
        content = "".join(lines)
        
        # For monthly pages, we don't want topic-level nesting if possible, 
        # but let's stick to the structure.
        
        if year_header not in content:
            with open(self.summary_path, "a", encoding="utf-8") as f:
                f.write(f"\n{year_header}\n")
                f.write(f"{new_entry}")
        else:
            # Simple append is safer for now given the current logic
            with open(self.summary_path, "a", encoding="utf-8") as f:
                f.write(f"{new_entry}")

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
            f"총 논문 수: {len(papers)}",
            "",
            "## 중요도별 주요 논문",
            "",
        ]
        
        if not papers:
             lines.append("오늘 이 주제에 대한 필수 논문을 찾지 못했습니다.")
             return "\n".join(lines)

        # Sort papers by rank
        sorted_papers = sorted(papers, key=lambda p: p.rank)

        for paper in sorted_papers:
            # Format paper entry using HTML for better control (optional, but MD is safer for mdbook)
            # We'll use the custom CSS classes we defined
            
            lines.append(f'<div class="paper-entry">')
            lines.append(f'<div class="paper-title">{paper.rank}. {paper.title}</div>')
            
            # Authors and Journal info
            authors = ', '.join(paper.authors)
            lines.append(f'<div class="paper-meta">')
            lines.append(f"저자: {authors}<br>")
            lines.append(f"저널: {paper.journal}")
            lines.append('</div>')
            lines.append("")
            
            # Meta info table
            doi_link = f"[{paper.doi}](https://doi.org/{paper.doi})" if paper.doi else "-"
            
            lines.append("| 날짜 | PMID | DOI | 점수 (인용/IF) |")
            lines.append("|:---:|:---:|:---:|:---:|")
            lines.append(f"| {paper.publication_date.strftime('%Y-%m-%d')} | [{paper.pmid}](https://pubmed.ncbi.nlm.nih.gov/{paper.pmid}/) | {doi_link} | {paper.score:.1f} ({paper.citation_count}/{paper.impact_factor:.1f}) |")
            lines.append("")

            # Abstract
            if paper.abstract:
                lines.append(f'<div class="paper-abstract">')
                lines.append(f"<details><summary>초록</summary>")
                lines.append(f"<p>{paper.abstract}</p>")
                lines.append(f"</details>")
                lines.append('</div>')

            lines.append('</div>') # End paper-entry
            lines.append("")

        return "\n".join(lines)
