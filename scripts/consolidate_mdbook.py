#!/usr/bin/env python3
"""
Consolidate daily mdBook pages into monthly tables.
This script migrates the old daily-page structure to the new monthly-table structure.
"""

import os
import re
import sys
from pathlib import Path
from datetime import datetime
import logging

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pubmed_miner.services.mdbook_manager import MdBookManager
from pubmed_miner.models import ScoredPaper

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def parse_daily_page(file_path: Path):
    """Parse a daily page to extract papers."""
    papers = []
    topic = ""
    date_str = ""
    
    # Try to extract topic and date from filename or content
    # Filename: 26_ai-drug-discovery.md
    # Path: book_src/2026/01/26_ai-drug-discovery.md
    
    match = re.search(r'(\d{4})/(\d{2})/(\d{2})_(.+)\.md$', str(file_path))
    if match:
        year, month, day, topic_slug = match.groups()
        date_str = f"{year}-{month}-{day}"
        topic = topic_slug.replace('_', ' ').title()
    
    if not file_path.exists():
        return papers, topic, date_str

    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
        
    # Extract papers using regex
    # Each paper is in a div with class "paper-entry"
    # We need title, journal, pmid, doi, score
    
    # This is a bit complex with regex, but let's try to find PMIDs and basic info
    # <div class="paper-title">1. Title</div>
    # 저널: Journal
    # | Date | PMID | DOI | 점수 (인용/IF) |
    # | ... | [PMID](.../12345/) | [DOI](.../12345) | 95.8 (...) |
    
    entries = content.split('<div class="paper-entry">')[1:]
    for entry in entries:
        try:
            title_match = re.search(r'<div class="paper-title">\d+\. (.+?)</div>', entry)
            journal_match = re.search(r'저널: (.+?)<br>', entry) or re.search(r'저널: (.+?)\n', entry)
            pmid_match = re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/', entry)
            doi_match = re.search(r'doi\.org/(10\.\d+/[\w\.\-/]+)', entry)
            score_match = re.search(r'\| [\d-]+ \| \[?\d+\]?\(.*?\) \| .*? \| ([\d\.]+) \(', entry)
            
            if title_match and pmid_match:
                paper = ScoredPaper(
                    pmid=pmid_match.group(1),
                    title=title_match.group(1),
                    authors=["Unknown"], # Authors list cannot be empty per model validation
                    journal=journal_match.group(1).strip() if journal_match else "Unknown",
                    publication_date=datetime.strptime(date_str, "%Y-%m-%d") if date_str else datetime.now(),
                    abstract="",
                    doi=doi_match.group(1) if doi_match else None,
                    topic=topic,
                    citation_count=0,
                    impact_factor=0.0,
                    score=float(score_match.group(1)) if score_match else 0.0,
                    rank=0
                )
                papers.append(paper)
        except Exception as e:
            logger.warning(f"Failed to parse an entry in {file_path}: {e}")
            
    return papers, topic, date_str

def main():
    book_src = Path("book_src")
    if not book_src.exists():
        logger.error("book_src directory not found")
        return

    mdbook_manager = MdBookManager()
    
    # Find all daily md files
    # Structure: book_src/YYYY/MM/DD_topic.md
    daily_files = sorted(list(book_src.glob("202[4-9]/[0-1][0-9]/[0-3][0-9]_*.md")))
    
    logger.info(f"Found {len(daily_files)} daily files to consolidate")
    
    consolidated_count = 0
    for daily_file in daily_files:
        papers, topic, date_str = parse_daily_page(daily_file)
        if papers:
            paper_date = papers[0].publication_date
            relative_path = mdbook_manager.update_monthly_page(topic, papers, date=paper_date)
            mdbook_manager.update_summary(relative_path, topic)
            consolidated_count += 1
            logger.info(f"Consolidated {daily_file} -> {relative_path}")
            
    logger.info(f"Successfully consolidated {consolidated_count} daily pages.")
    print("\nNext steps:")
    print("1. Review the new monthly pages in book_src/")
    print("2. You can manually delete the old daily directories (e.g., book_src/2026/01/) if everything looks good.")
    print("3. Clean up SUMMARY.md to remove old daily entries.")

if __name__ == "__main__":
    main()
