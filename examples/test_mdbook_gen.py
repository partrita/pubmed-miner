#!/usr/bin/env python3
"""
Local testing script for MdBook integration.
This demonstrates how to use the system locally to generate mdbook content.
"""

import os
import sys
from datetime import datetime
from pathlib import Path

# Add src to path
# Path(__file__).parent is examples/
# Path(__file__).parent.parent is root/
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pubmed_miner.models import ScoredPaper
from pubmed_miner.services.mdbook_manager import MdBookManager


def create_sample_papers():
    """Create sample papers for testing."""
    return [
        ScoredPaper(
            pmid="12345678",
            title="CRISPR-Cas9 gene editing: A revolutionary tool for precision medicine",
            authors=["Jennifer Doudna", "Emmanuelle Charpentier", "Feng Zhang"],
            journal="Nature",
            publication_date=datetime(2023, 3, 15),
            abstract="This study demonstrates the application of CRISPR-Cas9 technology for precise gene editing in human cells, showing promising results for treating genetic disorders.",
            doi="10.1038/nature12345",
            citation_count=1250,
            impact_factor=64.8,
            score=95.8,
            rank=1,
        ),
        ScoredPaper(
            pmid="87654321",
            title="Machine learning approaches for drug discovery and development",
            authors=["Alex Smith", "Maria Garcia", "John Chen", "Sarah Kim"],
            journal="Nature Biotechnology",
            publication_date=datetime(2022, 11, 20),
            abstract="We present novel machine learning algorithms that significantly accelerate the drug discovery process by predicting molecular interactions and optimizing compound design.",
            doi="10.1038/nbt.4567",
            citation_count=890,
            impact_factor=68.2,
            score=92.3,
            rank=2,
        ),
        ScoredPaper(
            pmid="11223344",
            title="COVID-19 vaccine efficacy and safety in real-world populations",
            authors=["Sarah Johnson", "Michael Brown", "Lisa Wang"],
            journal="New England Journal of Medicine",
            publication_date=datetime(2021, 8, 10),
            abstract="Comprehensive analysis of COVID-19 vaccine effectiveness across diverse populations, demonstrating high efficacy and acceptable safety profiles in real-world settings.",
            doi="10.1056/NEJMoa2021234",
            citation_count=2100,
            impact_factor=176.1,
            score=98.1,
            rank=3,
        ),
    ]


def test_mdbook_generation():
    """Test MdBook manager generation."""
    print("🧪 Testing PubMed Miner MdBook Integration")
    print("=" * 60)

    # Initialize MdBook manager
    mdbook_manager = MdBookManager()

    print("🔧 MdBook Manager Initialized")
    print(f"   Book Root: {mdbook_manager.book_root}")
    print(f"   Source Dir: {mdbook_manager.src_dir}")
    print()

    # Create sample papers
    papers = create_sample_papers()
    print(f"📄 Sample Papers Created: {len(papers)}")
    
    # Test creating monthly page
    print("🚀 Testing Monthly Page Creation...")
    try:
        topic = "Biomedical Research"
        relative_path = mdbook_manager.update_monthly_page(topic, papers)

        print("✅ Success! Page Updated/Created:")
        print(f"   Path: {relative_path}")
        
        full_path = mdbook_manager.src_dir / relative_path
        if full_path.exists():
            print(f"   File exists at: {full_path}")
            print(f"   Size: {full_path.stat().st_size} bytes")
        else:
            print(f"❌ File not found at expected path: {full_path}")
            return False

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

    # Test updating summary
    print("\n📚 Testing Summary Update...")
    try:
        mdbook_manager.update_summary(relative_path, topic)
        
        if mdbook_manager.summary_path.exists():
            print(f"✅ SUMMARY.md updated")
            with open(mdbook_manager.summary_path, 'r') as f:
                content = f.read()
                print("   Preview of SUMMARY.md:")
                print("-" * 40)
                print("\n".join(content.splitlines()[-5:])) # Show last 5 lines
                print("-" * 40)
                
            if relative_path in content:
                print("   Link successfully added to summary")
            else:
                print("❌ Link not found in summary")
                return False
        else:
            print("❌ SUMMARY.md not found")
            return False

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

    print("\n🎉 All tests passed!")
    print("\n💡 To view the book locally:")
    print("   1. Install mdbook (if not installed)")
    print("   2. Run 'mdbook serve'")
    print("   3. Open http://localhost:3000")

    return True


if __name__ == "__main__":
    test_mdbook_generation()
