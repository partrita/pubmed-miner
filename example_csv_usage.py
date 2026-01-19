"""
Example script showing how to save collected papers to CSV.
"""

from datetime import datetime
from pathlib import Path

from src.pubmed_miner.models import Paper, ScoredPaper
from src.pubmed_miner.utils import CSVManager


def example_save_basic_papers():
    """Example: Save basic paper collection to CSV."""
    
    # Create sample papers with topic information
    papers = [
        Paper(
            pmid="12345678",
            title="Sample Paper 1: Machine Learning Applications",
            authors=["John Doe", "Jane Smith"],
            journal="Nature Machine Intelligence",
            publication_date=datetime(2024, 1, 15),
            doi="10.1234/sample.2024.001",
            abstract="This paper discusses modern ML applications...",
            topic="machine-learning-antibody",
        ),
        Paper(
            pmid="87654321",
            title="Sample Paper 2: Deep Learning Methods",
            authors=["Alice Johnson", "Bob Wilson"],
            journal="IEEE Transactions on Neural Networks",
            publication_date=datetime(2023, 12, 10),
            doi="10.1234/sample.2023.002",
            abstract="We present novel deep learning techniques...",
            topic="ai-drug-discovery",
        ),
    ]
    
    # Save to CSV
    output_path = Path("data/collections.csv")
    CSVManager.save_papers(papers, str(output_path))
    print(f"✓ Saved {len(papers)} papers to {output_path}")


def example_save_scored_papers():
    """Example: Save scored papers to CSV with ranking information."""
    
    # Create sample scored papers with topic information
    papers = [
        ScoredPaper(
            pmid="11111111",
            title="Highly Cited Research",
            authors=["Dr. Einstein"],
            journal="Science",
            publication_date=datetime(2023, 6, 1),
            citation_count=150,
            impact_factor=42.5,
            score=95.5,
            rank=1,
            doi="10.1234/score.2023.001",
            topic="cancer-immunotherapy",
        ),
        ScoredPaper(
            pmid="22222222",
            title="Important Study",
            authors=["Prof. Newton"],
            journal="Nature",
            publication_date=datetime(2023, 7, 15),
            citation_count=87,
            impact_factor=39.8,
            score=78.3,
            rank=2,
            doi="10.1234/score.2023.002",
            topic="de-novo-protein-design",
        ),
    ]
    
    # Save with scoring information
    output_path = Path("data/collections.csv")
    CSVManager.update_collection(papers, str(output_path))
    print(f"✓ Saved {len(papers)} scored papers to {output_path}")


def example_append_papers():
    """Example: Append new papers to existing CSV."""
    
    papers = [
        Paper(
            pmid="99999999",
            title="Additional Paper",
            authors=["New Author"],
            journal="Journal of Examples",
            publication_date=datetime(2024, 1, 19),
            doi="10.1234/additional.2024.001",
            topic="bioinformatics",
        ),
    ]
    
    output_path = Path("data/collections.csv")
    CSVManager.append_papers(papers, str(output_path))
    print(f"✓ Appended {len(papers)} papers to {output_path}")


def example_load_papers():
    """Example: Load papers from existing CSV."""
    
    input_path = Path("data/collections.csv")
    papers = CSVManager.load_papers(str(input_path))
    print(f"✓ Loaded {len(papers)} papers from {input_path}")
    
    # Display first paper
    if papers:
        print(f"\nFirst paper:")
        for key, value in papers[0].items():
            print(f"  {key}: {value}")


if __name__ == "__main__":
    print("📊 Paper Collection CSV Examples\n")
    
    # Run examples
    print("1. Saving basic papers...")
    example_save_basic_papers()
    
    print("\n2. Saving scored papers with ranking...")
    example_save_scored_papers()
    
    print("\n3. Appending additional papers...")
    example_append_papers()
    
    print("\n4. Loading papers from CSV...")
    example_load_papers()
    
    print("\n✓ All examples completed!")
