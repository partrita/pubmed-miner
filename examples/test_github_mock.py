#!/usr/bin/env python3
"""
Local testing script that works without GitHub token.
This demonstrates how to use the system in mock mode for local development.
"""

import os
import sys
from datetime import datetime
from pathlib import Path

# Add src to path
# Path(__file__).parent is examples/
# Path(__file__).parent.parent is root/
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pubmed_miner.models import GitHubConfig, ScoredPaper
from pubmed_miner.services.github_manager import GitHubIssuesManager


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


def test_mock_mode():
    """Test GitHub manager in mock mode."""
    print("🧪 Testing PubMed Miner GitHub Integration (Mock Mode)")
    print("=" * 60)

    # Create configuration for mock mode
    config = GitHubConfig(
        token="mock_token_for_local_testing",  # This triggers mock mode
        repository="your-username/essential-papers-repo",
        issue_labels=["essential-papers", "automated", "research"],
    )

    print("📋 Configuration:")
    print(f"   Repository: {config.repository}")
    print(f"   Token: {config.token}")
    print(f"   Labels: {config.issue_labels}")
    print()

    # Initialize GitHub manager
    github_manager = GitHubIssuesManager(config)

    print("🔧 GitHub Manager Status:")
    print(f"   Mock Mode: {github_manager.mock_mode}")
    print(f"   Base URL: {github_manager.base_url}")
    print()

    # Create sample papers
    papers = create_sample_papers()
    print(f"📄 Sample Papers Created: {len(papers)}")
    for i, paper in enumerate(papers, 1):
        print(f"   {i}. {paper.title[:50]}... (Score: {paper.score:.1f})")
    print()

    # Test creating/updating issue
    print("🚀 Testing Issue Creation/Update...")
    try:
        result = github_manager.create_or_update_issue(
            topic="machine-learning-healthcare", papers=papers
        )

        print("✅ Success! Issue Details:")
        print(f"   Issue Number: #{result['number']}")
        print(f"   Title: {result['title']}")
        print(f"   URL: {result['html_url']}")
        print(f"   Created: {result['created_at']}")
        print()

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

    # Test finding existing issue
    print("🔍 Testing Issue Search...")
    try:
        existing = github_manager.find_existing_issue("machine-learning-healthcare")
        if existing:
            print(f"   Found existing issue: #{existing['number']}")
        else:
            print("   No existing issue found (expected in mock mode)")
        print()

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

    # Test paper list formatting
    print("📝 Testing Paper List Formatting...")
    try:
        formatted_body = github_manager.format_paper_list(papers)
        print(f"   Generated body length: {len(formatted_body)} characters")
        print(
            f"   Contains paper titles: {'✅' if papers[0].title in formatted_body else '❌'}"
        )
        print(f"   Contains rankings: {'✅' if 'Rank #1' in formatted_body else '❌'}")
        print(
            f"   Contains scores: {'✅' if str(papers[0].score) in formatted_body else '❌'}"
        )
        print()

        # Show a preview of the formatted content
        print("📋 Preview of Generated Issue Body:")
        print("-" * 40)
        preview_lines = formatted_body.split("\n")[:15]
        for line in preview_lines:
            print(f"   {line}")
        if len(formatted_body.split("\n")) > 15:
            remaining_lines = len(formatted_body.split("\n")) - 15
            print(f"   ... ({remaining_lines} more lines)")
        print("-" * 40)
        print()

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

    print("🎉 All tests passed! The system works in mock mode.")
    print()
    print("💡 To use with real GitHub:")
    print("   1. Set GITHUB_TOKEN environment variable")
    print("   2. Update repository in config/settings.yaml")
    print("   3. Run the system normally")

    return True


def test_environment_setup():
    """Test if environment is properly set up."""
    print("🔧 Environment Setup Check")
    print("=" * 30)

    # Check if GitHub token is available
    github_token = os.getenv("GITHUB_TOKEN")
    if github_token:
        print(f"✅ GITHUB_TOKEN found (length: {len(github_token)})")
        print("   → Real GitHub integration available")
    else:
        print("ℹ️  GITHUB_TOKEN not found")
        print("   → Will use mock mode for local testing")

    # Check config files
    config_dir = Path("config")
    if config_dir.exists():
        print(f"✅ Config directory found: {config_dir}")

        topics_file = config_dir / "topics.yaml"
        settings_file = config_dir / "settings.yaml"

        if topics_file.exists():
            print(f"✅ Topics config found: {topics_file}")
        else:
            print(f"⚠️  Topics config missing: {topics_file}")

        if settings_file.exists():
            print(f"✅ Settings config found: {settings_file}")
        else:
            print(f"⚠️  Settings config missing: {settings_file}")
    else:
        print(f"⚠️  Config directory missing: {config_dir}")

    print()


if __name__ == "__main__":
    print("🔬 PubMed Miner - Local Testing (No GitHub Token Required)")
    print("=" * 70)
    print()

    # Test environment
    test_environment_setup()

    # Test mock mode functionality
    success = test_mock_mode()

    if success:
        print("🎯 Summary: Local testing completed successfully!")
        print(
            "   The system is ready for development and testing without GitHub credentials."
        )
    else:
        print("❌ Summary: Some tests failed. Please check the error messages above.")
        sys.exit(1)
