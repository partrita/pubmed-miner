"""
Integration example: Collect papers and save to CSV.
This demonstrates the complete workflow from PubMed search to CSV storage.
"""


from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.utils import CSVManager


def collect_and_save_papers(
    query: str,
    topic: str | None = None,
    output_csv: str = "data/collections.csv",
    max_results: int = 100,
    email: str = "pubmed.miner@example.com",
):
    """
    Collect papers from PubMed based on a search query and save to CSV.

    Args:
        query: PubMed search query (e.g., "COVID-19 AND vaccine")
        topic: Topic name to associate with papers
        output_csv: Path to save the CSV file
        max_results: Maximum number of papers to collect (default: 100)
        email: Email for Entrez API (required by NCBI)

    Returns:
        Number of papers saved

    Example:
        >>> num_saved = collect_and_save_papers(
        ...     query="machine learning AND drug discovery",
        ...     topic="ai-drug-discovery",
        ...     max_results=50
        ... )
        >>> print(f"Saved {num_saved} papers to data/collections.csv")
    """
    print(f"🔍 Starting paper collection: '{query}'")
    if topic:
        print(f"📂 Topic: {topic}")
    print(f"📊 Maximum results: {max_results}")

    # Initialize service
    service = PaperCollectionService(email=email)

    # Step 1: Search PubMed
    print("\n1️⃣  Searching PubMed...")
    try:
        pmids = service.search_papers(query, max_results=max_results)
        print(f"   ✓ Found {len(pmids)} papers")
    except Exception as e:
        print(f"   ✗ Search failed: {e}")
        return 0

    if not pmids:
        print("   ✗ No papers found for this query")
        return 0

    # Step 2: Get paper details
    print("\n2️⃣  Retrieving paper details...")
    try:
        papers = service.get_paper_details(pmids, topic=topic)
        print(f"   ✓ Retrieved details for {len(papers)} papers")
    except Exception as e:
        print(f"   ✗ Failed to retrieve details: {e}")
        return 0

    if not papers:
        print("   ✗ Could not retrieve any paper details")
        return 0

    # Step 3: Save to CSV
    print("\n3️⃣  Saving to CSV...")
    try:
        CSVManager.save_papers(papers, output_csv)
        print(f"   ✓ Saved {len(papers)} papers to {output_csv}")
        return len(papers)
    except Exception as e:
        print(f"   ✗ Failed to save to CSV: {e}")
        return 0


def append_papers_to_collection(
    query: str,
    topic: str | None = None,
    output_csv: str = "data/collections.csv",
    max_results: int = 50,
    email: str = "pubmed.miner@example.com",
):
    """
    Append newly collected papers to existing CSV collection.

    Args:
        query: PubMed search query
        topic: Topic name to associate with papers
        output_csv: Path to CSV file to append to
        max_results: Maximum number of papers to collect
        email: Email for Entrez API

    Returns:
        Number of papers appended

    Example:
        >>> num_appended = append_papers_to_collection(
        ...     query="CRISPR gene therapy",
        ...     topic="cancer-immunotherapy"
        ... )
        >>> print(f"Appended {num_appended} papers")
    """
    print(f"🔍 Appending papers: '{query}'")
    if topic:
        print(f"📂 Topic: {topic}")

    # Initialize service
    service = PaperCollectionService(email=email)

    # Search and collect papers
    print("\n1️⃣  Searching PubMed...")
    try:
        pmids = service.search_papers(query, max_results=max_results)
        print(f"   ✓ Found {len(pmids)} papers")
    except Exception as e:
        print(f"   ✗ Search failed: {e}")
        return 0

    if not pmids:
        print("   ✗ No papers found")
        return 0

    # Get details
    print("\n2️⃣  Retrieving paper details...")
    try:
        papers = service.get_paper_details(pmids, topic=topic)
        print(f"   ✓ Retrieved {len(papers)} papers")
    except Exception as e:
        print(f"   ✗ Failed: {e}")
        return 0

    # Append to CSV
    print("\n3️⃣  Appending to CSV...")
    try:
        CSVManager.append_papers(papers, output_csv)
        print(f"   ✓ Appended {len(papers)} papers to {output_csv}")
        return len(papers)
    except Exception as e:
        print(f"   ✗ Failed: {e}")
        return 0


if __name__ == "__main__":
    print("=" * 60)
    print("📚 PubMed Paper Collection to CSV with Topics")
    print("=" * 60)

    # Example 1: Create new collection with topic
    print("\n📝 Example 1: Create new collection with topic")
    print("-" * 60)
    num_saved = collect_and_save_papers(
        query="artificial intelligence AND healthcare",
        topic="ai-healthcare",
        output_csv="data/collections.csv",
        max_results=10,  # Small number for testing
    )

    if num_saved > 0:
        print(f"\n✓ Successfully saved {num_saved} papers!")

        # Example 2: Append more papers with different topic
        print("\n\n📝 Example 2: Append additional papers with different topic")
        print("-" * 60)
        num_appended = append_papers_to_collection(
            query="machine learning AND drug discovery",
            topic="ai-drug-discovery",
            output_csv="data/collections.csv",
            max_results=10,
        )

        if num_appended > 0:
            print(f"\n✓ Successfully appended {num_appended} papers!")
        else:
            print("\n✗ Append operation failed")
    else:
        print("\n✗ Initial collection failed")

    print("\n" + "=" * 60)
    print("✓ Demo completed!")
    print("=" * 60)
