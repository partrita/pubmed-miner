import pytest
from datetime import datetime
from src.pubmed_miner.services.mdbook_manager import MdBookManager
from src.pubmed_miner.models.paper import ScoredPaper

def test_xss_mitigation():
    manager = MdBookManager(book_root="book_test")
    paper = ScoredPaper(
        pmid="123",
        title="<script>alert('title')</script>",
        authors=["<script>alert('author')</script>"],
        journal="<script>alert('journal')</script>",
        publication_date=datetime.now(),
        abstract="<script>alert('abstract')</script>",
        doi="10.123/456",
        score=1.0,
        rank=1
    )
    content = manager._format_page_content("Topic", [paper], datetime.now())
    assert "<script>" not in content
    assert "&lt;script&gt;alert(&#x27;title&#x27;)&lt;/script&gt;" in content
    assert "&lt;script&gt;alert(&#x27;author&#x27;)&lt;/script&gt;" in content
    assert "&lt;script&gt;alert(&#x27;journal&#x27;)&lt;/script&gt;" in content
    assert "&lt;script&gt;alert(&#x27;abstract&#x27;)&lt;/script&gt;" in content

def test_update_monthly_page_sanitizes_xss(tmp_path):
    """Verify that update_monthly_page sanitizes external input to prevent XSS."""
    manager = MdBookManager(book_root=str(tmp_path))

    # Create a paper with malicious HTML in fields
    paper = ScoredPaper(
        pmid="12345\"><script>alert('pmid')</script>",
        title="<script>alert('title')</script> & Normal Title | With Pipe",
        authors=["Normal Author"],
        journal="<script>alert('journal')</script> Malicious Journal",
        publication_date=datetime(2024, 1, 1),
        citation_count=100,
        impact_factor=5.0,
        score=90.0,
        rank=1,
        doi="10.1234/\"><script>alert('doi')</script>",
        topic="normal-topic"
    )

    malicious_topic = "<script>alert('topic')</script>"
    date = datetime(2024, 5, 15)

    # Call the method
    manager.update_monthly_page(topic=malicious_topic, papers=[paper], date=date)

    # Read the generated file
    file_path = tmp_path / "book_src" / "2024" / "05.md"
    assert file_path.exists()

    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Verify that malicious inputs were escaped
    assert "<script>" not in content
    assert "&lt;script&gt;alert(&#x27;pmid&#x27;)&lt;/script&gt;" in content
    assert "&lt;script&gt;alert(&#x27;title&#x27;)&lt;/script&gt; &amp; Normal Title \\| With Pipe" in content
    assert "&lt;script&gt;alert(&#x27;journal&#x27;)&lt;/script&gt; Malicious Journal" in content
    assert "&lt;script&gt;alert(&#x27;doi&#x27;)&lt;/script&gt;" in content
    assert "&lt;script&gt;alert(&#x27;topic&#x27;)&lt;/script&gt;" in content
