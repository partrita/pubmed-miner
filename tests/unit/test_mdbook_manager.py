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
