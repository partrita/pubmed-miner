import pytest
from datetime import datetime
from src.pubmed_miner.models import Paper
from src.pubmed_miner.utils.csv_manager import CSVManager

def test_paper_to_dict_sanitizes_csv_injection():
    """Verify that _paper_to_dict prepends a single quote to strings starting with special chars."""
    paper = Paper(
        pmid="12345",
        title="=cmd|' /C calc'!A0",  # CSV injection payload
        authors=["+Author One"],      # CSV injection payload
        journal="-Malicious Journal", # CSV injection payload
        publication_date=datetime(2024, 1, 1),
        doi="@doi123",                # CSV injection payload
        abstract="Normal abstract",
        topic="test-topic"
    )

    result = CSVManager._paper_to_dict(paper)

    assert result["title"] == "'=cmd|' /C calc'!A0"
    assert result["authors"] == "'+Author One"
    assert result["journal"] == "'-Malicious Journal"
    assert result["doi"] == "'@doi123"
    assert result["abstract"] == "Normal abstract"
    assert result["topic"] == "test-topic"
