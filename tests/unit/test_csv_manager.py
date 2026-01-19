"""
Tests for CSV Manager functionality.
"""

import pytest
import tempfile
from pathlib import Path
from datetime import datetime

from src.pubmed_miner.models import Paper, ScoredPaper
from src.pubmed_miner.utils import CSVManager


@pytest.fixture
def temp_csv_file():
    """Create a temporary CSV file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        temp_path = f.name
    yield temp_path
    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


@pytest.fixture
def sample_papers():
    """Create sample Paper objects for testing."""
    return [
        Paper(
            pmid="12345678",
            title="Test Paper 1",
            authors=["Author One", "Author Two"],
            journal="Test Journal",
            publication_date=datetime(2024, 1, 1),
            doi="10.1234/test.2024.001",
            abstract="Test abstract 1",
            topic="machine-learning",
        ),
        Paper(
            pmid="87654321",
            title="Test Paper 2",
            authors=["Author Three"],
            journal="Another Journal",
            publication_date=datetime(2023, 12, 1),
            doi="10.1234/test.2023.001",
            abstract="Test abstract 2",
            topic="drug-discovery",
        ),
    ]


@pytest.fixture
def sample_scored_papers():
    """Create sample ScoredPaper objects for testing."""
    return [
        ScoredPaper(
            pmid="11111111",
            title="Scored Paper 1",
            authors=["Scorer One"],
            journal="Top Journal",
            publication_date=datetime(2024, 1, 1),
            citation_count=100,
            impact_factor=40.5,
            score=95.0,
            rank=1,
            doi="10.1234/scored.2024.001",
            topic="ai-drug-discovery",
        ),
        ScoredPaper(
            pmid="22222222",
            title="Scored Paper 2",
            authors=["Scorer Two"],
            journal="Good Journal",
            publication_date=datetime(2023, 11, 1),
            citation_count=50,
            impact_factor=20.3,
            score=75.0,
            rank=2,
            doi="10.1234/scored.2023.001",
            topic="cancer-immunotherapy",
        ),
    ]


class TestCSVManager:
    """Test suite for CSVManager."""

    def test_save_papers_creates_file(self, temp_csv_file, sample_papers):
        """Test that save_papers creates a CSV file."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()  # Delete temp file first
        
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        assert csv_path.exists()
        
    def test_save_papers_with_correct_headers(self, temp_csv_file, sample_papers):
        """Test that saved CSV has correct headers."""
        with open(temp_csv_file, 'r') as f:
            header_line = f.readline().strip()
        
        CSVManager.save_papers(sample_papers, temp_csv_file)
        
        with open(temp_csv_file, 'r') as f:
            actual_headers = f.readline().strip()
        
        expected_headers = ",".join(CSVManager.HEADERS)
        assert actual_headers == expected_headers
    
    def test_save_papers_contains_data(self, temp_csv_file, sample_papers):
        """Test that saved CSV contains paper data."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()  # Delete temp file first
        
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        with open(csv_path, 'r') as f:
            lines = f.readlines()
        
        # Header + 2 papers
        assert len(lines) == 3
        assert "Test Paper 1" in lines[1]
        assert "Test Paper 2" in lines[2]
    
    def test_save_scored_papers_with_scoring_info(
        self, temp_csv_file, sample_scored_papers
    ):
        """Test that scored papers include scoring information."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()  # Delete temp file first
        
        CSVManager.save_papers(
            sample_scored_papers, str(csv_path), include_scoring=True
        )
        
        with open(csv_path, 'r') as f:
            header_line = f.readline().strip()
        
        expected_headers = ",".join(CSVManager.SCORED_HEADERS)
        assert header_line == expected_headers
    
    def test_append_papers(self, temp_csv_file, sample_papers):
        """Test appending papers to existing CSV."""
        csv_path = Path(temp_csv_file)
        
        # Save initial papers
        CSVManager.save_papers(sample_papers[:1], str(csv_path))
        
        # Append more papers
        CSVManager.append_papers(sample_papers[1:], str(csv_path))
        
        # Load and verify
        with open(csv_path, 'r') as f:
            lines = f.readlines()
        
        # Header + 2 papers
        assert len(lines) == 3
    
    def test_load_papers(self, temp_csv_file, sample_papers):
        """Test loading papers from CSV."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        
        # Save papers
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        # Load papers
        loaded_papers = CSVManager.load_papers(str(csv_path))
        
        assert len(loaded_papers) == 2
        assert loaded_papers[0]["pmid"] == "12345678"
        assert loaded_papers[0]["title"] == "Test Paper 1"
        assert "Author One" in loaded_papers[0]["authors"]
    
    def test_save_papers_empty_list_raises_error(self, temp_csv_file):
        """Test that saving empty list raises ValueError."""
        with pytest.raises(ValueError):
            CSVManager.save_papers([], temp_csv_file)
    
    def test_append_papers_empty_list_raises_error(self, temp_csv_file, sample_papers):
        """Test that appending empty list raises ValueError."""
        # Create initial file first
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        # Try to append empty list
        with pytest.raises(ValueError):
            CSVManager.append_papers([], temp_csv_file)
    
    def test_load_nonexistent_file_raises_error(self):
        """Test that loading non-existent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            CSVManager.load_papers("/nonexistent/path/to/file.csv")
    
    def test_authors_joined_with_semicolon(self, temp_csv_file, sample_papers):
        """Test that multiple authors are joined with semicolon."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        with open(csv_path, 'r') as f:
            lines = f.readlines()
        
        # Check that authors are separated by semicolon
        assert "Author One; Author Two" in lines[1]
    
    def test_publication_date_as_isoformat(self, temp_csv_file, sample_papers):
        """Test that publication date is stored as ISO format."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        with open(csv_path, 'r') as f:
            lines = f.readlines()
        
        # Check ISO format date
        assert "2024-01-01" in lines[1]
    
    def test_update_collection(self, temp_csv_file, sample_scored_papers):
        """Test update_collection method."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        
        CSVManager.update_collection(sample_scored_papers, str(csv_path))
        
        with open(csv_path, 'r') as f:
            header = f.readline().strip()
        
        expected_headers = ",".join(CSVManager.SCORED_HEADERS)
        assert header == expected_headers
    
    def test_topic_field_included(self, temp_csv_file, sample_papers):
        """Test that topic field is included in CSV."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        with open(csv_path, 'r') as f:
            header = f.readline().strip()
        
        assert "topic" in header
        
    def test_topic_value_saved_correctly(self, temp_csv_file, sample_papers):
        """Test that topic values are saved correctly in CSV."""
        csv_path = Path(temp_csv_file)
        csv_path.unlink()
        
        CSVManager.save_papers(sample_papers, str(csv_path))
        
        loaded_papers = CSVManager.load_papers(str(csv_path))
        
        assert loaded_papers[0]["topic"] == "machine-learning"
        assert loaded_papers[1]["topic"] == "drug-discovery"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
