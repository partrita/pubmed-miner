"""
Unit tests for PaperCollectionService.
"""

import pytest
from unittest.mock import Mock, patch
from datetime import datetime

from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.utils.error_handler import APIError


class TestPaperCollectionService:
    """Test cases for PaperCollectionService class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.service = PaperCollectionService(
            email="test@example.com",
            rate_limit=10.0,  # Higher rate limit for testing
        )

    def test_initialization(self):
        """Test service initialization."""
        assert self.service.email == "test@example.com"
        assert self.service.rate_limit == 10.0

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_success(self, mock_esearch):
        """Test successful paper search."""
        # Mock Entrez response
        mock_handle = Mock()
        mock_handle.read.return_value = """<?xml version="1.0" ?>
        <eSearchResult>
            <IdList>
                <Id>12345</Id>
                <Id>67890</Id>
                <Id>11111</Id>
            </IdList>
        </eSearchResult>"""
        mock_esearch.return_value = mock_handle

        pmids = self.service.search_papers("machine learning", max_results=100)

        assert len(pmids) == 3
        assert "12345" in pmids
        assert "67890" in pmids
        assert "11111" in pmids

        # Verify Entrez.esearch was called correctly
        mock_esearch.assert_called_once()
        call_args = mock_esearch.call_args[1]
        assert call_args["db"] == "pubmed"
        assert call_args["term"] == "machine learning"
        assert call_args["retmax"] == 100

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_empty_query(self, mock_esearch):
        """Test search with empty query."""
        with pytest.raises(ValueError, match="Query cannot be empty"):
            self.service.search_papers("", max_results=100)

        with pytest.raises(ValueError, match="Query cannot be empty"):
            self.service.search_papers("   ", max_results=100)

    def test_search_papers_invalid_max_results(self):
        """Test search with invalid max_results."""
        with pytest.raises(ValueError, match="Max results must be positive"):
            self.service.search_papers("test", max_results=0)

        with pytest.raises(ValueError, match="Max results must be positive"):
            self.service.search_papers("test", max_results=-10)

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_no_results(self, mock_esearch):
        """Test search with no results."""
        mock_handle = Mock()
        mock_handle.read.return_value = """<?xml version="1.0" ?>
        <eSearchResult>
            <IdList>
            </IdList>
        </eSearchResult>"""
        mock_esearch.return_value = mock_handle

        pmids = self.service.search_papers("nonexistent query")

        assert pmids == []

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_api_error(self, mock_esearch):
        """Test search with API error."""
        mock_esearch.side_effect = Exception("PubMed API error")

        with pytest.raises(APIError, match="PubMed search failed"):
            self.service.search_papers("test query")

    @patch("src.pubmed_miner.services.paper_collection.Entrez.efetch")
    def test_get_paper_details_success(self, mock_efetch):
        """Test successful paper details retrieval."""
        mock_handle = Mock()
        mock_handle.read.return_value = """<?xml version="1.0" ?>
        <PubmedArticleSet>
            <PubmedArticle>
                <MedlineCitation>
                    <PMID>12345</PMID>
                    <Article>
                        <ArticleTitle>Test Paper Title</ArticleTitle>
                        <AuthorList>
                            <Author>
                                <LastName>Doe</LastName>
                                <ForeName>John</ForeName>
                            </Author>
                            <Author>
                                <LastName>Smith</LastName>
                                <ForeName>Jane</ForeName>
                            </Author>
                        </AuthorList>
                        <Journal>
                            <Title>Nature Medicine</Title>
                        </Journal>
                        <ArticleDate>
                            <Year>2023</Year>
                            <Month>01</Month>
                            <Day>15</Day>
                        </ArticleDate>
                        <Abstract>
                            <AbstractText>This is a test abstract.</AbstractText>
                        </Abstract>
                    </Article>
                </MedlineCitation>
                <PubmedData>
                    <ArticleIdList>
                        <ArticleId IdType="doi">10.1038/test.2023.12345</ArticleId>
                    </ArticleIdList>
                </PubmedData>
            </PubmedArticle>
        </PubmedArticleSet>"""
        mock_efetch.return_value = mock_handle

        papers = self.service.get_paper_details(["12345"])

        assert len(papers) == 1
        paper = papers[0]
        assert paper.pmid == "12345"
        assert paper.title == "Test Paper Title"
        assert len(paper.authors) == 2
        assert "John Doe" in paper.authors
        assert "Jane Smith" in paper.authors
        assert paper.journal == "Nature Medicine"
        assert paper.abstract == "This is a test abstract."
        assert paper.doi == "10.1038/test.2023.12345"

    def test_get_paper_details_empty_pmids(self):
        """Test paper details with empty PMID list."""
        papers = self.service.get_paper_details([])
        assert papers == []

    def test_get_paper_details_invalid_pmids(self):
        """Test paper details with invalid PMIDs."""
        with pytest.raises(ValueError, match="All PMIDs must be non-empty strings"):
            self.service.get_paper_details(["12345", "", "67890"])

        with pytest.raises(ValueError, match="All PMIDs must be non-empty strings"):
            self.service.get_paper_details([None, "12345"])

    @patch("src.pubmed_miner.services.paper_collection.Entrez.efetch")
    def test_get_paper_details_api_error(self, mock_efetch):
        """Test paper details with API error."""
        mock_efetch.side_effect = Exception("PubMed API error")

        with pytest.raises(APIError, match="Failed to fetch paper details"):
            self.service.get_paper_details(["12345"])

    @patch("src.pubmed_miner.services.paper_collection.time.sleep")
    def test_rate_limiting(self, mock_sleep):
        """Test rate limiting functionality."""
        # Set a very low rate limit
        service = PaperCollectionService(rate_limit=1.0)

        with patch.object(service, "_make_request") as mock_request:
            mock_request.return_value = []

            # Make multiple requests
            service.search_papers("test1")
            service.search_papers("test2")

            # Should have slept between requests
            assert mock_sleep.call_count >= 1

    def test_parse_authors_various_formats(self):
        """Test author parsing with various formats."""
        # Test with full names
        authors_xml = """
        <AuthorList>
            <Author>
                <LastName>Doe</LastName>
                <ForeName>John A</ForeName>
            </Author>
            <Author>
                <LastName>Smith</LastName>
                <ForeName>Jane</ForeName>
            </Author>
        </AuthorList>"""

        authors = self.service._parse_authors(authors_xml)
        assert "John A Doe" in authors
        assert "Jane Smith" in authors

    def test_parse_date_various_formats(self):
        """Test date parsing with various formats."""
        # Test complete date
        date_xml = """
        <ArticleDate>
            <Year>2023</Year>
            <Month>01</Month>
            <Day>15</Day>
        </ArticleDate>"""

        date = self.service._parse_date(date_xml)
        assert date == datetime(2023, 1, 15)

        # Test year only
        date_xml_year = """
        <PubDate>
            <Year>2023</Year>
        </PubDate>"""

        date = self.service._parse_date(date_xml_year)
        assert date.year == 2023
        assert date.month == 1
        assert date.day == 1

    def test_batch_processing(self):
        """Test batch processing of PMIDs."""
        pmids = [str(i) for i in range(250)]  # More than batch size

        with patch.object(self.service, "_fetch_batch") as mock_fetch:
            mock_fetch.return_value = []

            self.service.get_paper_details(pmids)

            # Should have made multiple batch calls
            assert mock_fetch.call_count > 1

    def test_retry_mechanism(self):
        """Test retry mechanism for failed requests."""
        with patch(
            "src.pubmed_miner.services.paper_collection.Entrez.esearch"
        ) as mock_esearch:
            # First call fails, second succeeds
            mock_esearch.side_effect = [
                Exception("Temporary error"),
                Mock(
                    read=Mock(
                        return_value="<eSearchResult><IdList></IdList></eSearchResult>"
                    )
                ),
            ]

            # Should retry and succeed
            pmids = self.service.search_papers("test query")
            assert pmids == []
            assert mock_esearch.call_count == 2

    def test_validate_pmid_format(self):
        """Test PMID format validation."""
        # Valid PMIDs
        assert self.service._validate_pmid("12345") is True
        assert self.service._validate_pmid("1") is True
        assert self.service._validate_pmid("12345678") is True

        # Invalid PMIDs
        assert self.service._validate_pmid("") is False
        assert self.service._validate_pmid("abc123") is False
        assert self.service._validate_pmid("123456789") is False  # Too long
        assert self.service._validate_pmid(None) is False

    def test_clean_text(self):
        """Test text cleaning functionality."""
        # Test HTML entity removal
        dirty_text = "This &amp; that &lt;tag&gt; &quot;quoted&quot;"
        clean_text = self.service._clean_text(dirty_text)
        assert clean_text == 'This & that <tag> "quoted"'

        # Test whitespace normalization
        messy_text = "  Multiple   spaces\n\nand\t\ttabs  "
        clean_text = self.service._clean_text(messy_text)
        assert clean_text == "Multiple spaces and tabs"

    def test_get_statistics(self):
        """Test statistics collection."""
        stats = self.service.get_statistics()

        required_keys = [
            "total_searches",
            "total_papers_fetched",
            "total_api_calls",
            "average_response_time",
            "error_count",
            "rate_limit_hits",
        ]

        for key in required_keys:
            assert key in stats
            assert isinstance(stats[key], (int, float))

    def test_reset_statistics(self):
        """Test statistics reset."""
        # Make some calls to generate statistics
        with patch.object(self.service, "_make_request") as mock_request:
            mock_request.return_value = []
            self.service.search_papers("test")

        # Reset and verify
        self.service.reset_statistics()
        stats = self.service.get_statistics()

        assert stats["total_searches"] == 0
        assert stats["total_api_calls"] == 0

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_with_filters(self, mock_esearch):
        """Test search with additional filters."""
        mock_handle = Mock()
        mock_handle.read.return_value = (
            "<eSearchResult><IdList></IdList></eSearchResult>"
        )
        mock_esearch.return_value = mock_handle

        # Test with date range
        self.service.search_papers(
            "machine learning", date_from="2020/01/01", date_to="2023/12/31"
        )

        call_args = mock_esearch.call_args[1]
        assert "2020/01/01" in call_args["term"]
        assert "2023/12/31" in call_args["term"]

    def test_concurrent_requests(self):
        """Test handling of concurrent requests."""
        import threading

        results = []
        errors = []

        def make_request():
            try:
                with patch.object(self.service, "_make_request") as mock_request:
                    mock_request.return_value = ["12345"]
                    result = self.service.search_papers("test")
                    results.append(result)
            except Exception as e:
                errors.append(e)

        # Create multiple threads
        threads = []
        for i in range(5):
            thread = threading.Thread(target=make_request)
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # All requests should succeed
        assert len(errors) == 0
        assert len(results) == 5
