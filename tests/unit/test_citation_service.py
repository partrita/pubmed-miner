"""
Unit tests for CitationService.
"""

import pytest
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

from src.pubmed_miner.services.citation_service import CitationService
from src.pubmed_miner.models.cache import CitationCache


class TestCitationService:
    """Test cases for CitationService class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.service = CitationService()

    def test_initialization(self):
        """Test service initialization."""
        assert self.service.cache_manager is not None
        assert hasattr(self.service, "crossref_base_url")
        assert hasattr(self.service, "semantic_scholar_base_url")

    @patch("src.pubmed_miner.services.citation_service.requests.get")
    def test_get_citation_count_crossref_success(self, mock_get):
        """Test successful citation count retrieval from Crossref."""
        # Mock Crossref API response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"message": {"is-referenced-by-count": 150}}
        mock_get.return_value = mock_response

        # Mock cache miss
        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = None

            with patch.object(
                self.service.cache_manager, "save_citation"
            ) as mock_cache_save:
                count = self.service.get_citation_count(
                    "12345", doi="10.1038/nature12345"
                )

                assert count == 150
                mock_cache_save.assert_called_once()

    @patch("src.pubmed_miner.services.citation_service.requests.get")
    def test_get_citation_count_semantic_scholar_fallback(self, mock_get):
        """Test fallback to Semantic Scholar when Crossref fails."""
        # Mock Crossref failure and Semantic Scholar success
        crossref_response = Mock()
        crossref_response.status_code = 404

        semantic_response = Mock()
        semantic_response.status_code = 200
        semantic_response.json.return_value = {"citationCount": 75}

        mock_get.side_effect = [crossref_response, semantic_response]

        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = None

            count = self.service.get_citation_count("12345", doi="10.1038/nature12345")
            assert count == 75

    def test_get_citation_count_from_cache(self):
        """Test citation count retrieval from cache."""
        # Mock cache hit
        cached_citation = CitationCache(
            pmid="12345",
            citation_count=100,
            last_updated=datetime.now(),
            source="crossref",
        )

        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = cached_citation

            count = self.service.get_citation_count("12345")
            assert count == 100

    def test_get_citation_count_expired_cache(self):
        """Test citation count with expired cache entry."""
        # Mock expired cache entry
        expired_citation = CitationCache(
            pmid="12345",
            citation_count=100,
            last_updated=datetime.now() - timedelta(days=10),  # Expired
            source="crossref",
        )

        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = expired_citation

            with patch.object(self.service, "_fetch_from_apis") as mock_fetch:
                mock_fetch.return_value = (150, "crossref")

                count = self.service.get_citation_count("12345")
                assert count == 150  # Should fetch new data

    def test_get_citation_count_invalid_pmid(self):
        """Test citation count with invalid PMID."""
        with pytest.raises(ValueError, match="PMID cannot be empty"):
            self.service.get_citation_count("")

        with pytest.raises(ValueError, match="PMID must be numeric"):
            self.service.get_citation_count("abc123")

    @patch("src.pubmed_miner.services.citation_service.requests.get")
    def test_get_citation_count_api_failure(self, mock_get):
        """Test citation count when all APIs fail."""
        # Mock all API failures
        mock_response = Mock()
        mock_response.status_code = 500
        mock_get.return_value = mock_response

        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = None

            count = self.service.get_citation_count("12345")
            assert count == 0  # Should return 0 when all APIs fail

    def test_batch_get_citation_counts(self):
        """Test batch citation count retrieval."""
        pmids = ["12345", "67890", "11111"]

        # Mock individual calls
        with patch.object(self.service, "get_citation_count") as mock_get:
            mock_get.side_effect = [100, 50, 200]

            counts = self.service.batch_get_citation_counts(pmids)

            assert len(counts) == 3
            assert counts["12345"] == 100
            assert counts["67890"] == 50
            assert counts["11111"] == 200

    def test_batch_get_citation_counts_with_errors(self):
        """Test batch citation count retrieval with some errors."""
        pmids = ["12345", "67890", "invalid"]

        def mock_get_side_effect(pmid, **kwargs):
            if pmid == "invalid":
                raise ValueError("Invalid PMID")
            return {"12345": 100, "67890": 50}[pmid]

        with patch.object(self.service, "get_citation_count") as mock_get:
            mock_get.side_effect = mock_get_side_effect

            counts = self.service.batch_get_citation_counts(pmids, skip_errors=True)

            assert len(counts) == 2
            assert counts["12345"] == 100
            assert counts["67890"] == 50
            assert "invalid" not in counts

    @patch("src.pubmed_miner.services.citation_service.requests.get")
    def test_crossref_api_rate_limiting(self, mock_get):
        """Test Crossref API rate limiting handling."""
        # Mock rate limit response
        mock_response = Mock()
        mock_response.status_code = 429
        mock_response.headers = {"Retry-After": "60"}
        mock_get.return_value = mock_response

        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = None

            with patch("time.sleep") as mock_sleep:
                count = self.service.get_citation_count("12345", doi="10.1038/test")

                # Should have slept due to rate limiting
                mock_sleep.assert_called_with(60)
                assert count == 0  # Should return 0 after rate limit

    def test_doi_validation(self):
        """Test DOI format validation."""
        # Valid DOIs
        assert self.service._validate_doi("10.1038/nature12345") is True
        assert self.service._validate_doi("10.1001/jama.2023.12345") is True

        # Invalid DOIs
        assert self.service._validate_doi("") is False
        assert self.service._validate_doi("invalid-doi") is False
        assert self.service._validate_doi("10.invalid") is False
        assert self.service._validate_doi(None) is False

    def test_pmid_to_doi_lookup(self):
        """Test PMID to DOI lookup functionality."""
        with patch(
            "src.pubmed_miner.services.citation_service.Entrez.efetch"
        ) as mock_efetch:
            mock_handle = Mock()
            mock_handle.read.return_value = """<?xml version="1.0" ?>
            <PubmedArticleSet>
                <PubmedArticle>
                    <PubmedData>
                        <ArticleIdList>
                            <ArticleId IdType="doi">10.1038/nature12345</ArticleId>
                        </ArticleIdList>
                    </PubmedData>
                </PubmedArticle>
            </PubmedArticleSet>"""
            mock_efetch.return_value = mock_handle

            doi = self.service._get_doi_from_pmid("12345")
            assert doi == "10.1038/nature12345"

    def test_pmid_to_doi_lookup_no_doi(self):
        """Test PMID to DOI lookup when no DOI is found."""
        with patch(
            "src.pubmed_miner.services.citation_service.Entrez.efetch"
        ) as mock_efetch:
            mock_handle = Mock()
            mock_handle.read.return_value = """<?xml version="1.0" ?>
            <PubmedArticleSet>
                <PubmedArticle>
                    <PubmedData>
                        <ArticleIdList>
                        </ArticleIdList>
                    </PubmedData>
                </PubmedArticle>
            </PubmedArticleSet>"""
            mock_efetch.return_value = mock_handle

            doi = self.service._get_doi_from_pmid("12345")
            assert doi is None

    def test_cache_expiry_check(self):
        """Test cache expiry checking."""
        # Recent cache entry (not expired)
        recent_cache = CitationCache(
            pmid="12345",
            citation_count=100,
            last_updated=datetime.now() - timedelta(hours=1),
            source="crossref",
        )
        assert not self.service._is_cache_expired(recent_cache)

        # Old cache entry (expired)
        old_cache = CitationCache(
            pmid="12345",
            citation_count=100,
            last_updated=datetime.now() - timedelta(days=10),
            source="crossref",
        )
        assert self.service._is_cache_expired(old_cache)

    def test_get_statistics(self):
        """Test statistics collection."""
        stats = self.service.get_statistics()

        required_keys = [
            "total_requests",
            "cache_hits",
            "cache_misses",
            "api_calls",
            "crossref_calls",
            "semantic_scholar_calls",
            "error_count",
            "average_response_time",
        ]

        for key in required_keys:
            assert key in stats
            assert isinstance(stats[key], (int, float))

    def test_clear_expired_cache(self):
        """Test clearing expired cache entries."""
        with patch.object(
            self.service.cache_manager, "clear_expired_citations"
        ) as mock_clear:
            mock_clear.return_value = 5

            cleared_count = self.service.clear_expired_cache()
            assert cleared_count == 5
            mock_clear.assert_called_once()

    @patch("src.pubmed_miner.services.citation_service.requests.get")
    def test_semantic_scholar_api_format(self, mock_get):
        """Test Semantic Scholar API response format handling."""
        # Mock Semantic Scholar response with different format
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "paperId": "test123",
            "citationCount": 42,
            "title": "Test Paper",
        }
        mock_get.return_value = mock_response

        with patch.object(self.service.cache_manager, "get_citation") as mock_cache_get:
            mock_cache_get.return_value = None

            count = self.service._fetch_from_semantic_scholar("12345")
            assert count == 42

    def test_concurrent_requests_thread_safety(self):
        """Test thread safety of concurrent requests."""
        import threading

        results = []
        errors = []

        def make_request(pmid):
            try:
                with patch.object(self.service, "_fetch_from_apis") as mock_fetch:
                    mock_fetch.return_value = (100, "test")
                    result = self.service.get_citation_count(pmid)
                    results.append(result)
            except Exception as e:
                errors.append(e)

        # Create multiple threads
        threads = []
        for i in range(10):
            thread = threading.Thread(target=make_request, args=(f"1234{i}",))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # All requests should succeed
        assert len(errors) == 0
        assert len(results) == 10

    def test_retry_mechanism(self):
        """Test retry mechanism for failed API calls."""
        with patch(
            "src.pubmed_miner.services.citation_service.requests.get"
        ) as mock_get:
            # First call fails, second succeeds
            failure_response = Mock()
            failure_response.status_code = 500

            success_response = Mock()
            success_response.status_code = 200
            success_response.json.return_value = {
                "message": {"is-referenced-by-count": 100}
            }

            mock_get.side_effect = [failure_response, success_response]

            with patch.object(
                self.service.cache_manager, "get_citation"
            ) as mock_cache_get:
                mock_cache_get.return_value = None

                count = self.service.get_citation_count("12345", doi="10.1038/test")
                assert count == 100
                assert mock_get.call_count == 2

    def test_malformed_api_response(self):
        """Test handling of malformed API responses."""
        with patch(
            "src.pubmed_miner.services.citation_service.requests.get"
        ) as mock_get:
            # Mock malformed JSON response
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.side_effect = ValueError("Invalid JSON")
            mock_get.return_value = mock_response

            with patch.object(
                self.service.cache_manager, "get_citation"
            ) as mock_cache_get:
                mock_cache_get.return_value = None

                count = self.service.get_citation_count("12345", doi="10.1038/test")
                assert count == 0  # Should handle gracefully

    def test_update_cache_settings(self):
        """Test updating cache settings."""
        new_settings = {"cache_expiry_days": 14, "max_cache_size": 10000}

        self.service.update_cache_settings(new_settings)

        # Verify settings were updated
        assert hasattr(self.service, "cache_expiry_days")
        assert self.service.cache_expiry_days == 14
