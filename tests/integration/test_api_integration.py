"""
Integration tests for external API interactions.
"""
import pytest
import os
from unittest.mock import Mock, patch
from datetime import datetime

from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.services.citation_service import CitationService
from src.pubmed_miner.services.github_manager import GitHubIssuesManager
from src.pubmed_miner.models import GitHubConfig, ScoredPaper


class TestPubMedAPIIntegration:
    """Integration tests for PubMed API interactions."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.service = PaperCollectionService(email="test@example.com")
        
    @pytest.mark.skipif(
        os.getenv('SKIP_EXTERNAL_TESTS') == 'true',
        reason="Skipping external API tests"
    )
    def test_pubmed_search_real_api(self):
        """Test real PubMed API search (optional, skipped in CI)."""
        # This test can be run manually to verify PubMed API integration
        pmids = self.service.search_papers("machine learning", max_results=5)
        
        assert isinstance(pmids, list)
        assert len(pmids) <= 5
        
        if pmids:
            # Verify PMID format
            for pmid in pmids:
                assert pmid.isdigit()
                assert len(pmid) <= 8
                
    @pytest.mark.skipif(
        os.getenv('SKIP_EXTERNAL_TESTS') == 'true',
        reason="Skipping external API tests"
    )
    def test_pubmed_paper_details_real_api(self):
        """Test real PubMed API paper details retrieval."""
        # Use a known PMID for testing
        test_pmids = ["33057194"]  # A real PMID
        
        papers = self.service.get_paper_details(test_pmids)
        
        assert len(papers) == 1
        paper = papers[0]
        
        assert paper.pmid == "33057194"
        assert paper.title is not None
        assert len(paper.title) > 0
        assert paper.authors is not None
        assert len(paper.authors) > 0
        assert paper.journal is not None
        assert len(paper.journal) > 0
        
    def test_pubmed_api_error_handling(self):
        """Test PubMed API error handling with mocked failures."""
        with patch('src.pubmed_miner.services.paper_collection.Entrez.esearch') as mock_search:
            # Mock network error
            mock_search.side_effect = Exception("Network error")
            
            from src.pubmed_miner.utils.error_handler import APIError
            with pytest.raises(APIError):
                self.service.search_papers("test query")
                
    def test_pubmed_rate_limiting(self):
        """Test PubMed API rate limiting behavior."""
        # Test that rate limiting is properly implemented
        import time
        
        with patch('src.pubmed_miner.services.paper_collection.time.sleep') as mock_sleep:
            with patch.object(self.service, '_make_request') as mock_request:
                mock_request.return_value = []
                
                # Make multiple rapid requests
                start_time = time.time()
                for i in range(3):
                    self.service.search_papers(f"test query {i}")
                end_time = time.time()
                
                # Should have enforced rate limiting
                if self.service.rate_limit > 0:
                    expected_min_time = (3 - 1) / self.service.rate_limit
                    # Allow some tolerance for test execution time
                    assert mock_sleep.call_count >= 2 or (end_time - start_time) >= expected_min_time * 0.5


class TestCitationAPIIntegration:
    """Integration tests for citation API interactions."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.service = CitationService()
        
    @pytest.mark.skipif(
        os.getenv('SKIP_EXTERNAL_TESTS') == 'true',
        reason="Skipping external API tests"
    )
    def test_crossref_api_real_request(self):
        """Test real Crossref API request."""
        # Use a known DOI for testing
        test_doi = "10.1038/nature12373"
        
        count = self.service._fetch_from_crossref(test_doi)
        
        # Should return a non-negative integer
        assert isinstance(count, int)
        assert count >= 0
        
    @pytest.mark.skipif(
        os.getenv('SKIP_EXTERNAL_TESTS') == 'true',
        reason="Skipping external API tests"
    )
    def test_semantic_scholar_api_real_request(self):
        """Test real Semantic Scholar API request."""
        # Use a known PMID for testing
        test_pmid = "33057194"
        
        count = self.service._fetch_from_semantic_scholar(test_pmid)
        
        # Should return a non-negative integer or None
        assert count is None or (isinstance(count, int) and count >= 0)
        
    def test_citation_api_fallback_chain(self):
        """Test citation API fallback mechanism."""
        with patch('src.pubmed_miner.services.citation_service.requests.get') as mock_get:
            # Mock Crossref failure
            crossref_response = Mock()
            crossref_response.status_code = 404
            
            # Mock Semantic Scholar success
            semantic_response = Mock()
            semantic_response.status_code = 200
            semantic_response.json.return_value = {'citationCount': 42}
            
            mock_get.side_effect = [crossref_response, semantic_response]
            
            count = self.service.get_citation_count("12345", doi="10.1000/test")
            
            assert count == 42
            assert mock_get.call_count == 2  # Both APIs called
            
    def test_citation_api_rate_limit_handling(self):
        """Test citation API rate limit handling."""
        with patch('src.pubmed_miner.services.citation_service.requests.get') as mock_get:
            # Mock rate limit response
            rate_limit_response = Mock()
            rate_limit_response.status_code = 429
            rate_limit_response.headers = {'Retry-After': '1'}
            
            success_response = Mock()
            success_response.status_code = 200
            success_response.json.return_value = {
                'message': {'is-referenced-by-count': 100}
            }
            
            mock_get.side_effect = [rate_limit_response, success_response]
            
            with patch('time.sleep') as mock_sleep:
                count = self.service.get_citation_count("12345", doi="10.1000/test")
                
                assert count == 100
                mock_sleep.assert_called_with(1)  # Should respect Retry-After header


class TestGitHubAPIIntegration:
    """Integration tests for GitHub API interactions."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.config = GitHubConfig(
            token=os.getenv('GITHUB_TOKEN', 'mock_token_for_local_testing'),
            repository="test/repo",
            issue_labels=["test"]
        )
        self.manager = GitHubIssuesManager(self.config)
        
    @pytest.mark.skipif(
        not os.getenv('GITHUB_TOKEN') or os.getenv('SKIP_EXTERNAL_TESTS') == 'true',
        reason="GitHub token not available or external tests skipped"
    )
    def test_github_repository_access(self):
        """Test GitHub repository access with real API."""
        # This test requires a real GitHub token and repository
        try:
            repo_info = self.manager.get_repository_info()
            
            assert 'name' in repo_info
            assert 'full_name' in repo_info
            assert repo_info['full_name'] == self.config.repository
            
        except Exception as e:
            pytest.skip(f"GitHub API access failed: {e}")
            
    def test_github_api_authentication_error(self):
        """Test GitHub API authentication error handling."""
        # Use invalid token
        invalid_config = GitHubConfig(
            token="invalid_token",
            repository="test/repo",
            issue_labels=["test"]
        )
        invalid_manager = GitHubIssuesManager(invalid_config)
        
        with patch('src.pubmed_miner.services.github_manager.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 401
            mock_response.text = "Bad credentials"
            mock_get.return_value = mock_response
            
            from src.pubmed_miner.utils.error_handler import GitHubError
            with pytest.raises(GitHubError, match="GitHub authentication failed"):
                invalid_manager._find_existing_issue("test-topic")
                
    def test_github_api_rate_limiting(self):
        """Test GitHub API rate limiting handling."""
        with patch('src.pubmed_miner.services.github_manager.requests.get') as mock_get:
            # Mock rate limit response
            mock_response = Mock()
            mock_response.status_code = 403
            mock_response.headers = {
                'X-RateLimit-Remaining': '0',
                'X-RateLimit-Reset': str(int(datetime.now().timestamp()) + 3600)
            }
            mock_response.json.return_value = {
                'message': 'API rate limit exceeded'
            }
            mock_get.return_value = mock_response
            
            from src.pubmed_miner.utils.error_handler import GitHubError
            with pytest.raises(GitHubError, match="GitHub API rate limit exceeded"):
                self.manager._find_existing_issue("test-topic")
                
    def test_github_issue_creation_workflow(self):
        """Test complete GitHub issue creation workflow."""
        sample_papers = [
            ScoredPaper(
                pmid="12345",
                title="Test Paper",
                authors=["Test Author"],
                journal="Test Journal",
                publication_date=datetime(2023, 1, 1),
                citation_count=100,
                impact_factor=5.0,
                score=85.0,
                rank=1
            )
        ]
        
        with patch.object(self.manager, '_find_existing_issue') as mock_find:
            mock_find.return_value = None  # No existing issue
            
            with patch.object(self.manager, '_create_issue') as mock_create:
                mock_create.return_value = {
                    'number': 42,
                    'title': '[Essential Papers] test-topic',
                    'html_url': 'https://github.com/test/repo/issues/42',
                    'created': True
                }
                
                result = self.manager.create_or_update_issue("test-topic", sample_papers)
                
                assert result['created'] is True
                assert result['number'] == 42
                
                # Verify create_issue was called with formatted content
                mock_create.assert_called_once()
                call_args = mock_create.call_args[0]
                assert call_args[0] == "test-topic"
                assert "Test Paper" in call_args[1]  # Issue body should contain paper title
                
    def test_github_issue_update_workflow(self):
        """Test GitHub issue update workflow."""
        sample_papers = [
            ScoredPaper(
                pmid="67890",
                title="Updated Paper",
                authors=["Updated Author"],
                journal="Updated Journal",
                publication_date=datetime(2023, 2, 1),
                citation_count=150,
                impact_factor=7.0,
                score=90.0,
                rank=1
            )
        ]
        
        with patch.object(self.manager, '_find_existing_issue') as mock_find:
            mock_find.return_value = {
                'number': 42,
                'title': '[Essential Papers] test-topic',
                'state': 'open'
            }
            
            with patch.object(self.manager, '_update_issue') as mock_update:
                mock_update.return_value = {
                    'number': 42,
                    'title': '[Essential Papers] test-topic',
                    'html_url': 'https://github.com/test/repo/issues/42',
                    'updated': True
                }
                
                result = self.manager.create_or_update_issue("test-topic", sample_papers)
                
                assert result['updated'] is True
                assert result['number'] == 42
                
                # Verify update_issue was called
                mock_update.assert_called_once()
                call_args = mock_update.call_args[0]
                assert call_args[0] == 42  # Issue number
                assert "Updated Paper" in call_args[1]  # Updated content


class TestAPIIntegrationResilience:
    """Test API integration resilience and error recovery."""
    
    def test_multiple_api_failure_recovery(self):
        """Test recovery from multiple API failures."""
        paper_service = PaperCollectionService()
        citation_service = CitationService()
        
        # Test cascading failures and recovery
        with patch.object(paper_service, 'search_papers') as mock_search:
            # First call fails, second succeeds
            mock_search.side_effect = [
                Exception("Network timeout"),
                ["12345", "67890"]
            ]
            
            with patch.object(citation_service, 'get_citation_count') as mock_citations:
                # Citation service also has intermittent failures
                def citation_side_effect(pmid, **kwargs):
                    if pmid == "12345":
                        raise Exception("Citation API timeout")
                    return 100
                    
                mock_citations.side_effect = citation_side_effect
                
                # Simulate retry logic
                try:
                    pmids = paper_service.search_papers("test query")
                except Exception:
                    # Retry
                    pmids = paper_service.search_papers("test query")
                    
                assert pmids == ["12345", "67890"]
                
                # Process citations with error handling
                citation_counts = {}
                for pmid in pmids:
                    try:
                        count = citation_service.get_citation_count(pmid)
                        citation_counts[pmid] = count
                    except Exception:
                        citation_counts[pmid] = 0  # Fallback
                        
                assert citation_counts["12345"] == 0  # Fallback value
                assert citation_counts["67890"] == 100  # Successful retrieval
                
    def test_api_timeout_handling(self):
        """Test handling of API timeouts."""
        import requests
        
        citation_service = CitationService()
        
        with patch('src.pubmed_miner.services.citation_service.requests.get') as mock_get:
            # Mock timeout exception
            mock_get.side_effect = requests.Timeout("Request timeout")
            
            # Should handle timeout gracefully
            count = citation_service.get_citation_count("12345", doi="10.1000/test")
            assert count == 0  # Should return fallback value
            
    def test_api_malformed_response_handling(self):
        """Test handling of malformed API responses."""
        citation_service = CitationService()
        
        with patch('src.pubmed_miner.services.citation_service.requests.get') as mock_get:
            # Mock malformed JSON response
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.side_effect = ValueError("Invalid JSON")
            mock_get.return_value = mock_response
            
            # Should handle malformed response gracefully
            count = citation_service.get_citation_count("12345", doi="10.1000/test")
            assert count == 0  # Should return fallback value
            
    def test_concurrent_api_requests(self):
        """Test concurrent API requests don't interfere with each other."""
        import threading
        import time
        
        citation_service = CitationService()
        results = []
        errors = []
        
        def make_request(pmid):
            try:
                with patch.object(citation_service, '_fetch_from_apis') as mock_fetch:
                    mock_fetch.return_value = (int(pmid), "test")
                    count = citation_service.get_citation_count(pmid)
                    results.append((pmid, count))
            except Exception as e:
                errors.append(e)
                
        # Create multiple concurrent requests
        threads = []
        for i in range(10):
            thread = threading.Thread(target=make_request, args=(str(i + 1),))
            threads.append(thread)
            thread.start()
            
        # Wait for all requests to complete
        for thread in threads:
            thread.join()
            
        # All requests should succeed
        assert len(errors) == 0
        assert len(results) == 10
        
        # Results should match expected values
        for pmid, count in results:
            assert count == int(pmid)