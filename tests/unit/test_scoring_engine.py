"""
Unit tests for the scoring engine.
"""
import pytest
from datetime import datetime, timedelta
from unittest.mock import Mock

from src.pubmed_miner.models import Paper, ScoredPaper, ScoringWeights
from src.pubmed_miner.scoring import ScoringEngine


class TestScoringEngine:
    """Test cases for ScoringEngine class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.engine = ScoringEngine()
        
        # Sample papers for testing
        self.sample_paper = Paper(
            pmid="12345",
            title="Machine Learning in Healthcare: A Comprehensive Review",
            authors=["John Doe", "Jane Smith"],
            journal="Nature Medicine",
            publication_date=datetime(2023, 1, 15),
            abstract="This paper reviews the applications of machine learning in healthcare...",
            doi="10.1038/s41591-023-12345"
        )
        
        self.old_paper = Paper(
            pmid="67890",
            title="Traditional Methods in Medical Research",
            authors=["Bob Johnson"],
            journal="Journal of Medicine",
            publication_date=datetime(2015, 6, 1),
            abstract="This paper discusses traditional research methods...",
            doi="10.1001/jama.2015.67890"
        )
        
    def test_initialization_default_weights(self):
        """Test engine initialization with default weights."""
        engine = ScoringEngine()
        assert engine.weights.citation_weight == 0.4
        assert engine.weights.impact_factor_weight == 0.3
        assert engine.weights.recency_weight == 0.2
        assert engine.weights.relevance_weight == 0.1
        
    def test_initialization_custom_weights(self):
        """Test engine initialization with custom weights."""
        custom_weights = ScoringWeights(
            citation_weight=0.5,
            impact_factor_weight=0.3,
            recency_weight=0.1,
            relevance_weight=0.1
        )
        engine = ScoringEngine(custom_weights)
        assert engine.weights.citation_weight == 0.5
        
    def test_calculate_citation_score(self):
        """Test citation score calculation."""
        # Test zero citations
        assert self.engine._calculate_citation_score(0) == 0.0
        
        # Test low citations
        score_10 = self.engine._calculate_citation_score(10)
        assert 0 < score_10 < 50
        
        # Test medium citations
        score_100 = self.engine._calculate_citation_score(100)
        assert 50 < score_100 < 90
        
        # Test high citations
        score_1000 = self.engine._calculate_citation_score(1000)
        assert score_1000 == 100.0
        
        # Test very high citations
        score_2000 = self.engine._calculate_citation_score(2000)
        assert score_2000 == 100.0
        
        # Verify monotonic increase
        assert score_10 < score_100 < score_1000
        
    def test_calculate_impact_factor_score(self):
        """Test impact factor score calculation."""
        # Test zero impact factor
        assert self.engine._calculate_impact_factor_score(0) == 0.0
        
        # Test low impact factor
        score_2 = self.engine._calculate_impact_factor_score(2.0)
        assert 0 < score_2 < 50
        
        # Test medium impact factor
        score_10 = self.engine._calculate_impact_factor_score(10.0)
        assert 50 < score_10 < 90
        
        # Test high impact factor
        score_50 = self.engine._calculate_impact_factor_score(50.0)
        assert score_50 == 100.0
        
        # Test very high impact factor
        score_100 = self.engine._calculate_impact_factor_score(100.0)
        assert score_100 == 100.0
        
        # Verify monotonic increase
        assert score_2 < score_10 < score_50
        
    def test_calculate_recency_score(self):
        """Test recency score calculation."""
        current_date = datetime.now()
        
        # Test very recent paper (1 month old)
        recent_date = current_date - timedelta(days=30)
        recent_score = self.engine._calculate_recency_score(recent_date)
        assert 90 < recent_score <= 100
        
        # Test moderately old paper (2 years old)
        medium_date = current_date - timedelta(days=730)
        medium_score = self.engine._calculate_recency_score(medium_date)
        assert 40 < medium_score < 80
        
        # Test old paper (6 years old)
        old_date = current_date - timedelta(days=2190)
        old_score = self.engine._calculate_recency_score(old_date)
        assert old_score == 10.0
        
        # Test future paper (shouldn't happen but handle gracefully)
        future_date = current_date + timedelta(days=30)
        future_score = self.engine._calculate_recency_score(future_date)
        assert future_score == 100.0
        
        # Verify monotonic decrease with age
        assert recent_score > medium_score > old_score
        
    def test_calculate_relevance_score(self):
        """Test relevance score calculation."""
        query = "machine learning healthcare"
        
        # Test with matching paper
        relevance_score = self.engine._calculate_relevance_score(self.sample_paper, query)
        assert relevance_score > 50  # Should have good relevance
        
        # Test with non-matching paper
        non_matching_score = self.engine._calculate_relevance_score(self.old_paper, query)
        assert non_matching_score < relevance_score  # Should be lower
        
        # Test with no query
        no_query_score = self.engine._calculate_relevance_score(self.sample_paper, None)
        assert no_query_score == 50.0  # Neutral score
        
        # Test with empty query
        empty_query_score = self.engine._calculate_relevance_score(self.sample_paper, "")
        assert empty_query_score == 50.0
        
    def test_calculate_paper_score(self):
        """Test overall paper score calculation."""
        citations = 150
        impact_factor = 25.0
        query = "machine learning healthcare"
        
        score = self.engine.calculate_paper_score(
            self.sample_paper, citations, impact_factor, query
        )
        
        # Score should be between 0 and 100
        assert 0 <= score <= 100
        
        # Recent paper with good citations and IF should score well
        assert score > 70
        
    def test_calculate_paper_score_components(self):
        """Test that score components are properly weighted."""
        citations = 100
        impact_factor = 10.0
        
        # Calculate individual components
        citation_score = self.engine._calculate_citation_score(citations)
        impact_score = self.engine._calculate_impact_factor_score(impact_factor)
        recency_score = self.engine._calculate_recency_score(self.sample_paper.publication_date)
        
        # Calculate expected total
        expected_total = (
            citation_score * self.engine.weights.citation_weight +
            impact_score * self.engine.weights.impact_factor_weight +
            recency_score * self.engine.weights.recency_weight +
            50.0 * self.engine.weights.relevance_weight  # No query = 50.0
        )
        
        actual_score = self.engine.calculate_paper_score(
            self.sample_paper, citations, impact_factor
        )
        
        assert abs(actual_score - expected_total) < 0.1
        
    def test_score_paper_batch(self):
        """Test batch scoring of papers."""
        papers = [self.sample_paper, self.old_paper]
        citations = {"12345": 150, "67890": 50}
        impact_factors = {"Nature Medicine": 87.2, "Journal of Medicine": 5.0}
        query = "machine learning"
        
        scored_papers = self.engine.score_paper_batch(
            papers, citations, impact_factors, query
        )
        
        assert len(scored_papers) == 2
        assert all(isinstance(p, ScoredPaper) for p in scored_papers)
        assert all(0 <= p.score <= 100 for p in scored_papers)
        
        # Recent paper should score higher
        recent_paper = next(p for p in scored_papers if p.pmid == "12345")
        old_paper = next(p for p in scored_papers if p.pmid == "67890")
        assert recent_paper.score > old_paper.score
        
    def test_rank_papers(self):
        """Test paper ranking functionality."""
        # Create sample scored papers
        scored_papers = [
            ScoredPaper(
                pmid="1", title="Paper 1", authors=["A"], journal="J1",
                publication_date=datetime.now(), citation_count=100,
                impact_factor=10.0, score=85.0, rank=0
            ),
            ScoredPaper(
                pmid="2", title="Paper 2", authors=["B"], journal="J2",
                publication_date=datetime.now(), citation_count=50,
                impact_factor=5.0, score=65.0, rank=0
            ),
            ScoredPaper(
                pmid="3", title="Paper 3", authors=["C"], journal="J3",
                publication_date=datetime.now(), citation_count=200,
                impact_factor=20.0, score=95.0, rank=0
            )
        ]
        
        ranked_papers = self.engine.rank_papers(scored_papers)
        
        # Check ranking order (highest score first)
        assert ranked_papers[0].pmid == "3"  # Score 95.0
        assert ranked_papers[1].pmid == "1"  # Score 85.0
        assert ranked_papers[2].pmid == "2"  # Score 65.0
        
        # Check rank assignment
        assert ranked_papers[0].rank == 1
        assert ranked_papers[1].rank == 2
        assert ranked_papers[2].rank == 3
        
    def test_select_essential_papers(self):
        """Test essential paper selection."""
        # Create sample scored papers
        scored_papers = []
        for i in range(20):
            paper = ScoredPaper(
                pmid=str(i), title=f"Paper {i}", authors=[f"Author {i}"],
                journal=f"Journal {i}", publication_date=datetime.now(),
                citation_count=i * 10, impact_factor=i * 2.0,
                score=100 - i * 2, rank=i + 1
            )
            scored_papers.append(paper)
            
        # Select top 10
        essential_papers = self.engine.select_essential_papers(scored_papers, 10)
        
        assert len(essential_papers) == 10
        assert all(p.rank <= 10 for p in essential_papers)
        
        # Should be in rank order
        for i, paper in enumerate(essential_papers):
            assert paper.rank == i + 1
            
    def test_get_score_breakdown(self):
        """Test score breakdown functionality."""
        scored_paper = ScoredPaper(
            pmid="12345", title="Test Paper", authors=["Author"],
            journal="Test Journal", publication_date=datetime(2023, 1, 1),
            citation_count=100, impact_factor=10.0, score=75.0, rank=1
        )
        
        breakdown = self.engine.get_score_breakdown(scored_paper, "test query")
        
        # Check all components are present
        required_keys = [
            'citation_score', 'impact_factor_score', 'recency_score',
            'relevance_score', 'weighted_citation', 'weighted_impact_factor',
            'weighted_recency', 'weighted_relevance', 'total_score'
        ]
        
        for key in required_keys:
            assert key in breakdown
            assert isinstance(breakdown[key], (int, float))
            
        # Check that weighted components sum to total (approximately)
        weighted_sum = (
            breakdown['weighted_citation'] +
            breakdown['weighted_impact_factor'] +
            breakdown['weighted_recency'] +
            breakdown['weighted_relevance']
        )
        
        assert abs(weighted_sum - breakdown['total_score']) < 0.1
        
    def test_update_weights(self):
        """Test weight updating functionality."""
        new_weights = ScoringWeights(
            citation_weight=0.6,
            impact_factor_weight=0.2,
            recency_weight=0.1,
            relevance_weight=0.1
        )
        
        self.engine.update_weights(new_weights)
        assert self.engine.weights.citation_weight == 0.6
        assert self.engine.weights.impact_factor_weight == 0.2
        
    def test_get_scoring_statistics(self):
        """Test scoring statistics calculation."""
        scored_papers = [
            ScoredPaper(
                pmid="1", title="Paper 1", authors=["A"], journal="J1",
                publication_date=datetime.now(), citation_count=100,
                impact_factor=10.0, score=85.0, rank=1
            ),
            ScoredPaper(
                pmid="2", title="Paper 2", authors=["B"], journal="J2",
                publication_date=datetime.now(), citation_count=0,
                impact_factor=0.0, score=45.0, rank=2
            )
        ]
        
        stats = self.engine.get_scoring_statistics(scored_papers)
        
        assert stats['total_papers'] == 2
        assert stats['average_score'] == 65.0
        assert stats['max_score'] == 85.0
        assert stats['min_score'] == 45.0
        assert stats['average_citations'] == 50.0
        assert stats['max_citations'] == 100
        assert stats['papers_with_citations'] == 1
        assert stats['papers_with_impact_factor'] == 1
        
    def test_extract_query_terms(self):
        """Test query term extraction."""
        query = "machine learning AND healthcare OR (artificial intelligence)"
        terms = self.engine._extract_query_terms(query)
        
        expected_terms = ["machine", "learning", "healthcare", "artificial", "intelligence"]
        assert all(term in terms for term in expected_terms)
        assert len(terms) == len(expected_terms)
        
    def test_normalize_journal_name(self):
        """Test journal name normalization."""
        # Test basic normalization
        assert self.engine._normalize_journal_name("Nature Medicine") == "nature medicine"
        
        # Test with prefixes/suffixes
        assert self.engine._normalize_journal_name("The Journal of Medicine (Online)") == "journal of medicine"
        
        # Test with extra spaces
        assert self.engine._normalize_journal_name("  Nature   Genetics  ") == "nature genetics"
        
    def test_edge_cases(self):
        """Test edge cases and error conditions."""
        # Empty papers list
        assert self.engine.rank_papers([]) == []
        assert self.engine.select_essential_papers([]) == []
        
        # Zero values
        score = self.engine.calculate_paper_score(self.sample_paper, 0, 0.0)
        assert 0 <= score <= 100
        
        # Negative values (should be handled gracefully)
        score = self.engine._calculate_citation_score(-10)
        assert score == 0.0
        
        score = self.engine._calculate_impact_factor_score(-5.0)
        assert score == 0.0