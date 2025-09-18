"""
Paper scoring and ranking engine.
"""
import math
import logging
from datetime import datetime, timedelta
from typing import List, Dict, Optional, Tuple
from dataclasses import replace

from ..models import Paper, ScoredPaper, ScoringWeights

logger = logging.getLogger(__name__)


class ScoringEngine:
    """Engine for calculating paper importance scores and ranking papers."""
    
    def __init__(self, weights: Optional[ScoringWeights] = None):
        """Initialize scoring engine.
        
        Args:
            weights: Scoring weights configuration
        """
        self.weights = weights or ScoringWeights()
        logger.info(f"Initialized ScoringEngine with weights: {self.weights}")
        
    def calculate_paper_score(
        self, 
        paper: Paper, 
        citations: int, 
        impact_factor: float,
        query: Optional[str] = None
    ) -> float:
        """Calculate comprehensive score for a paper.
        
        Args:
            paper: Paper object
            citations: Citation count
            impact_factor: Journal impact factor
            query: Original search query for relevance calculation
            
        Returns:
            Calculated score (0-100)
        """
        # Calculate individual component scores
        citation_score = self._calculate_citation_score(citations)
        impact_score = self._calculate_impact_factor_score(impact_factor)
        recency_score = self._calculate_recency_score(paper.publication_date)
        relevance_score = self._calculate_relevance_score(paper, query) if query else 50.0
        
        # Calculate weighted total score
        total_score = (
            citation_score * self.weights.citation_weight +
            impact_score * self.weights.impact_factor_weight +
            recency_score * self.weights.recency_weight +
            relevance_score * self.weights.relevance_weight
        )
        
        logger.debug(f"Score for {paper.pmid}: C:{citation_score:.1f} IF:{impact_score:.1f} "
                    f"R:{recency_score:.1f} Rel:{relevance_score:.1f} Total:{total_score:.1f}")
        
        return min(100.0, max(0.0, total_score))
        
    def rank_papers(self, papers: List[ScoredPaper]) -> List[ScoredPaper]:
        """Rank papers by score and assign rank numbers.
        
        Args:
            papers: List of scored papers
            
        Returns:
            List of papers sorted by score with rank assigned
        """
        if not papers:
            return papers
            
        # Sort by score (descending)
        sorted_papers = sorted(papers, key=lambda p: p.score, reverse=True)
        
        # Assign ranks
        ranked_papers = []
        for i, paper in enumerate(sorted_papers):
            ranked_paper = replace(paper, rank=i + 1)
            ranked_papers.append(ranked_paper)
            
        logger.info(f"Ranked {len(ranked_papers)} papers")
        return ranked_papers
        
    def select_essential_papers(
        self, 
        papers: List[ScoredPaper], 
        count: int = 15
    ) -> List[ScoredPaper]:
        """Select top essential papers based on scores.
        
        Args:
            papers: List of scored and ranked papers
            count: Number of papers to select
            
        Returns:
            List of top essential papers
        """
        if not papers:
            return []
            
        # Ensure papers are ranked
        if not all(p.rank > 0 for p in papers):
            papers = self.rank_papers(papers)
            
        # Select top papers
        essential_papers = papers[:count]
        
        logger.info(f"Selected {len(essential_papers)} essential papers from {len(papers)} total")
        return essential_papers
        
    def score_paper_batch(
        self,
        papers: List[Paper],
        citations: Dict[str, int],
        impact_factors: Dict[str, float],
        query: Optional[str] = None
    ) -> List[ScoredPaper]:
        """Score a batch of papers efficiently.
        
        Args:
            papers: List of Paper objects
            citations: Dictionary mapping PMID to citation count
            impact_factors: Dictionary mapping journal name to impact factor
            query: Original search query
            
        Returns:
            List of ScoredPaper objects
        """
        scored_papers = []
        
        for paper in papers:
            try:
                # Get citation count
                citation_count = citations.get(paper.pmid, 0)
                
                # Get impact factor (try exact match first, then normalized)
                journal_if = impact_factors.get(paper.journal, 0.0)
                if journal_if == 0.0:
                    # Try normalized journal name
                    normalized_journal = self._normalize_journal_name(paper.journal)
                    journal_if = impact_factors.get(normalized_journal, 0.0)
                
                # Calculate score
                score = self.calculate_paper_score(
                    paper, citation_count, journal_if, query
                )
                
                # Create scored paper
                scored_paper = ScoredPaper(
                    pmid=paper.pmid,
                    title=paper.title,
                    authors=paper.authors,
                    journal=paper.journal,
                    publication_date=paper.publication_date,
                    abstract=paper.abstract,
                    doi=paper.doi,
                    citation_count=citation_count,
                    impact_factor=journal_if,
                    score=score,
                    rank=0  # Will be assigned during ranking
                )
                
                scored_papers.append(scored_paper)
                
            except Exception as e:
                logger.warning(f"Error scoring paper {paper.pmid}: {e}")
                continue
                
        logger.info(f"Scored {len(scored_papers)} papers successfully")
        return scored_papers
        
    def _calculate_citation_score(self, citations: int) -> float:
        """Calculate normalized citation score (0-100).
        
        Args:
            citations: Citation count
            
        Returns:
            Normalized citation score
        """
        if citations <= 0:
            return 0.0
            
        # Use logarithmic scaling to handle wide range of citation counts
        # Most papers have 0-100 citations, exceptional papers have 1000+
        max_citations = 1000  # Papers with 1000+ citations get max score
        
        if citations >= max_citations:
            return 100.0
            
        # Logarithmic scaling: log(citations + 1) / log(max_citations + 1) * 100
        score = (math.log(citations + 1) / math.log(max_citations + 1)) * 100
        return min(100.0, score)
        
    def _calculate_impact_factor_score(self, impact_factor: float) -> float:
        """Calculate normalized impact factor score (0-100).
        
        Args:
            impact_factor: Journal impact factor
            
        Returns:
            Normalized impact factor score
        """
        if impact_factor <= 0:
            return 0.0
            
        # Use logarithmic scaling for impact factors
        # Most journals have IF 0-10, top journals have IF 50+
        max_if = 50.0  # Journals with IF 50+ get max score
        
        if impact_factor >= max_if:
            return 100.0
            
        # Logarithmic scaling
        score = (math.log(impact_factor + 1) / math.log(max_if + 1)) * 100
        return min(100.0, score)
        
    def _calculate_recency_score(self, publication_date: datetime) -> float:
        """Calculate recency score based on publication date (0-100).
        
        Args:
            publication_date: Paper publication date
            
        Returns:
            Recency score (newer papers get higher scores)
        """
        current_date = datetime.now()
        age_days = (current_date - publication_date).days
        
        # Papers published in the last 5 years get higher scores
        max_age_days = 5 * 365  # 5 years
        
        if age_days <= 0:
            return 100.0  # Future papers (shouldn't happen, but handle gracefully)
        elif age_days >= max_age_days:
            return 10.0   # Minimum score for old papers
        else:
            # Linear decay from 100 to 10 over 5 years
            score = 100.0 - (90.0 * age_days / max_age_days)
            return max(10.0, score)
            
    def _calculate_relevance_score(self, paper: Paper, query: Optional[str]) -> float:
        """Calculate relevance score based on query match (0-100).
        
        Args:
            paper: Paper object
            query: Search query string
            
        Returns:
            Relevance score
        """
        if not query:
            return 50.0  # Neutral score if no query
            
        # Extract query terms
        query_terms = self._extract_query_terms(query)
        if not query_terms:
            return 50.0
            
        # Create searchable text
        searchable_text = self._create_searchable_text(paper)
        
        # Calculate term matches
        title_matches = self._count_term_matches(paper.title, query_terms)
        abstract_matches = self._count_term_matches(paper.abstract or "", query_terms)
        
        # Calculate weighted relevance score
        title_score = (title_matches / len(query_terms)) * 100 if query_terms else 0
        abstract_score = (abstract_matches / len(query_terms)) * 100 if query_terms else 0
        
        # Weighted combination (title is more important than abstract)
        relevance_score = (title_score * 0.7) + (abstract_score * 0.3)
        
        return min(100.0, relevance_score)
        
    def _extract_query_terms(self, query: str) -> List[str]:
        """Extract meaningful terms from search query.
        
        Args:
            query: Search query string
            
        Returns:
            List of query terms
        """
        import re
        
        # Remove PubMed operators
        operators = ['AND', 'OR', 'NOT', '(', ')', '[', ']', '"']
        clean_query = query
        for op in operators:
            clean_query = clean_query.replace(op, ' ')
            
        # Extract terms (minimum 3 characters)
        terms = []
        for term in clean_query.split():
            term = re.sub(r'[^\w\s]', '', term).strip().lower()
            if len(term) >= 3:
                terms.append(term)
                
        return terms
        
    def _create_searchable_text(self, paper: Paper) -> str:
        """Create searchable text from paper data.
        
        Args:
            paper: Paper object
            
        Returns:
            Combined searchable text
        """
        parts = [
            paper.title or "",
            paper.abstract or "",
            " ".join(paper.authors),
            paper.journal
        ]
        
        return " ".join(part for part in parts if part).lower()
        
    def _count_term_matches(self, text: str, terms: List[str]) -> int:
        """Count how many query terms appear in text.
        
        Args:
            text: Text to search
            terms: List of terms to find
            
        Returns:
            Number of matching terms
        """
        if not text or not terms:
            return 0
            
        text_lower = text.lower()
        matches = 0
        
        for term in terms:
            if term in text_lower:
                matches += 1
                
        return matches
        
    def _normalize_journal_name(self, journal_name: str) -> str:
        """Normalize journal name for matching.
        
        Args:
            journal_name: Raw journal name
            
        Returns:
            Normalized journal name
        """
        import re
        
        if not journal_name:
            return ""
            
        # Convert to lowercase and clean
        normalized = journal_name.lower().strip()
        
        # Remove common patterns
        patterns = [
            r'^the\s+',
            r'\s*\(online\)$',
            r'\s*\(print\)$',
            r'\s*:\s*official.*$'
        ]
        
        for pattern in patterns:
            normalized = re.sub(pattern, '', normalized, flags=re.IGNORECASE)
            
        # Clean up spaces
        normalized = re.sub(r'\s+', ' ', normalized).strip()
        
        return normalized
        
    def get_score_breakdown(self, scored_paper: ScoredPaper, query: Optional[str] = None) -> Dict[str, float]:
        """Get detailed score breakdown for a paper.
        
        Args:
            scored_paper: ScoredPaper object
            query: Original search query
            
        Returns:
            Dictionary with score component breakdown
        """
        citation_score = self._calculate_citation_score(scored_paper.citation_count)
        impact_score = self._calculate_impact_factor_score(scored_paper.impact_factor)
        recency_score = self._calculate_recency_score(scored_paper.publication_date)
        relevance_score = self._calculate_relevance_score(scored_paper, query) if query else 50.0
        
        return {
            'citation_score': citation_score,
            'impact_factor_score': impact_score,
            'recency_score': recency_score,
            'relevance_score': relevance_score,
            'weighted_citation': citation_score * self.weights.citation_weight,
            'weighted_impact_factor': impact_score * self.weights.impact_factor_weight,
            'weighted_recency': recency_score * self.weights.recency_weight,
            'weighted_relevance': relevance_score * self.weights.relevance_weight,
            'total_score': scored_paper.score
        }
        
    def update_weights(self, new_weights: ScoringWeights) -> None:
        """Update scoring weights.
        
        Args:
            new_weights: New scoring weights
        """
        self.weights = new_weights
        logger.info(f"Updated scoring weights: {self.weights}")
        
    def get_scoring_statistics(self, scored_papers: List[ScoredPaper]) -> Dict[str, float]:
        """Get statistics about scored papers.
        
        Args:
            scored_papers: List of scored papers
            
        Returns:
            Dictionary with scoring statistics
        """
        if not scored_papers:
            return {}
            
        scores = [p.score for p in scored_papers]
        citations = [p.citation_count for p in scored_papers]
        impact_factors = [p.impact_factor for p in scored_papers if p.impact_factor > 0]
        
        stats = {
            'total_papers': len(scored_papers),
            'average_score': sum(scores) / len(scores),
            'max_score': max(scores),
            'min_score': min(scores),
            'average_citations': sum(citations) / len(citations) if citations else 0,
            'max_citations': max(citations) if citations else 0,
            'average_impact_factor': sum(impact_factors) / len(impact_factors) if impact_factors else 0,
            'max_impact_factor': max(impact_factors) if impact_factors else 0,
            'papers_with_impact_factor': len(impact_factors),
            'papers_with_citations': len([c for c in citations if c > 0])
        }
        
        return stats