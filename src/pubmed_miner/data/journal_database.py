"""
Extended journal database with comprehensive impact factor data.
"""
import json
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class JournalDatabase:
    """Extended database of journal impact factors and metadata."""
    
    def __init__(self):
        """Initialize journal database with comprehensive data."""
        self.journals = self._load_comprehensive_journal_data()
        self.abbreviations = self._load_journal_abbreviations()
        
        logger.info(f"Loaded {len(self.journals)} journals and {len(self.abbreviations)} abbreviations")
        
    def _load_comprehensive_journal_data(self) -> Dict[str, Dict]:
        """Load comprehensive journal data including impact factors, categories, and metadata.
        
        Returns:
            Dictionary mapping normalized journal names to journal data
        """
        # Comprehensive journal database with 2023 impact factors
        journal_data = {
            # Top-tier journals
            'nature': {
                'impact_factor': 64.8,
                'category': 'multidisciplinary',
                'publisher': 'Nature Publishing Group',
                'issn': '0028-0836',
                'country': 'UK',
                'open_access': False
            },
            'science': {
                'impact_factor': 56.9,
                'category': 'multidisciplinary',
                'publisher': 'American Association for the Advancement of Science',
                'issn': '0036-8075',
                'country': 'USA',
                'open_access': False
            },
            'cell': {
                'impact_factor': 64.5,
                'category': 'cell biology',
                'publisher': 'Cell Press',
                'issn': '0092-8674',
                'country': 'USA',
                'open_access': False
            },
            'new england journal of medicine': {
                'impact_factor': 176.1,
                'category': 'medicine',
                'publisher': 'Massachusetts Medical Society',
                'issn': '0028-4793',
                'country': 'USA',
                'open_access': False
            },
            'lancet': {
                'impact_factor': 168.9,
                'category': 'medicine',
                'publisher': 'Elsevier',
                'issn': '0140-6736',
                'country': 'UK',
                'open_access': False
            },
            
            # Nature family journals
            'nature medicine': {
                'impact_factor': 87.2,
                'category': 'medicine',
                'publisher': 'Nature Publishing Group',
                'issn': '1078-8956',
                'country': 'UK',
                'open_access': False
            },
            'nature genetics': {
                'impact_factor': 41.3,
                'category': 'genetics',
                'publisher': 'Nature Publishing Group',
                'issn': '1061-4036',
                'country': 'UK',
                'open_access': False
            },
            'nature biotechnology': {
                'impact_factor': 68.2,
                'category': 'biotechnology',
                'publisher': 'Nature Publishing Group',
                'issn': '1087-0156',
                'country': 'UK',
                'open_access': False
            },
            'nature communications': {
                'impact_factor': 17.7,
                'category': 'multidisciplinary',
                'publisher': 'Nature Publishing Group',
                'issn': '2041-1723',
                'country': 'UK',
                'open_access': True
            },
            'nature methods': {
                'impact_factor': 47.9,
                'category': 'biochemistry',
                'publisher': 'Nature Publishing Group',
                'issn': '1548-7091',
                'country': 'UK',
                'open_access': False
            },
            
            # High-impact specialized journals
            'proceedings of the national academy of sciences': {
                'impact_factor': 12.8,
                'category': 'multidisciplinary',
                'publisher': 'National Academy of Sciences',
                'issn': '0027-8424',
                'country': 'USA',
                'open_access': False
            },
            'journal of biological chemistry': {
                'impact_factor': 4.8,
                'category': 'biochemistry',
                'publisher': 'American Society for Biochemistry and Molecular Biology',
                'issn': '0021-9258',
                'country': 'USA',
                'open_access': True
            },
            'nucleic acids research': {
                'impact_factor': 14.9,
                'category': 'biochemistry',
                'publisher': 'Oxford University Press',
                'issn': '0305-1048',
                'country': 'UK',
                'open_access': True
            },
            'genome biology': {
                'impact_factor': 17.9,
                'category': 'genomics',
                'publisher': 'BioMed Central',
                'issn': '1474-760X',
                'country': 'UK',
                'open_access': True
            },
            'genome research': {
                'impact_factor': 7.0,
                'category': 'genomics',
                'publisher': 'Cold Spring Harbor Laboratory Press',
                'issn': '1088-9051',
                'country': 'USA',
                'open_access': True
            },
            
            # Bioinformatics and computational biology
            'bioinformatics': {
                'impact_factor': 6.9,
                'category': 'bioinformatics',
                'publisher': 'Oxford University Press',
                'issn': '1367-4803',
                'country': 'UK',
                'open_access': False
            },
            'bmc bioinformatics': {
                'impact_factor': 3.3,
                'category': 'bioinformatics',
                'publisher': 'BioMed Central',
                'issn': '1471-2105',
                'country': 'UK',
                'open_access': True
            },
            'plos computational biology': {
                'impact_factor': 4.3,
                'category': 'computational biology',
                'publisher': 'Public Library of Science',
                'issn': '1553-7358',
                'country': 'USA',
                'open_access': True
            },
            'molecular biology and evolution': {
                'impact_factor': 16.2,
                'category': 'evolutionary biology',
                'publisher': 'Oxford University Press',
                'issn': '0737-4038',
                'country': 'UK',
                'open_access': False
            },
            
            # PLOS journals
            'plos one': {
                'impact_factor': 3.7,
                'category': 'multidisciplinary',
                'publisher': 'Public Library of Science',
                'issn': '1932-6203',
                'country': 'USA',
                'open_access': True
            },
            'plos biology': {
                'impact_factor': 9.8,
                'category': 'biology',
                'publisher': 'Public Library of Science',
                'issn': '1544-9173',
                'country': 'USA',
                'open_access': True
            },
            'plos medicine': {
                'impact_factor': 13.8,
                'category': 'medicine',
                'publisher': 'Public Library of Science',
                'issn': '1549-1277',
                'country': 'USA',
                'open_access': True
            },
            'plos genetics': {
                'impact_factor': 4.5,
                'category': 'genetics',
                'publisher': 'Public Library of Science',
                'issn': '1553-7404',
                'country': 'USA',
                'open_access': True
            },
            
            # Other major journals
            'scientific reports': {
                'impact_factor': 4.6,
                'category': 'multidisciplinary',
                'publisher': 'Nature Publishing Group',
                'issn': '2045-2322',
                'country': 'UK',
                'open_access': True
            },
            'elife': {
                'impact_factor': 8.7,
                'category': 'life sciences',
                'publisher': 'eLife Sciences Publications',
                'issn': '2050-084X',
                'country': 'UK',
                'open_access': True
            },
            'current biology': {
                'impact_factor': 10.9,
                'category': 'biology',
                'publisher': 'Cell Press',
                'issn': '0960-9822',
                'country': 'USA',
                'open_access': False
            },
            'molecular cell': {
                'impact_factor': 16.0,
                'category': 'cell biology',
                'publisher': 'Cell Press',
                'issn': '1097-2765',
                'country': 'USA',
                'open_access': False
            },
            
            # Medical journals
            'jama': {
                'impact_factor': 157.3,
                'category': 'medicine',
                'publisher': 'American Medical Association',
                'issn': '0098-7484',
                'country': 'USA',
                'open_access': False
            },
            'british medical journal': {
                'impact_factor': 105.7,
                'category': 'medicine',
                'publisher': 'BMJ Publishing Group',
                'issn': '0959-8138',
                'country': 'UK',
                'open_access': False
            },
            'annals of internal medicine': {
                'impact_factor': 51.6,
                'category': 'medicine',
                'publisher': 'American College of Physicians',
                'issn': '0003-4819',
                'country': 'USA',
                'open_access': False
            },
            
            # Specialized high-impact journals
            'immunity': {
                'impact_factor': 43.5,
                'category': 'immunology',
                'publisher': 'Cell Press',
                'issn': '1074-7613',
                'country': 'USA',
                'open_access': False
            },
            'cancer cell': {
                'impact_factor': 38.5,
                'category': 'oncology',
                'publisher': 'Cell Press',
                'issn': '1535-6108',
                'country': 'USA',
                'open_access': False
            },
            'neuron': {
                'impact_factor': 16.2,
                'category': 'neuroscience',
                'publisher': 'Cell Press',
                'issn': '0896-6273',
                'country': 'USA',
                'open_access': False
            }
        }
        
        return journal_data
        
    def _load_journal_abbreviations(self) -> Dict[str, str]:
        """Load common journal abbreviations and their full names.
        
        Returns:
            Dictionary mapping abbreviations to full journal names
        """
        abbreviations = {
            # Nature family
            'nat': 'nature',
            'nat med': 'nature medicine',
            'nat genet': 'nature genetics',
            'nat biotechnol': 'nature biotechnology',
            'nat commun': 'nature communications',
            'nat methods': 'nature methods',
            
            # Common abbreviations
            'nejm': 'new england journal of medicine',
            'jama': 'jama',
            'bmj': 'british medical journal',
            'pnas': 'proceedings of the national academy of sciences',
            'jbc': 'journal of biological chemistry',
            'nar': 'nucleic acids research',
            'curr biol': 'current biology',
            'mol cell': 'molecular cell',
            'sci rep': 'scientific reports',
            
            # PLOS abbreviations
            'plos one': 'plos one',
            'plos biol': 'plos biology',
            'plos med': 'plos medicine',
            'plos genet': 'plos genetics',
            'plos comput biol': 'plos computational biology',
            
            # Bioinformatics
            'bioinformatics': 'bioinformatics',
            'bmc bioinformatics': 'bmc bioinformatics',
            'genome res': 'genome research',
            'genome biol': 'genome biology',
            'mol biol evol': 'molecular biology and evolution',
            
            # Medical abbreviations
            'ann intern med': 'annals of internal medicine',
            'j clin invest': 'journal of clinical investigation',
            'blood': 'blood',
            'cancer res': 'cancer research',
            'j immunol': 'journal of immunology'
        }
        
        return abbreviations
        
    def get_journal_info(self, journal_name: str) -> Optional[Dict]:
        """Get comprehensive journal information.
        
        Args:
            journal_name: Journal name (can be abbreviated)
            
        Returns:
            Dictionary with journal information or None
        """
        normalized_name = self._normalize_name(journal_name)
        
        # Direct lookup
        if normalized_name in self.journals:
            return self.journals[normalized_name]
            
        # Try abbreviation lookup
        if normalized_name in self.abbreviations:
            full_name = self.abbreviations[normalized_name]
            if full_name in self.journals:
                return self.journals[full_name]
                
        # Fuzzy matching
        best_match = self._find_best_match(normalized_name)
        if best_match:
            return self.journals[best_match]
            
        return None
        
    def get_impact_factor(self, journal_name: str) -> Optional[float]:
        """Get impact factor for a journal.
        
        Args:
            journal_name: Journal name
            
        Returns:
            Impact factor or None
        """
        journal_info = self.get_journal_info(journal_name)
        return journal_info.get('impact_factor') if journal_info else None
        
    def get_journals_by_category(self, category: str) -> List[Tuple[str, Dict]]:
        """Get all journals in a specific category.
        
        Args:
            category: Journal category
            
        Returns:
            List of (journal_name, journal_info) tuples
        """
        category_lower = category.lower()
        matching_journals = []
        
        for journal_name, journal_info in self.journals.items():
            if journal_info.get('category', '').lower() == category_lower:
                matching_journals.append((journal_name, journal_info))
                
        # Sort by impact factor (descending)
        matching_journals.sort(key=lambda x: x[1].get('impact_factor', 0), reverse=True)
        return matching_journals
        
    def get_open_access_journals(self) -> List[Tuple[str, Dict]]:
        """Get all open access journals.
        
        Returns:
            List of (journal_name, journal_info) tuples for open access journals
        """
        open_access_journals = []
        
        for journal_name, journal_info in self.journals.items():
            if journal_info.get('open_access', False):
                open_access_journals.append((journal_name, journal_info))
                
        # Sort by impact factor (descending)
        open_access_journals.sort(key=lambda x: x[1].get('impact_factor', 0), reverse=True)
        return open_access_journals
        
    def get_top_journals(self, limit: int = 10, category: Optional[str] = None) -> List[Tuple[str, Dict]]:
        """Get top journals by impact factor.
        
        Args:
            limit: Maximum number of journals to return
            category: Optional category filter
            
        Returns:
            List of (journal_name, journal_info) tuples
        """
        if category:
            journals = self.get_journals_by_category(category)
        else:
            journals = [(name, info) for name, info in self.journals.items()]
            journals.sort(key=lambda x: x[1].get('impact_factor', 0), reverse=True)
            
        return journals[:limit]
        
    def search_journals(self, query: str, limit: int = 10) -> List[Tuple[str, Dict, float]]:
        """Search journals by name with similarity scores.
        
        Args:
            query: Search query
            limit: Maximum number of results
            
        Returns:
            List of (journal_name, journal_info, similarity_score) tuples
        """
        from difflib import SequenceMatcher
        
        normalized_query = self._normalize_name(query)
        results = []
        
        for journal_name, journal_info in self.journals.items():
            similarity = SequenceMatcher(None, normalized_query, journal_name).ratio()
            if similarity > 0.3:  # Minimum similarity threshold
                results.append((journal_name, journal_info, similarity))
                
        # Sort by similarity (descending)
        results.sort(key=lambda x: x[2], reverse=True)
        return results[:limit]
        
    def add_journal(self, name: str, impact_factor: float, category: str = 'unknown', 
                   publisher: str = 'unknown', **kwargs) -> None:
        """Add a new journal to the database.
        
        Args:
            name: Journal name
            impact_factor: Impact factor
            category: Journal category
            publisher: Publisher name
            **kwargs: Additional journal metadata
        """
        normalized_name = self._normalize_name(name)
        
        journal_info = {
            'impact_factor': impact_factor,
            'category': category,
            'publisher': publisher,
            **kwargs
        }
        
        self.journals[normalized_name] = journal_info
        logger.info(f"Added journal: {name} (IF: {impact_factor})")
        
    def export_to_csv(self, filepath: str) -> None:
        """Export journal database to CSV file.
        
        Args:
            filepath: Output CSV file path
        """
        with open(filepath, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['name', 'impact_factor', 'category', 'publisher', 'issn', 'country', 'open_access']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for journal_name, journal_info in self.journals.items():
                row = {'name': journal_name}
                row.update(journal_info)
                writer.writerow(row)
                
        logger.info(f"Exported {len(self.journals)} journals to {filepath}")
        
    def _normalize_name(self, name: str) -> str:
        """Normalize journal name for consistent lookup."""
        if not name:
            return ""
            
        # Convert to lowercase and strip
        normalized = name.lower().strip()
        
        # Remove common prefixes and suffixes
        patterns_to_remove = [
            r'^the\s+',
            r'\s*\(online\)$',
            r'\s*\(print\)$',
            r'\s*:\s*official.*$'
        ]
        
        import re
        for pattern in patterns_to_remove:
            normalized = re.sub(pattern, '', normalized, flags=re.IGNORECASE)
            
        # Clean up spaces
        normalized = re.sub(r'\s+', ' ', normalized).strip()
        
        return normalized
        
    def _find_best_match(self, normalized_name: str, threshold: float = 0.6) -> Optional[str]:
        """Find best matching journal name using fuzzy matching.
        
        Args:
            normalized_name: Normalized journal name to match
            threshold: Minimum similarity threshold
            
        Returns:
            Best matching journal name or None
        """
        from difflib import SequenceMatcher
        
        best_match = None
        best_score = 0.0
        
        for journal_name in self.journals.keys():
            score = SequenceMatcher(None, normalized_name, journal_name).ratio()
            if score > best_score and score >= threshold:
                best_score = score
                best_match = journal_name
                
        return best_match
        
    def get_statistics(self) -> Dict[str, any]:
        """Get database statistics.
        
        Returns:
            Dictionary with database statistics
        """
        categories = {}
        publishers = {}
        open_access_count = 0
        impact_factors = []
        
        for journal_info in self.journals.values():
            # Count categories
            category = journal_info.get('category', 'unknown')
            categories[category] = categories.get(category, 0) + 1
            
            # Count publishers
            publisher = journal_info.get('publisher', 'unknown')
            publishers[publisher] = publishers.get(publisher, 0) + 1
            
            # Count open access
            if journal_info.get('open_access', False):
                open_access_count += 1
                
            # Collect impact factors
            if 'impact_factor' in journal_info:
                impact_factors.append(journal_info['impact_factor'])
                
        stats = {
            'total_journals': len(self.journals),
            'total_abbreviations': len(self.abbreviations),
            'categories': len(categories),
            'publishers': len(publishers),
            'open_access_journals': open_access_count,
            'average_impact_factor': sum(impact_factors) / len(impact_factors) if impact_factors else 0,
            'max_impact_factor': max(impact_factors) if impact_factors else 0,
            'min_impact_factor': min(impact_factors) if impact_factors else 0,
            'top_categories': sorted(categories.items(), key=lambda x: x[1], reverse=True)[:5],
            'top_publishers': sorted(publishers.items(), key=lambda x: x[1], reverse=True)[:5]
        }
        
        return stats