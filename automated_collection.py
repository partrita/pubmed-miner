#!/usr/bin/env python3
"""
Automated Essential Papers Collection Script

This script orchestrates the complete workflow for collecting, scoring, and 
publishing essential papers to GitHub Issues. It's designed to be run by 
GitHub Actions on a daily schedule.

Requirements addressed:
- 4.1: Read topic configurations and execute searches
- 4.2: Daily automated execution with proper error handling
- 4.3: GitHub Issues creation and updates
"""

import os
import sys
import logging
import traceback
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from pubmed_miner.utils.config_manager import ConfigurationManager
from pubmed_miner.services.paper_collection import PaperCollectionService
from pubmed_miner.services.paper_details import PaperDetailsService
from pubmed_miner.services.citation_service import CitationService
from pubmed_miner.services.impact_factor_service import ImpactFactorService
from pubmed_miner.scoring.engine import ScoringEngine
from pubmed_miner.services.github_manager import GitHubIssuesManager
from pubmed_miner.utils.change_tracker import ChangeTracker
from pubmed_miner.utils.error_handler import ErrorHandler, APIError, DataError
from pubmed_miner.models import TopicConfig, ScoredPaper


class AutomatedCollectionOrchestrator:
    """Main orchestrator for the automated paper collection workflow."""
    
    def __init__(self):
        """Initialize the orchestrator with all required services."""
        self.setup_logging()
        self.logger = logging.getLogger(__name__)
        
        try:
            # Initialize configuration manager
            self.config_manager = ConfigurationManager()
            self.system_config = self.config_manager.load_system_config()
            
            # Initialize services
            pubmed_email = self.config_manager.get_pubmed_email()
            self.paper_collector = PaperCollectionService(
                email=pubmed_email,
                rate_limit=3.0
            )
            self.paper_details = PaperDetailsService()
            self.citation_service = CitationService()
            self.impact_factor_service = ImpactFactorService()
            self.scoring_engine = ScoringEngine(self.system_config.scoring_weights)
            self.github_manager = GitHubIssuesManager(self.system_config.github)
            self.change_tracker = ChangeTracker()
            self.error_handler = ErrorHandler()
            
            self.logger.info("Automated collection orchestrator initialized successfully")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize orchestrator: {e}")
            raise
    
    def setup_logging(self) -> None:
        """Configure logging for the automation script."""
        # Create logs directory if it doesn't exist
        log_dir = Path("logs")
        log_dir.mkdir(exist_ok=True)
        
        # Configure logging
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        log_level = os.getenv('LOG_LEVEL', 'INFO').upper()
        
        logging.basicConfig(
            level=getattr(logging, log_level),
            format=log_format,
            handlers=[
                logging.FileHandler(log_dir / "automated_collection.log"),
                logging.StreamHandler(sys.stdout)
            ]
        )
    
    def run_complete_workflow(self) -> Dict[str, Any]:
        """Execute the complete automated workflow.
        
        Returns:
            Dictionary with execution results and statistics
        """
        start_time = datetime.now()
        results = {
            'start_time': start_time.isoformat(),
            'topics_processed': 0,
            'papers_collected': 0,
            'issues_created': 0,
            'issues_updated': 0,
            'errors': [],
            'success': False
        }
        
        try:
            self.logger.info("Starting automated essential papers collection workflow")
            
            # Load enabled topics
            enabled_topics = [topic for topic in self.system_config.topics if topic.enabled]
            self.logger.info(f"Found {len(enabled_topics)} enabled topics to process")
            
            if not enabled_topics:
                self.logger.warning("No enabled topics found, exiting")
                return results
            
            # Process each topic
            for topic in enabled_topics:
                try:
                    topic_result = self.process_topic(topic)
                    results['topics_processed'] += 1
                    results['papers_collected'] += topic_result.get('papers_collected', 0)
                    
                    if topic_result.get('issue_created'):
                        results['issues_created'] += 1
                    elif topic_result.get('issue_updated'):
                        results['issues_updated'] += 1
                        
                except Exception as e:
                    error_msg = f"Failed to process topic '{topic.name}': {str(e)}"
                    self.logger.error(error_msg)
                    results['errors'].append(error_msg)
                    
                    # Continue with other topics even if one fails
                    continue
            
            # Mark as successful if at least one topic was processed without critical errors
            results['success'] = results['topics_processed'] > 0
            
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            results['end_time'] = end_time.isoformat()
            results['duration_seconds'] = duration
            
            self.logger.info(f"Workflow completed in {duration:.2f} seconds")
            self.logger.info(f"Results: {results}")
            
            return results
            
        except Exception as e:
            error_msg = f"Critical error in workflow execution: {str(e)}"
            self.logger.error(error_msg)
            self.logger.error(traceback.format_exc())
            results['errors'].append(error_msg)
            results['success'] = False
            return results
    
    def process_topic(self, topic: TopicConfig) -> Dict[str, Any]:
        """Process a single topic through the complete pipeline.
        
        Args:
            topic: Topic configuration
            
        Returns:
            Dictionary with processing results
        """
        self.logger.info(f"Processing topic: {topic.name}")
        
        result = {
            'topic_name': topic.name,
            'papers_collected': 0,
            'essential_papers': 0,
            'issue_created': False,
            'issue_updated': False,
            'changes_detected': False
        }
        
        try:
            # Step 1: Search and collect papers
            self.logger.info(f"Searching PubMed for: {topic.query}")
            pmids = self.paper_collector.search_papers(topic.query, topic.max_papers)
            
            if not pmids:
                self.logger.warning(f"No papers found for topic: {topic.name}")
                return result
            
            self.logger.info(f"Found {len(pmids)} papers for topic: {topic.name}")
            result['papers_collected'] = len(pmids)
            
            # Step 2: Get detailed paper information
            self.logger.info("Fetching detailed paper information...")
            papers = self.paper_collector.get_paper_details(pmids)
            
            if not papers:
                self.logger.warning(f"No detailed paper information retrieved for topic: {topic.name}")
                return result
            
            # Step 3: Collect citation and impact factor data
            self.logger.info("Collecting citation and impact factor data...")
            scored_papers = []
            
            for paper in papers:
                try:
                    # Get citation count
                    citations = self.citation_service.get_citation_count(paper.pmid)
                    
                    # Get impact factor
                    impact_factor = self.impact_factor_service.get_impact_factor(paper.journal)
                    
                    # Calculate score
                    score = self.scoring_engine.calculate_paper_score(
                        paper, citations, impact_factor, topic.query
                    )
                    
                    scored_paper = ScoredPaper(
                        pmid=paper.pmid,
                        title=paper.title,
                        authors=paper.authors,
                        journal=paper.journal,
                        publication_date=paper.publication_date,
                        abstract=paper.abstract,
                        doi=paper.doi,
                        citation_count=citations,
                        impact_factor=impact_factor,
                        score=score
                    )
                    scored_papers.append(scored_paper)
                    
                except Exception as e:
                    self.logger.warning(f"Failed to score paper {paper.pmid}: {e}")
                    # Continue with other papers
                    continue
            
            if not scored_papers:
                self.logger.warning(f"No papers could be scored for topic: {topic.name}")
                return result
            
            # Step 4: Select essential papers
            essential_papers = self.scoring_engine.select_essential_papers(
                scored_papers, topic.essential_count
            )
            
            self.logger.info(f"Selected {len(essential_papers)} essential papers for topic: {topic.name}")
            result['essential_papers'] = len(essential_papers)
            
            # Step 5: Check for changes
            changes = self.change_tracker.detect_changes(topic.name, essential_papers)
            result['changes_detected'] = changes.has_changes()
            
            # Step 6: Create or update GitHub issue
            if essential_papers:
                issue_result = self.github_manager.create_or_update_issue(
                    topic.name, essential_papers
                )
                
                if issue_result.get('created'):
                    result['issue_created'] = True
                    self.logger.info(f"Created new issue for topic: {topic.name}")
                elif issue_result.get('updated'):
                    result['issue_updated'] = True
                    self.logger.info(f"Updated existing issue for topic: {topic.name}")
                
                # Update change tracking
                self.change_tracker.save_current_state(topic.name, essential_papers)
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error processing topic {topic.name}: {e}")
            raise


def main():
    """Main entry point for the automated collection script."""
    try:
        # Check if GitHub token is available (warn if not, but don't fail)
        if not os.getenv('GITHUB_TOKEN'):
            print("Warning: GITHUB_TOKEN not found - GitHub features will run in mock mode for local testing")
        
        # Initialize and run orchestrator
        orchestrator = AutomatedCollectionOrchestrator()
        results = orchestrator.run_complete_workflow()
        
        # Exit with appropriate code
        if results['success']:
            print("Automated collection completed successfully")
            sys.exit(0)
        else:
            print("Automated collection completed with errors")
            sys.exit(1)
            
    except Exception as e:
        print(f"Critical error: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()