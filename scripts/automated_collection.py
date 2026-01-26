#!/usr/bin/env python3
"""
Automated Essential Papers Collection Script

This script orchestrates the complete workflow for collecting, scoring, and
publishing essential papers to MdBook. It's designed to be run by
GitHub Actions on a daily schedule.

Requirements addressed:
- 4.1: Read topic configurations and execute searches
- 4.2: Daily automated execution with proper error handling
- 4.3: MdBook page creation and updates
"""

import os
import sys
import logging
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

# Add src to path for imports
# Path(__file__).parent is scripts/
# Path(__file__).parent.parent is root/
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pubmed_miner.utils.config_manager import ConfigurationManager
from pubmed_miner.utils.csv_manager import CSVManager
from pubmed_miner.services.paper_collection import PaperCollectionService
from pubmed_miner.services.paper_details import PaperDetailsService
from pubmed_miner.services.citation_service import CitationService
from pubmed_miner.services.impact_factor_service import ImpactFactorService
from pubmed_miner.scoring.engine import ScoringEngine
from pubmed_miner.services.mdbook_manager import MdBookManager
from pubmed_miner.utils.change_tracker import ChangeTracker
from pubmed_miner.utils.error_handler import ErrorHandler
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
                email=pubmed_email, rate_limit=3.0
            )
            self.paper_details = PaperDetailsService()
            self.citation_service = CitationService()
            self.impact_factor_service = ImpactFactorService()
            self.scoring_engine = ScoringEngine(self.system_config.scoring_weights)
            self.mdbook_manager = MdBookManager()
            self.change_tracker = ChangeTracker()
            self.error_handler = ErrorHandler()

            self.logger.info(
                "Automated collection orchestrator initialized successfully"
            )

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
        log_level = os.getenv("LOG_LEVEL", "INFO").upper()

        logging.basicConfig(
            level=getattr(logging, log_level),
            format=log_format,
            handlers=[
                logging.FileHandler(log_dir / "automated_collection.log"),
                logging.StreamHandler(sys.stdout),
            ],
        )

    def run_complete_workflow(self) -> Dict[str, Any]:
        """Execute the complete automated workflow.

        Returns:
            Dictionary with execution results and statistics
        """
        start_time = datetime.now()
        results = {
            "start_time": start_time.isoformat(),
            "topics_processed": 0,
            "papers_collected": 0,
            "pages_created": 0,
            "errors": [],
            "success": False,
        }

        try:
            self.logger.info("Starting automated essential papers collection workflow")

            # Load enabled topics
            enabled_topics = [
                topic for topic in self.system_config.topics if topic.enabled
            ]
            self.logger.info(f"Found {len(enabled_topics)} enabled topics to process")

            if not enabled_topics:
                self.logger.warning("No enabled topics found, exiting")
                return results

            # Process each topic
            for topic in enabled_topics:
                try:
                    topic_result = self.process_topic(topic)
                    results["topics_processed"] += 1
                    results["papers_collected"] += topic_result.get(
                        "papers_collected", 0
                    )

                    if topic_result.get("page_created"):
                        results["pages_created"] += 1

                except Exception as e:
                    error_msg = f"Failed to process topic '{topic.name}': {str(e)}"
                    self.logger.error(error_msg)
                    results["errors"].append(error_msg)

                    # Continue with other topics even if one fails
                    continue

            # Mark as successful if at least one topic was processed without critical errors
            results["success"] = results["topics_processed"] > 0

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            results["end_time"] = end_time.isoformat()
            results["duration_seconds"] = duration

            self.logger.info(f"Workflow completed in {duration:.2f} seconds")
            self.logger.info(f"Results: {results}")

            return results

        except Exception as e:
            error_msg = f"Critical error in workflow execution: {str(e)}"
            self.logger.error(error_msg)
            self.logger.error(traceback.format_exc())
            results["errors"].append(error_msg)
            results["success"] = False
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
            "topic_name": topic.name,
            "papers_collected": 0,
            "essential_papers": 0,
            "page_created": False,
            "changes_detected": False,
        }

        try:
            # Step 1: Search and collect papers
            self.logger.info(f"Searching PubMed for: {topic.query}")
            pmids = self.paper_collector.search_papers(topic.query, topic.max_papers)

            if not pmids:
                self.logger.warning(f"No papers found for topic: {topic.name}")
                return result

            self.logger.info(f"Found {len(pmids)} papers for topic: {topic.name}")
            result["papers_collected"] = len(pmids)

            # Step 2: Get detailed paper information
            self.logger.info("Fetching detailed paper information...")
            papers = self.paper_collector.get_paper_details(pmids)

            if not papers:
                self.logger.warning(
                    f"No detailed paper information retrieved for topic: {topic.name}"
                )
                return result

            # Step 3: Collect citation and impact factor data
            self.logger.info("Collecting citation and impact factor data...")
            scored_papers = []

            for paper in papers:
                try:
                    # Get citation count (safely handle None values)
                    citations = self.citation_service.get_citation_count(paper.pmid)
                    safe_citations = self._safe_int_value(citations, default=0)

                    # Get impact factor (safely handle None values)
                    impact_factor = self.impact_factor_service.get_impact_factor(
                        paper.journal
                    )
                    safe_impact_factor = self._safe_float_value(
                        impact_factor, default=0.0
                    )

                    # Calculate score
                    score = self.scoring_engine.calculate_paper_score(
                        paper, safe_citations, safe_impact_factor, topic.query
                    )

                    scored_paper = ScoredPaper(
                        pmid=paper.pmid,
                        title=paper.title,
                        authors=paper.authors,
                        journal=paper.journal,
                        publication_date=paper.publication_date,
                        abstract=paper.abstract,
                        doi=paper.doi,
                        topic=topic.name,
                        citation_count=safe_citations,
                        impact_factor=safe_impact_factor,
                        score=score,
                    )
                    scored_papers.append(scored_paper)

                except Exception as e:
                    self.logger.warning(f"Failed to score paper {paper.pmid}: {e}")
                    # Continue with other papers
                    continue

            if not scored_papers:
                self.logger.warning(
                    f"No papers could be scored for topic: {topic.name}"
                )
                return result

            # Step 4: Select essential papers
            essential_papers = self.scoring_engine.select_essential_papers(
                scored_papers, topic.essential_count
            )

            self.logger.info(
                f"Selected {len(essential_papers)} essential papers for topic: {topic.name}"
            )
            result["essential_papers"] = len(essential_papers)

            # Save to CSV
            if essential_papers:
                try:
                    CSVManager.save_papers(
                        essential_papers,
                        "data/collections.csv",
                        include_scoring=True,
                        append=True
                    )
                    self.logger.info(f"Saved {len(essential_papers)} papers to CSV")
                except Exception as e:
                    self.logger.error(f"Failed to save papers to CSV: {e}")

            # Step 5: Check for changes
            changes = self.change_tracker.detect_changes(topic.name, essential_papers)
            result["changes_detected"] = changes["has_changes"]

            # Step 6: Create MdBook page
            if essential_papers:
                relative_path = self.mdbook_manager.create_daily_page(
                    topic.name, essential_papers
                )
                self.mdbook_manager.update_summary(relative_path, topic.name)
                
                result["page_created"] = True
                self.logger.info(f"Created mdbook page for topic: {topic.name}")

                # Update change tracking
                self.change_tracker.save_paper_snapshot(topic.name, essential_papers)

            return result

        except Exception as e:
            self.logger.error(f"Error processing topic {topic.name}: {e}")
            raise

    def _safe_int_value(self, value, default=0):
        """Safely convert value to integer, handling None and invalid values.

        Args:
            value: Value to convert
            default: Default value if conversion fails

        Returns:
            Safe integer value
        """
        if value is None:
            return default

        try:
            int_value = int(value)
            return max(0, int_value)  # Ensure non-negative
        except (ValueError, TypeError):
            self.logger.debug(
                f"Failed to convert {value} to int, using default {default}"
            )
            return default

    def _safe_float_value(self, value, default=0.0):
        """Safely convert value to float, handling None and invalid values.

        Args:
            value: Value to convert
            default: Default value if conversion fails

        Returns:
            Safe float value
        """
        if value is None:
            return default

        try:
            float_value = float(value)
            return max(0.0, float_value)  # Ensure non-negative
        except (ValueError, TypeError):
            self.logger.debug(
                f"Failed to convert {value} to float, using default {default}"
            )
            return default


def main():
    """Main entry point for the automated collection script."""
    try:
        # Initialize and run orchestrator
        orchestrator = AutomatedCollectionOrchestrator()
        results = orchestrator.run_complete_workflow()

        # Exit with appropriate code
        if results["success"]:
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