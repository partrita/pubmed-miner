"""
PubMed Miner - Essential Papers Recommender

A tool for automatically collecting and ranking essential papers from PubMed
based on citation counts, journal impact factors, and other metrics.
"""

__version__ = "0.2.0"
__author__ = "taeyoon kim"

# Import main models for easy access
from .models import (
    GitHubConfig,
    Paper,
    ScoredPaper,
    ScoringWeights,
    SystemConfig,
    TopicConfig,
)

__all__ = [
    "GitHubConfig",
    "Paper",
    "ScoredPaper",
    "ScoringWeights",
    "SystemConfig",
    "TopicConfig",
]
