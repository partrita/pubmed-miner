import pytest
import logging
import io
from src.pubmed_miner.models.config import GitHubConfig, SystemConfig, ScoringWeights, TopicConfig
from src.pubmed_miner.utils.logging_config import setup_logging

def test_github_token_not_in_repr():
    """Verify that GitHub token is not included in the repr of GitHubConfig."""
    secret = "ghp_VERY_SECRET_TOKEN"
    config = GitHubConfig(token=secret, repository="owner/repo")

    repr_str = repr(config)
    assert secret not in repr_str
    assert "token" not in repr_str
    assert "owner/repo" in repr_str

def test_github_token_not_in_system_config_repr():
    """Verify that GitHub token is not included in the repr of SystemConfig."""
    secret = "ghp_VERY_SECRET_TOKEN"
    config = GitHubConfig(token=secret, repository="owner/repo")
    system_config = SystemConfig(
        topics=[TopicConfig(name="test", query="test")],
        github=config,
        scoring_weights=ScoringWeights()
    )

    repr_str = repr(system_config)
    assert secret not in repr_str
    assert "token" not in repr_str

def test_github_token_not_in_logs():
    """Verify that GitHub token is not leaked in logs when logging the config object."""
    # Capture logs
    log_capture = io.StringIO()
    logger = logging.getLogger("test_security")
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(log_capture)
    logger.addHandler(handler)

    secret = "ghp_LOG_SECRET_TOKEN"
    config = GitHubConfig(token=secret, repository="owner/repo")

    logger.info("Config state: %s", config)

    log_output = log_capture.getvalue()
    assert secret not in log_output
    assert "token" not in log_output
