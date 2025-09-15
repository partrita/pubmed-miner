"""
Pytest configuration and shared fixtures.
"""
import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    
    # Cleanup
    import shutil
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def sample_config_dir(temp_dir):
    """Create a temporary directory with sample configuration files."""
    config_dir = Path(temp_dir) / "config"
    config_dir.mkdir()
    
    # Create sample topics.yaml
    topics_content = """
topics:
  - name: "test-topic"
    query: "test query"
    max_papers: 100
    essential_count: 10
    enabled: true
  - name: "disabled-topic"
    query: "disabled query"
    max_papers: 50
    essential_count: 5
    enabled: false
"""
    
    with open(config_dir / "topics.yaml", "w") as f:
        f.write(topics_content)
    
    # Create sample settings.yaml
    settings_content = """
github:
  repository: "test/repo"
  issue_labels:
    - "test"
    - "automated"

scoring_weights:
  citation_weight: 0.4
  impact_factor_weight: 0.3
  recency_weight: 0.2
  relevance_weight: 0.1

cache_settings:
  citation_cache_days: 7
  impact_factor_cache_days: 365
  paper_metadata_cache_days: 30
"""
    
    with open(config_dir / "settings.yaml", "w") as f:
        f.write(settings_content)
    
    return config_dir


@pytest.fixture
def mock_github_token():
    """Mock GitHub token environment variable."""
    original_token = os.environ.get('GITHUB_TOKEN')
    os.environ['GITHUB_TOKEN'] = 'mock_token_for_local_testing'
    
    yield 'mock_token_for_local_testing'
    
    # Restore original value
    if original_token:
        os.environ['GITHUB_TOKEN'] = original_token
    else:
        os.environ.pop('GITHUB_TOKEN', None)


@pytest.fixture
def no_github_token():
    """Remove GitHub token for testing without token."""
    original_token = os.environ.get('GITHUB_TOKEN')
    if 'GITHUB_TOKEN' in os.environ:
        del os.environ['GITHUB_TOKEN']
    
    yield
    
    # Restore original value
    if original_token:
        os.environ['GITHUB_TOKEN'] = original_token


@pytest.fixture
def sample_papers():
    """Create sample Paper objects for testing."""
    from datetime import datetime
    from src.pubmed_miner.models import Paper
    
    papers = [
        Paper(
            pmid="12345",
            title="Sample Paper 1",
            authors=["John Doe", "Jane Smith"],
            journal="Nature",
            publication_date=datetime(2023, 1, 15),
            abstract="This is a sample abstract for paper 1.",
            doi="10.1038/nature.2023.12345"
        ),
        Paper(
            pmid="67890",
            title="Sample Paper 2",
            authors=["Bob Johnson"],
            journal="Science",
            publication_date=datetime(2023, 2, 10),
            abstract="This is a sample abstract for paper 2.",
            doi="10.1126/science.2023.67890"
        ),
        Paper(
            pmid="11111",
            title="Sample Paper 3",
            authors=["Alice Brown", "Charlie Wilson"],
            journal="Cell",
            publication_date=datetime(2023, 3, 5),
            abstract="This is a sample abstract for paper 3.",
            doi="10.1016/j.cell.2023.11111"
        )
    ]
    
    return papers


@pytest.fixture
def sample_scored_papers(sample_papers):
    """Create sample ScoredPaper objects for testing."""
    from src.pubmed_miner.models import ScoredPaper
    
    scored_papers = []
    for i, paper in enumerate(sample_papers):
        scored_paper = ScoredPaper(
            pmid=paper.pmid,
            title=paper.title,
            authors=paper.authors,
            journal=paper.journal,
            publication_date=paper.publication_date,
            abstract=paper.abstract,
            doi=paper.doi,
            citation_count=(i + 1) * 50,
            impact_factor=(i + 1) * 10.0,
            score=90.0 - (i * 5),
            rank=i + 1
        )
        scored_papers.append(scored_paper)
    
    return scored_papers


@pytest.fixture
def mock_requests():
    """Mock requests library for API testing."""
    with pytest.mock.patch('requests.get') as mock_get:
        with pytest.mock.patch('requests.post') as mock_post:
            with pytest.mock.patch('requests.patch') as mock_patch:
                yield {
                    'get': mock_get,
                    'post': mock_post,
                    'patch': mock_patch
                }


@pytest.fixture
def mock_entrez():
    """Mock BioPython Entrez for PubMed API testing."""
    with pytest.mock.patch('src.pubmed_miner.services.paper_collection.Entrez') as mock_entrez:
        yield mock_entrez


@pytest.fixture(autouse=True)
def setup_test_environment():
    """Set up test environment variables and cleanup."""
    # Set test environment variables
    test_env = {
        'LOG_LEVEL': 'DEBUG',
        'PUBMED_EMAIL': 'test@example.com'
    }
    
    original_env = {}
    for key, value in test_env.items():
        original_env[key] = os.environ.get(key)
        os.environ[key] = value
    
    yield
    
    # Restore original environment
    for key, original_value in original_env.items():
        if original_value is not None:
            os.environ[key] = original_value
        else:
            os.environ.pop(key, None)


# Pytest markers for different test categories
def pytest_configure(config):
    """Configure pytest markers."""
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "e2e: mark test as an end-to-end test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "external: mark test as requiring external APIs"
    )


# Skip external tests by default
def pytest_collection_modifyitems(config, items):
    """Modify test collection to handle external test skipping."""
    if config.getoption("--skip-external"):
        skip_external = pytest.mark.skip(reason="--skip-external option given")
        for item in items:
            if "external" in item.keywords:
                item.add_marker(skip_external)


def pytest_addoption(parser):
    """Add custom command line options."""
    parser.addoption(
        "--skip-external",
        action="store_true",
        default=False,
        help="Skip tests that require external API access"
    )
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run slow tests"
    )


# Custom test result reporting
@pytest.fixture(scope="session", autouse=True)
def test_session_setup():
    """Set up test session."""
    print("\n" + "="*50)
    print("Essential Papers Recommender - Test Suite")
    print("="*50)
    
    yield
    
    print("\n" + "="*50)
    print("Test Session Complete")
    print("="*50)