import pytest
from src.pubmed_miner.utils.cache import CacheManager
from src.pubmed_miner.utils.error_handler import CacheError

def test_export_cache_data_sql_injection():
    """Verify that export_cache_data prevents SQL injection in table names."""
    manager = CacheManager(enable_memory_cache=False)

    with pytest.raises(CacheError, match="Invalid table name: citations; DROP TABLE citations;"):
        manager.export_cache_data("test.json", table_name="citations; DROP TABLE citations;")

    with pytest.raises(CacheError, match="Invalid table name: sqlite_master"):
        manager.export_cache_data("test.json", table_name="sqlite_master")
