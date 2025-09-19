"""
Cache management utilities for storing and retrieving data.

This module provides comprehensive caching functionality for the PubMed Miner system,
including SQLite-based storage, automatic cleanup, and performance monitoring.

Requirements addressed:
- 2.5: Citation and impact factor caching
- 7.5: Journal data caching and management
"""

import json
import sqlite3
import logging
import threading
import hashlib
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any, Callable
from contextlib import contextmanager

from ..models.cache import CitationCache, ImpactFactorCache, PaperMetadataCache
from .error_handler import CacheError, handle_exceptions

logger = logging.getLogger(__name__)


class CacheManager:
    """Manages local caching of paper data, citations, and impact factors with enhanced features."""

    def __init__(self, cache_dir: str = "cache", enable_memory_cache: bool = True):
        """Initialize cache manager.

        Args:
            cache_dir: Directory to store cache files
            enable_memory_cache: Whether to enable in-memory caching for performance
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        self.db_path = self.cache_dir / "pubmed_cache.db"
        self.enable_memory_cache = enable_memory_cache

        # Thread-safe in-memory cache
        self._memory_cache: Dict[str, Any] = {}
        self._cache_lock = threading.RLock()
        self._cache_stats = {"hits": 0, "misses": 0, "memory_hits": 0, "db_hits": 0}

        self._init_database()

        logger.info(f"Initialized CacheManager with database: {self.db_path}")
        logger.info(f"Memory cache enabled: {enable_memory_cache}")

    def _init_database(self) -> None:
        """Initialize SQLite database with required tables."""
        with self._get_connection() as conn:
            # Citations table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS citations (
                    pmid TEXT PRIMARY KEY,
                    citation_count INTEGER NOT NULL,
                    last_updated TIMESTAMP NOT NULL,
                    source TEXT NOT NULL
                )
            """)

            # Impact factors table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS impact_factors (
                    journal_name TEXT PRIMARY KEY,
                    impact_factor REAL NOT NULL,
                    year INTEGER NOT NULL,
                    last_updated TIMESTAMP NOT NULL,
                    source TEXT NOT NULL
                )
            """)

            # Paper metadata table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS paper_metadata (
                    pmid TEXT PRIMARY KEY,
                    title TEXT NOT NULL,
                    authors_json TEXT NOT NULL,
                    journal TEXT NOT NULL,
                    publication_date TIMESTAMP NOT NULL,
                    abstract TEXT,
                    doi TEXT,
                    last_updated TIMESTAMP NOT NULL
                )
            """)

            # Create indexes for better performance
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_citations_updated ON citations(last_updated)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_impact_factors_updated ON impact_factors(last_updated)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_metadata_updated ON paper_metadata(last_updated)"
            )

            conn.commit()

    @contextmanager
    def _get_connection(self):
        """Get database connection with proper error handling."""
        conn = None
        try:
            conn = sqlite3.connect(
                self.db_path,
                detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
            )
            conn.row_factory = sqlite3.Row
            yield conn
        except Exception as e:
            if conn:
                conn.rollback()
            logger.error(f"Database error: {e}")
            raise
        finally:
            if conn:
                conn.close()

    # Citation cache methods
    def get_citation(self, pmid: str) -> Optional[CitationCache]:
        """Get citation data from cache.

        Args:
            pmid: PubMed ID

        Returns:
            CitationCache object or None if not found
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute("SELECT * FROM citations WHERE pmid = ?", (pmid,))
                row = cursor.fetchone()

                if row:
                    return CitationCache(
                        pmid=row["pmid"],
                        citation_count=row["citation_count"],
                        last_updated=row["last_updated"],
                        source=row["source"],
                    )

        except Exception as e:
            logger.error(f"Error getting citation cache for {pmid}: {e}")

        return None

    def save_citation(self, citation: CitationCache) -> None:
        """Save citation data to cache.

        Args:
            citation: CitationCache object to save
        """
        try:
            with self._get_connection() as conn:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO citations 
                    (pmid, citation_count, last_updated, source)
                    VALUES (?, ?, ?, ?)
                """,
                    (
                        citation.pmid,
                        citation.citation_count,
                        citation.last_updated,
                        citation.source,
                    ),
                )
                conn.commit()

        except Exception as e:
            logger.error(f"Error saving citation cache for {citation.pmid}: {e}")

    def batch_save_citations(self, citations: List[CitationCache]) -> None:
        """Save multiple citations to cache.

        Args:
            citations: List of CitationCache objects
        """
        if not citations:
            return

        try:
            with self._get_connection() as conn:
                data = [
                    (c.pmid, c.citation_count, c.last_updated, c.source)
                    for c in citations
                ]

                conn.executemany(
                    """
                    INSERT OR REPLACE INTO citations 
                    (pmid, citation_count, last_updated, source)
                    VALUES (?, ?, ?, ?)
                """,
                    data,
                )
                conn.commit()

                logger.info(f"Saved {len(citations)} citations to cache")

        except Exception as e:
            logger.error(f"Error batch saving citations: {e}")

    def clear_expired_citations(self, max_age_days: int = 7) -> int:
        """Clear expired citation cache entries.

        Args:
            max_age_days: Maximum age in days

        Returns:
            Number of entries cleared
        """
        try:
            cutoff_date = datetime.now() - timedelta(days=max_age_days)

            with self._get_connection() as conn:
                cursor = conn.execute(
                    "DELETE FROM citations WHERE last_updated < ?", (cutoff_date,)
                )
                cleared_count = cursor.rowcount
                conn.commit()

                return cleared_count

        except Exception as e:
            logger.error(f"Error clearing expired citations: {e}")
            return 0

    # Impact factor cache methods
    def get_impact_factor(self, journal_name: str) -> Optional[ImpactFactorCache]:
        """Get impact factor from cache.

        Args:
            journal_name: Journal name

        Returns:
            ImpactFactorCache object or None if not found
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute(
                    "SELECT * FROM impact_factors WHERE journal_name = ?",
                    (journal_name,),
                )
                row = cursor.fetchone()

                if row:
                    return ImpactFactorCache(
                        journal_name=row["journal_name"],
                        impact_factor=row["impact_factor"],
                        year=row["year"],
                        last_updated=row["last_updated"],
                        source=row["source"],
                    )

        except Exception as e:
            logger.error(f"Error getting impact factor cache for {journal_name}: {e}")

        return None

    def save_impact_factor(self, impact_factor: ImpactFactorCache) -> None:
        """Save impact factor to cache.

        Args:
            impact_factor: ImpactFactorCache object to save
        """
        try:
            with self._get_connection() as conn:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO impact_factors 
                    (journal_name, impact_factor, year, last_updated, source)
                    VALUES (?, ?, ?, ?, ?)
                """,
                    (
                        impact_factor.journal_name,
                        impact_factor.impact_factor,
                        impact_factor.year,
                        impact_factor.last_updated,
                        impact_factor.source,
                    ),
                )
                conn.commit()

        except Exception as e:
            logger.error(
                f"Error saving impact factor cache for {impact_factor.journal_name}: {e}"
            )

    def search_similar_journals(
        self, journal_name: str, limit: int = 5
    ) -> List[ImpactFactorCache]:
        """Search for journals with similar names.

        Args:
            journal_name: Journal name to search for
            limit: Maximum number of results

        Returns:
            List of similar journal impact factors
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute(
                    """
                    SELECT * FROM impact_factors 
                    WHERE journal_name LIKE ? 
                    ORDER BY impact_factor DESC
                    LIMIT ?
                """,
                    (f"%{journal_name}%", limit),
                )

                results = []
                for row in cursor.fetchall():
                    results.append(
                        ImpactFactorCache(
                            journal_name=row["journal_name"],
                            impact_factor=row["impact_factor"],
                            year=row["year"],
                            last_updated=row["last_updated"],
                            source=row["source"],
                        )
                    )

                return results

        except Exception as e:
            logger.error(f"Error searching similar journals for {journal_name}: {e}")
            return []

    # Paper metadata cache methods
    def get_paper_metadata(self, pmid: str) -> Optional[PaperMetadataCache]:
        """Get paper metadata from cache.

        Args:
            pmid: PubMed ID

        Returns:
            PaperMetadataCache object or None if not found
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute(
                    "SELECT * FROM paper_metadata WHERE pmid = ?", (pmid,)
                )
                row = cursor.fetchone()

                if row:
                    return PaperMetadataCache(
                        pmid=row["pmid"],
                        title=row["title"],
                        authors_json=row["authors_json"],
                        journal=row["journal"],
                        publication_date=row["publication_date"],
                        abstract=row["abstract"],
                        doi=row["doi"],
                        last_updated=row["last_updated"],
                    )

        except Exception as e:
            logger.error(f"Error getting paper metadata cache for {pmid}: {e}")

        return None

    def save_paper_metadata(self, metadata: PaperMetadataCache) -> None:
        """Save paper metadata to cache.

        Args:
            metadata: PaperMetadataCache object to save
        """
        try:
            with self._get_connection() as conn:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO paper_metadata 
                    (pmid, title, authors_json, journal, publication_date, abstract, doi, last_updated)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                    (
                        metadata.pmid,
                        metadata.title,
                        metadata.authors_json,
                        metadata.journal,
                        metadata.publication_date,
                        metadata.abstract,
                        metadata.doi,
                        metadata.last_updated,
                    ),
                )
                conn.commit()

        except Exception as e:
            logger.error(f"Error saving paper metadata cache for {metadata.pmid}: {e}")

    # Statistics and maintenance methods
    def get_citation_cache_stats(self) -> Dict[str, int]:
        """Get citation cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        try:
            with self._get_connection() as conn:
                # Total citations
                cursor = conn.execute("SELECT COUNT(*) as total FROM citations")
                total = cursor.fetchone()["total"]

                # Expired citations (older than 7 days)
                cutoff_date = datetime.now() - timedelta(days=7)
                cursor = conn.execute(
                    "SELECT COUNT(*) as expired FROM citations WHERE last_updated < ?",
                    (cutoff_date,),
                )
                expired = cursor.fetchone()["expired"]

                return {"total": total, "expired": expired, "valid": total - expired}

        except Exception as e:
            logger.error(f"Error getting citation cache stats: {e}")
            return {"total": 0, "expired": 0, "valid": 0}

    def get_impact_factor_cache_stats(self) -> Dict[str, int]:
        """Get impact factor cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        try:
            with self._get_connection() as conn:
                # Total impact factors
                cursor = conn.execute("SELECT COUNT(*) as total FROM impact_factors")
                total = cursor.fetchone()["total"]

                # Expired impact factors (older than 365 days)
                cutoff_date = datetime.now() - timedelta(days=365)
                cursor = conn.execute(
                    "SELECT COUNT(*) as expired FROM impact_factors WHERE last_updated < ?",
                    (cutoff_date,),
                )
                expired = cursor.fetchone()["expired"]

                return {"total": total, "expired": expired, "valid": total - expired}

        except Exception as e:
            logger.error(f"Error getting impact factor cache stats: {e}")
            return {"total": 0, "expired": 0, "valid": 0}

    def cleanup_cache(
        self, citation_max_age: int = 7, impact_factor_max_age: int = 365
    ) -> Dict[str, int]:
        """Clean up expired cache entries.

        Args:
            citation_max_age: Maximum age for citations in days
            impact_factor_max_age: Maximum age for impact factors in days

        Returns:
            Dictionary with cleanup statistics
        """
        stats = {
            "citations_cleared": 0,
            "impact_factors_cleared": 0,
            "metadata_cleared": 0,
        }

        try:
            # Clear expired citations
            stats["citations_cleared"] = self.clear_expired_citations(citation_max_age)

            # Clear expired impact factors
            cutoff_date = datetime.now() - timedelta(days=impact_factor_max_age)
            with self._get_connection() as conn:
                cursor = conn.execute(
                    "DELETE FROM impact_factors WHERE last_updated < ?", (cutoff_date,)
                )
                stats["impact_factors_cleared"] = cursor.rowcount

                # Clear expired metadata (30 days)
                metadata_cutoff = datetime.now() - timedelta(days=30)
                cursor = conn.execute(
                    "DELETE FROM paper_metadata WHERE last_updated < ?",
                    (metadata_cutoff,),
                )
                stats["metadata_cleared"] = cursor.rowcount

                conn.commit()

            logger.info(f"Cache cleanup completed: {stats}")

        except Exception as e:
            logger.error(f"Error during cache cleanup: {e}")

        return stats

    def vacuum_database(self) -> None:
        """Vacuum the database to reclaim space."""
        try:
            with self._get_connection() as conn:
                conn.execute("VACUUM")
                conn.commit()
            logger.info("Database vacuum completed")
        except Exception as e:
            logger.error(f"Error vacuuming database: {e}")

    def get_database_size(self) -> int:
        """Get database file size in bytes.

        Returns:
            Database file size in bytes
        """
        try:
            return self.db_path.stat().st_size
        except Exception as e:
            logger.error(f"Error getting database size: {e}")
            return 0

    # Enhanced utility methods
    def _generate_cache_key(self, prefix: str, identifier: str) -> str:
        """Generate a cache key for memory cache.

        Args:
            prefix: Cache type prefix
            identifier: Unique identifier

        Returns:
            Cache key string
        """
        return f"{prefix}:{identifier}"

    def _get_from_memory_cache(self, key: str) -> Optional[Any]:
        """Get item from memory cache.

        Args:
            key: Cache key

        Returns:
            Cached item or None
        """
        if not self.enable_memory_cache:
            return None

        with self._cache_lock:
            if key in self._memory_cache:
                self._cache_stats["memory_hits"] += 1
                return self._memory_cache[key]
            return None

    def _set_memory_cache(self, key: str, value: Any, max_size: int = 1000) -> None:
        """Set item in memory cache with size limit.

        Args:
            key: Cache key
            value: Value to cache
            max_size: Maximum number of items in memory cache
        """
        if not self.enable_memory_cache:
            return

        with self._cache_lock:
            # Simple LRU: remove oldest items if cache is full
            if len(self._memory_cache) >= max_size:
                # Remove first 10% of items (simple cleanup)
                items_to_remove = list(self._memory_cache.keys())[: max_size // 10]
                for item_key in items_to_remove:
                    del self._memory_cache[item_key]

            self._memory_cache[key] = value

    def clear_memory_cache(self) -> None:
        """Clear all items from memory cache."""
        with self._cache_lock:
            self._memory_cache.clear()
            logger.info("Memory cache cleared")

    def get_cache_statistics(self) -> Dict[str, Any]:
        """Get comprehensive cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        db_stats = {
            "citations": self.get_citation_cache_stats(),
            "impact_factors": self.get_impact_factor_cache_stats(),
            "database_size_bytes": self.get_database_size(),
            "database_size_mb": round(self.get_database_size() / (1024 * 1024), 2),
        }

        memory_stats = {
            "memory_cache_enabled": self.enable_memory_cache,
            "memory_cache_size": len(self._memory_cache),
            "cache_hits": self._cache_stats["hits"],
            "cache_misses": self._cache_stats["misses"],
            "memory_hits": self._cache_stats["memory_hits"],
            "db_hits": self._cache_stats["db_hits"],
            "hit_rate": (
                self._cache_stats["hits"]
                / (self._cache_stats["hits"] + self._cache_stats["misses"])
                if (self._cache_stats["hits"] + self._cache_stats["misses"]) > 0
                else 0
            ),
        }

        return {
            "database": db_stats,
            "memory": memory_stats,
            "timestamp": datetime.now().isoformat(),
        }

    def export_cache_data(
        self, output_file: str, table_name: Optional[str] = None
    ) -> None:
        """Export cache data to JSON file.

        Args:
            output_file: Output file path
            table_name: Specific table to export (None for all tables)
        """
        try:
            export_data = {}

            with self._get_connection() as conn:
                tables = (
                    [table_name]
                    if table_name
                    else ["citations", "impact_factors", "paper_metadata"]
                )

                for table in tables:
                    cursor = conn.execute(f"SELECT * FROM {table}")
                    rows = cursor.fetchall()
                    export_data[table] = [dict(row) for row in rows]

            # Convert datetime objects to strings for JSON serialization
            def datetime_converter(obj):
                if isinstance(obj, datetime):
                    return obj.isoformat()
                raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(
                    export_data,
                    f,
                    indent=2,
                    default=datetime_converter,
                    ensure_ascii=False,
                )

            logger.info(f"Cache data exported to {output_file}")

        except Exception as e:
            logger.error(f"Error exporting cache data: {e}")
            raise CacheError(f"Failed to export cache data: {e}")

    def import_cache_data(
        self, input_file: str, overwrite: bool = False
    ) -> Dict[str, int]:
        """Import cache data from JSON file.

        Args:
            input_file: Input file path
            overwrite: Whether to overwrite existing data

        Returns:
            Dictionary with import statistics
        """
        try:
            with open(input_file, "r", encoding="utf-8") as f:
                import_data = json.load(f)

            stats = {"citations": 0, "impact_factors": 0, "paper_metadata": 0}

            with self._get_connection() as conn:
                for table_name, rows in import_data.items():
                    if table_name not in stats:
                        continue

                    for row in rows:
                        # Convert ISO datetime strings back to datetime objects
                        for key, value in row.items():
                            if key.endswith("_updated") or key == "publication_date":
                                if isinstance(value, str):
                                    try:
                                        row[key] = datetime.fromisoformat(
                                            value.replace("Z", "+00:00")
                                        )
                                    except ValueError:
                                        pass

                        # Insert or replace based on overwrite setting
                        operation = (
                            "INSERT OR REPLACE" if overwrite else "INSERT OR IGNORE"
                        )

                        if table_name == "citations":
                            conn.execute(
                                f"""
                                {operation} INTO citations 
                                (pmid, citation_count, last_updated, source)
                                VALUES (?, ?, ?, ?)
                            """,
                                (
                                    row["pmid"],
                                    row["citation_count"],
                                    row["last_updated"],
                                    row["source"],
                                ),
                            )

                        elif table_name == "impact_factors":
                            conn.execute(
                                f"""
                                {operation} INTO impact_factors 
                                (journal_name, impact_factor, year, last_updated, source)
                                VALUES (?, ?, ?, ?, ?)
                            """,
                                (
                                    row["journal_name"],
                                    row["impact_factor"],
                                    row["year"],
                                    row["last_updated"],
                                    row["source"],
                                ),
                            )

                        elif table_name == "paper_metadata":
                            conn.execute(
                                f"""
                                {operation} INTO paper_metadata 
                                (pmid, title, authors_json, journal, publication_date, abstract, doi, last_updated)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                            """,
                                (
                                    row["pmid"],
                                    row["title"],
                                    row["authors_json"],
                                    row["journal"],
                                    row["publication_date"],
                                    row["abstract"],
                                    row["doi"],
                                    row["last_updated"],
                                ),
                            )

                        stats[table_name] += 1

                conn.commit()

            logger.info(f"Cache data imported: {stats}")
            return stats

        except Exception as e:
            logger.error(f"Error importing cache data: {e}")
            raise CacheError(f"Failed to import cache data: {e}")

    def optimize_database(self) -> None:
        """Optimize database performance by analyzing and rebuilding indexes."""
        try:
            with self._get_connection() as conn:
                # Analyze tables for query optimization
                conn.execute("ANALYZE")

                # Rebuild indexes
                conn.execute("REINDEX")

                conn.commit()

            logger.info("Database optimization completed")

        except Exception as e:
            logger.error(f"Error optimizing database: {e}")
            raise CacheError(f"Failed to optimize database: {e}")

    def backup_database(self, backup_path: str) -> None:
        """Create a backup of the cache database.

        Args:
            backup_path: Path for the backup file
        """
        try:
            import shutil

            shutil.copy2(self.db_path, backup_path)
            logger.info(f"Database backed up to {backup_path}")

        except Exception as e:
            logger.error(f"Error backing up database: {e}")
            raise CacheError(f"Failed to backup database: {e}")


# Utility functions for cache management
def create_cache_key(prefix: str, *args: Any) -> str:
    """Create a consistent cache key from arguments.

    Args:
        prefix: Key prefix
        *args: Arguments to include in key

    Returns:
        Cache key string
    """
    key_parts = [str(prefix)]
    for arg in args:
        if isinstance(arg, (dict, list)):
            # Create hash for complex objects
            key_parts.append(
                hashlib.md5(json.dumps(arg, sort_keys=True).encode()).hexdigest()[:8]
            )
        else:
            key_parts.append(str(arg))
    return ":".join(key_parts)


def cache_result(cache_manager: CacheManager, key_prefix: str, ttl_hours: int = 24):
    """Decorator to cache function results.

    Args:
        cache_manager: CacheManager instance
        key_prefix: Prefix for cache keys
        ttl_hours: Time to live in hours
    """

    def decorator(func: Callable) -> Callable:
        def wrapper(*args, **kwargs) -> Any:
            # Create cache key from function arguments
            cache_key = create_cache_key(key_prefix, func.__name__, args, kwargs)

            # Try to get from memory cache first
            cached_result = cache_manager._get_from_memory_cache(cache_key)
            if cached_result is not None:
                return cached_result

            # Execute function and cache result
            result = func(*args, **kwargs)

            # Cache in memory for quick access
            cache_manager._set_memory_cache(cache_key, result)

            return result

        return wrapper

    return decorator


@handle_exceptions(default_return=None, log_errors=True, context="cache_operation")
def safe_cache_operation(operation: Callable[[], Any]) -> Any:
    """Safely execute a cache operation with error handling.

    Args:
        operation: Cache operation to execute

    Returns:
        Operation result or None if failed
    """
    return operation()
