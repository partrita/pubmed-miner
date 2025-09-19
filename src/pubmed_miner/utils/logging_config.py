"""
Logging configuration and utilities for the PubMed Miner system.

This module provides centralized logging configuration with support for
different log levels, file rotation, and structured logging.
"""

import os
import sys
import logging
import logging.handlers
import functools
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds colors to console output."""

    # ANSI color codes
    COLORS = {
        "DEBUG": "\033[36m",  # Cyan
        "INFO": "\033[32m",  # Green
        "WARNING": "\033[33m",  # Yellow
        "ERROR": "\033[31m",  # Red
        "CRITICAL": "\033[35m",  # Magenta
        "RESET": "\033[0m",  # Reset
    }

    def format(self, record):
        # Add color to levelname
        if record.levelname in self.COLORS:
            record.levelname = f"{self.COLORS[record.levelname]}{record.levelname}{self.COLORS['RESET']}"
        return super().format(record)


class StructuredFormatter(logging.Formatter):
    """Formatter that outputs structured JSON logs."""

    def format(self, record):
        log_entry = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }

        # Add exception info if present
        if record.exc_info:
            log_entry["exception"] = self.formatException(record.exc_info)

        # Add extra fields if present
        for key, value in record.__dict__.items():
            if key not in [
                "name",
                "msg",
                "args",
                "levelname",
                "levelno",
                "pathname",
                "filename",
                "module",
                "lineno",
                "funcName",
                "created",
                "msecs",
                "relativeCreated",
                "thread",
                "threadName",
                "processName",
                "process",
                "getMessage",
                "exc_info",
                "exc_text",
                "stack_info",
            ]:
                log_entry[key] = value

        import json

        return json.dumps(log_entry, ensure_ascii=False)


def setup_logging(
    log_level: str = "INFO",
    log_dir: str = "logs",
    log_file: str = "pubmed_miner.log",
    console_output: bool = True,
    file_output: bool = True,
    structured_logging: bool = False,
    max_file_size: int = 10 * 1024 * 1024,  # 10MB
    backup_count: int = 5,
) -> logging.Logger:
    """Set up comprehensive logging configuration.

    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_dir: Directory for log files
        log_file: Name of the log file
        console_output: Whether to output to console
        file_output: Whether to output to file
        structured_logging: Whether to use structured JSON logging
        max_file_size: Maximum size of log file before rotation
        backup_count: Number of backup files to keep

    Returns:
        Configured logger instance
    """
    # Create logs directory
    log_path = Path(log_dir)
    log_path.mkdir(exist_ok=True)

    # Get numeric log level
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)

    # Create root logger
    logger = logging.getLogger()
    logger.setLevel(numeric_level)

    # Clear existing handlers
    logger.handlers.clear()

    # Console handler
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(numeric_level)

        if structured_logging:
            console_formatter = StructuredFormatter()
        else:
            console_formatter = ColoredFormatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

    # File handler with rotation
    if file_output:
        file_path = log_path / log_file
        file_handler = logging.handlers.RotatingFileHandler(
            file_path,
            maxBytes=max_file_size,
            backupCount=backup_count,
            encoding="utf-8",
        )
        file_handler.setLevel(numeric_level)

        if structured_logging:
            file_formatter = StructuredFormatter()
        else:
            file_formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s"
            )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    # Error file handler (always enabled for ERROR and CRITICAL)
    if file_output:
        error_file_path = log_path / f"error_{log_file}"
        error_handler = logging.handlers.RotatingFileHandler(
            error_file_path,
            maxBytes=max_file_size,
            backupCount=backup_count,
            encoding="utf-8",
        )
        error_handler.setLevel(logging.ERROR)

        error_formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s\n"
            "Exception: %(exc_info)s\n"
        )
        error_handler.setFormatter(error_formatter)
        logger.addHandler(error_handler)

    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance with the specified name.

    Args:
        name: Logger name (typically __name__)

    Returns:
        Logger instance
    """
    return logging.getLogger(name)


def log_function_call(logger: Optional[logging.Logger] = None):
    """Decorator to log function calls with parameters and execution time.

    Args:
        logger: Logger instance to use (defaults to function's module logger)
    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            import time

            # Get logger
            log = logger or logging.getLogger(func.__module__)

            # Log function entry
            log.debug(f"Entering {func.__name__} with args={args}, kwargs={kwargs}")

            start_time = time.time()
            try:
                result = func(*args, **kwargs)
                execution_time = time.time() - start_time
                log.debug(f"Exiting {func.__name__} after {execution_time:.3f}s")
                return result
            except Exception as e:
                execution_time = time.time() - start_time
                log.error(
                    f"Exception in {func.__name__} after {execution_time:.3f}s: {e}"
                )
                raise

        return functools.wraps(func)(wrapper)

    return decorator


def log_performance_metrics(
    operation_name: str,
    metrics: Dict[str, Any],
    logger: Optional[logging.Logger] = None,
) -> None:
    """Log performance metrics for monitoring.

    Args:
        operation_name: Name of the operation
        metrics: Dictionary of metrics to log
        logger: Logger instance to use
    """
    log = logger or logging.getLogger(__name__)

    # Format metrics for logging
    metrics_str = ", ".join([f"{k}={v}" for k, v in metrics.items()])
    log.info(f"METRICS [{operation_name}]: {metrics_str}")


def configure_third_party_loggers(level: str = "WARNING") -> None:
    """Configure third-party library loggers to reduce noise.

    Args:
        level: Log level for third-party loggers
    """
    third_party_loggers = ["urllib3", "requests", "biopython", "github", "sqlite3"]

    numeric_level = getattr(logging, level.upper(), logging.WARNING)

    for logger_name in third_party_loggers:
        logging.getLogger(logger_name).setLevel(numeric_level)


def setup_github_actions_logging() -> None:
    """Configure logging for GitHub Actions environment."""
    if os.getenv("GITHUB_ACTIONS"):
        # GitHub Actions specific configuration
        setup_logging(
            log_level=os.getenv("LOG_LEVEL", "INFO"),
            console_output=True,
            file_output=True,
            structured_logging=False,  # GitHub Actions prefers plain text
        )

        # Add GitHub Actions specific formatters
        logger = logging.getLogger()
        for handler in logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                # Use GitHub Actions workflow commands
                formatter = GitHubActionsFormatter()
                handler.setFormatter(formatter)


class GitHubActionsFormatter(logging.Formatter):
    """Formatter that outputs GitHub Actions workflow commands."""

    def format(self, record):
        message = record.getMessage()

        if record.levelno >= logging.ERROR:
            return f"::error::{message}"
        elif record.levelno >= logging.WARNING:
            return f"::warning::{message}"
        elif record.levelno >= logging.INFO:
            return f"::notice::{message}"
        else:
            return f"::debug::{message}"


class LogContext:
    """Context manager for adding context to log messages."""

    def __init__(self, logger: logging.Logger, **context):
        """Initialize log context.

        Args:
            logger: Logger instance
            **context: Context key-value pairs
        """
        self.logger = logger
        self.context = context
        self.old_factory = logging.getLogRecordFactory()

    def __enter__(self):
        def record_factory(*args, **kwargs):
            record = self.old_factory(*args, **kwargs)
            for key, value in self.context.items():
                setattr(record, key, value)
            return record

        logging.setLogRecordFactory(record_factory)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logging.setLogRecordFactory(self.old_factory)


# Initialize logging on module import
def initialize_default_logging():
    """Initialize default logging configuration."""
    log_level = os.getenv("LOG_LEVEL", "INFO")

    # Check if we're in GitHub Actions
    if os.getenv("GITHUB_ACTIONS"):
        setup_github_actions_logging()
    else:
        setup_logging(log_level=log_level)

    # Configure third-party loggers
    configure_third_party_loggers()


# Auto-initialize if not already configured
if not logging.getLogger().handlers:
    initialize_default_logging()
