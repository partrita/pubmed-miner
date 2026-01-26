"""
Error handling utilities and custom exceptions.

This module provides comprehensive error handling for the PubMed Miner system,
including custom exceptions, retry logic, and fallback strategies.

Requirements addressed:
- 1.4: API error handling for PubMed requests
- 4.4: Workflow error handling and logging
"""

import time
import logging
import traceback
from typing import Callable, Any, Optional, Dict, List
from functools import wraps
from datetime import datetime
from enum import Enum

logger = logging.getLogger(__name__)


class PubMedMinerError(Exception):
    """Base exception for PubMed Miner errors."""

    pass


class APIError(PubMedMinerError):
    """Exception raised for API-related errors."""

    pass


class DataError(PubMedMinerError):
    """Exception raised for data validation or processing errors."""

    pass


class ConfigurationError(PubMedMinerError):
    """Exception raised for configuration-related errors."""

    pass


class RateLimitError(APIError):
    """Exception raised when API rate limits are exceeded."""

    def __init__(self, message: str, retry_after: Optional[int] = None):
        super().__init__(message)
        self.retry_after = retry_after


class NetworkError(APIError):
    """Exception raised for network connectivity issues."""

    pass


class AuthenticationError(APIError):
    """Exception raised for authentication/authorization issues."""

    pass


class ValidationError(DataError):
    """Exception raised for data validation failures."""

    def __init__(self, message: str, field: Optional[str] = None, value: Any = None):
        super().__init__(message)
        self.field = field
        self.value = value


class CacheError(PubMedMinerError):
    """Exception raised for cache-related errors."""

    pass


class GitHubError(APIError):
    """Exception raised for GitHub API-related errors."""

    def __init__(self, message: str, status_code: Optional[int] = None):
        super().__init__(message)
        self.status_code = status_code


class PubMedError(APIError):
    """Exception raised for PubMed API-related errors."""

    pass


class ScoringError(DataError):
    """Exception raised for scoring calculation errors."""

    pass


class ErrorSeverity(Enum):
    """Enumeration of error severity levels."""

    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class ErrorHandler:
    """Comprehensive error handler with retry logic, fallback strategies, and monitoring."""

    def __init__(self):
        """Initialize error handler with tracking capabilities."""
        self.error_counts: Dict[str, int] = {}
        self.last_errors: Dict[str, datetime] = {}
        self.circuit_breakers: Dict[str, bool] = {}

    def handle_api_error(
        self,
        error: APIError,
        context: str = "",
        severity: ErrorSeverity = ErrorSeverity.MEDIUM,
    ) -> None:
        """Handle API-related errors with enhanced logging and tracking.

        Args:
            error: The API error that occurred
            context: Additional context about where the error occurred
            severity: Severity level of the error
        """
        error_key = f"{type(error).__name__}:{context}"
        self._track_error(error_key)

        log_level = self._get_log_level(severity)
        logger.log(
            log_level, f"API Error{' in ' + context if context else ''}: {error}"
        )

        # Log additional details for debugging
        if hasattr(error, "__cause__") and error.__cause__:
            logger.debug(f"Underlying cause: {error.__cause__}")

        # Handle specific API error types
        if isinstance(error, RateLimitError):
            self._handle_rate_limit_error(error)
        elif isinstance(error, AuthenticationError):
            self._handle_auth_error(error, context)
        elif isinstance(error, NetworkError):
            self._handle_network_error(error, context)
        elif isinstance(error, GitHubError):
            self._handle_github_error(error, context)
        elif isinstance(error, PubMedError):
            self._handle_pubmed_error(error, context)

    def handle_data_error(
        self,
        error: DataError,
        context: str = "",
        severity: ErrorSeverity = ErrorSeverity.MEDIUM,
    ) -> None:
        """Handle data validation or processing errors.

        Args:
            error: The data error that occurred
            context: Additional context about where the error occurred
            severity: Severity level of the error
        """
        error_key = f"{type(error).__name__}:{context}"
        self._track_error(error_key)

        log_level = self._get_log_level(severity)
        logger.log(
            log_level, f"Data Error{' in ' + context if context else ''}: {error}"
        )

        # Handle specific data error types
        if isinstance(error, ValidationError):
            self._handle_validation_error(error, context)
        elif isinstance(error, ScoringError):
            self._handle_scoring_error(error, context)

    def handle_configuration_error(
        self, error: ConfigurationError, context: str = ""
    ) -> None:
        """Handle configuration-related errors.

        Args:
            error: The configuration error that occurred
            context: Additional context about where the error occurred
        """
        error_key = f"ConfigurationError:{context}"
        self._track_error(error_key)

        logger.critical(
            f"Configuration Error{' in ' + context if context else ''}: {error}"
        )
        logger.critical("System cannot continue with invalid configuration")

    def handle_cache_error(self, error: CacheError, context: str = "") -> None:
        """Handle cache-related errors.

        Args:
            error: The cache error that occurred
            context: Additional context about where the error occurred
        """
        error_key = f"CacheError:{context}"
        self._track_error(error_key)

        logger.warning(f"Cache Error{' in ' + context if context else ''}: {error}")
        logger.info("Continuing without cache functionality")

    def _handle_rate_limit_error(self, error: RateLimitError) -> None:
        """Handle rate limiting errors with intelligent waiting."""
        wait_time = error.retry_after if error.retry_after else 60.0

        logger.warning(f"Rate limit exceeded: {error}")
        logger.info(f"Waiting {wait_time} seconds before retrying...")

        # Implement exponential backoff for repeated rate limits
        error_key = "RateLimitError"
        if error_key in self.error_counts and self.error_counts[error_key] > 1:
            wait_time *= min(2 ** (self.error_counts[error_key] - 1), 8)  # Cap at 8x
            logger.info(f"Applying exponential backoff: {wait_time} seconds")

        time.sleep(wait_time)

    def _handle_auth_error(self, error: AuthenticationError, context: str) -> None:
        """Handle authentication errors."""
        logger.error(f"Authentication failed in {context}: {error}")
        logger.error("Please check your API credentials and permissions")

        # Set circuit breaker for this context
        self.circuit_breakers[context] = True

    def _handle_network_error(self, error: NetworkError, context: str) -> None:
        """Handle network connectivity errors."""
        logger.warning(f"Network error in {context}: {error}")
        logger.info("This may be a temporary connectivity issue")

    def _handle_github_error(self, error: GitHubError, context: str) -> None:
        """Handle GitHub-specific errors."""
        if error.status_code:
            logger.error(f"GitHub API error {error.status_code} in {context}: {error}")

            if error.status_code == 403:
                logger.error("GitHub API rate limit or permission issue")
            elif error.status_code == 404:
                logger.error("GitHub resource not found - check repository settings")
            elif error.status_code >= 500:
                logger.warning("GitHub server error - may be temporary")
        else:
            logger.error(f"GitHub error in {context}: {error}")

    def _handle_pubmed_error(self, error: PubMedError, context: str) -> None:
        """Handle PubMed-specific errors."""
        logger.error(f"PubMed API error in {context}: {error}")
        logger.info("This may be due to PubMed server issues or query problems")

    def _handle_validation_error(self, error: ValidationError, context: str) -> None:
        """Handle data validation errors."""
        if error.field and error.value is not None:
            logger.error(
                f"Validation failed for field '{error.field}' with value '{error.value}': {error}"
            )
        else:
            logger.error(f"Data validation failed in {context}: {error}")

    def _handle_scoring_error(self, error: ScoringError, context: str) -> None:
        """Handle scoring calculation errors."""
        logger.error(f"Scoring calculation failed in {context}: {error}")
        logger.info("Paper may be excluded from results due to scoring issues")

    def _track_error(self, error_key: str) -> None:
        """Track error occurrences for monitoring and circuit breaking."""
        self.error_counts[error_key] = self.error_counts.get(error_key, 0) + 1
        self.last_errors[error_key] = datetime.now()

    def _get_log_level(self, severity: ErrorSeverity) -> int:
        """Get appropriate log level for error severity."""
        severity_map = {
            ErrorSeverity.LOW: logging.INFO,
            ErrorSeverity.MEDIUM: logging.WARNING,
            ErrorSeverity.HIGH: logging.ERROR,
            ErrorSeverity.CRITICAL: logging.CRITICAL,
        }
        return severity_map.get(severity, logging.WARNING)

    def is_circuit_open(self, context: str) -> bool:
        """Check if circuit breaker is open for a given context."""
        return self.circuit_breakers.get(context, False)

    def reset_circuit_breaker(self, context: str) -> None:
        """Reset circuit breaker for a given context."""
        self.circuit_breakers[context] = False

    def get_error_summary(self) -> Dict[str, Any]:
        """Get summary of tracked errors."""
        return {
            "error_counts": dict(self.error_counts),
            "last_errors": {k: v.isoformat() for k, v in self.last_errors.items()},
            "circuit_breakers": dict(self.circuit_breakers),
            "total_errors": sum(self.error_counts.values()),
        }

    def reset_error_tracking(self) -> None:
        """Reset all error tracking data."""
        self.error_counts.clear()
        self.last_errors.clear()
        self.circuit_breakers.clear()

    @staticmethod
    def retry_on_failure(
        max_attempts: int = 3,
        delay: float = 1.0,
        backoff_factor: float = 2.0,
        exceptions: tuple = (APIError, RateLimitError),
        jitter: bool = True,
    ):
        """Decorator to retry function calls on specific exceptions with enhanced strategies.

        Args:
            max_attempts: Maximum number of retry attempts
            delay: Initial delay between retries (seconds)
            backoff_factor: Factor to multiply delay by after each failure
            exceptions: Tuple of exception types to retry on
            jitter: Whether to add random jitter to delay times
        """

        def decorator(func: Callable) -> Callable:
            @wraps(func)
            def wrapper(*args, **kwargs) -> Any:
                import random

                last_exception = None
                current_delay = delay

                for attempt in range(max_attempts):
                    try:
                        return func(*args, **kwargs)
                    except exceptions as e:
                        last_exception = e

                        if attempt < max_attempts - 1:  # Don't sleep on last attempt
                            logger.warning(
                                f"Attempt {attempt + 1}/{max_attempts} failed for {func.__name__}: {e}"
                            )

                            # Handle rate limit errors specially
                            if isinstance(e, RateLimitError) and e.retry_after:
                                sleep_time = e.retry_after
                                logger.info(
                                    f"Rate limited, waiting {sleep_time} seconds as requested..."
                                )
                            else:
                                sleep_time = current_delay
                                if jitter:
                                    # Add ±25% jitter to prevent thundering herd
                                    jitter_range = sleep_time * 0.25
                                    sleep_time += random.uniform(
                                        -jitter_range, jitter_range
                                    )

                                logger.info(f"Retrying in {sleep_time:.2f} seconds...")

                            time.sleep(max(0, sleep_time))
                            current_delay *= backoff_factor
                        else:
                            logger.error(
                                f"All {max_attempts} attempts failed for {func.__name__}"
                            )

                # Re-raise the last exception if all attempts failed
                raise last_exception

            return wrapper

        return decorator

    @staticmethod
    def circuit_breaker(
        failure_threshold: int = 5,
        recovery_timeout: int = 60,
        expected_exception: type = APIError,
    ):
        """Decorator implementing circuit breaker pattern.

        Args:
            failure_threshold: Number of failures before opening circuit
            recovery_timeout: Seconds to wait before attempting recovery
            expected_exception: Exception type that triggers circuit breaker
        """

        def decorator(func: Callable) -> Callable:
            # Circuit breaker state
            state = {"failures": 0, "last_failure": None, "is_open": False}

            @wraps(func)
            def wrapper(*args, **kwargs) -> Any:
                now = datetime.now()

                # Check if circuit should be reset
                if (
                    state["is_open"]
                    and state["last_failure"]
                    and (now - state["last_failure"]).seconds >= recovery_timeout
                ):
                    logger.info(
                        f"Circuit breaker for {func.__name__} attempting recovery"
                    )
                    state["is_open"] = False
                    state["failures"] = 0

                # If circuit is open, fail fast
                if state["is_open"]:
                    raise expected_exception(
                        f"Circuit breaker open for {func.__name__}"
                    )

                try:
                    result = func(*args, **kwargs)
                    # Reset failure count on success
                    state["failures"] = 0
                    return result
                except expected_exception:
                    state["failures"] += 1
                    state["last_failure"] = now

                    if state["failures"] >= failure_threshold:
                        state["is_open"] = True
                        logger.error(
                            f"Circuit breaker opened for {func.__name__} after {failure_threshold} failures"
                        )

                    raise

            return wrapper

        return decorator

    @staticmethod
    def safe_execute(
        func: Callable,
        default_return: Any = None,
        log_errors: bool = True,
        context: str = "",
        allowed_exceptions: tuple = (),
    ) -> Any:
        """Safely execute a function with comprehensive error handling.

        Args:
            func: Function to execute
            default_return: Value to return if function fails
            log_errors: Whether to log errors
            context: Additional context for error logging
            allowed_exceptions: Exceptions that should be re-raised instead of handled

        Returns:
            Function result or default_return if function fails
        """
        try:
            return func()
        except allowed_exceptions:
            # Re-raise allowed exceptions
            raise
        except Exception as e:
            if log_errors:
                func_name = getattr(func, "__name__", str(func))
                logger.error(
                    f"Error in {func_name}{' (' + context + ')' if context else ''}: {e}"
                )
                logger.debug(f"Full traceback: {traceback.format_exc()}")
            return default_return

    @staticmethod
    def validate_and_handle(
        data: Any,
        validator: Callable[[Any], bool],
        error_message: str = "Data validation failed",
        field_name: Optional[str] = None,
    ) -> None:
        """Validate data and raise ValidationError if validation fails.

        Args:
            data: Data to validate
            validator: Function that returns True if data is valid
            error_message: Error message to use if validation fails
            field_name: Name of the field being validated

        Raises:
            ValidationError: If validation fails
        """
        try:
            if not validator(data):
                raise ValidationError(error_message, field=field_name, value=data)
        except ValidationError:
            raise
        except Exception as e:
            raise ValidationError(f"{error_message}: {e}", field=field_name, value=data)

    @staticmethod
    def batch_execute(
        items: List[Any],
        func: Callable[[Any], Any],
        max_failures: Optional[int] = None,
        continue_on_error: bool = True,
        context: str = "",
    ) -> Dict[str, Any]:
        """Execute a function on a batch of items with error tracking.

        Args:
            items: List of items to process
            func: Function to execute on each item
            max_failures: Maximum failures before stopping (None = no limit)
            continue_on_error: Whether to continue processing after errors
            context: Context for error logging

        Returns:
            Dictionary with results, errors, and statistics
        """
        results = []
        errors = []
        processed = 0

        for i, item in enumerate(items):
            try:
                result = func(item)
                results.append({"index": i, "item": item, "result": result})
                processed += 1
            except Exception as e:
                error_info = {
                    "index": i,
                    "item": item,
                    "error": str(e),
                    "error_type": type(e).__name__,
                }
                errors.append(error_info)

                logger.warning(f"Error processing item {i} in {context}: {e}")

                # Check if we should stop due to too many failures
                if max_failures and len(errors) >= max_failures:
                    logger.error(
                        f"Stopping batch processing after {len(errors)} failures"
                    )
                    break

                if not continue_on_error:
                    logger.error(f"Stopping batch processing due to error: {e}")
                    break

        return {
            "results": results,
            "errors": errors,
            "total_items": len(items),
            "processed": processed,
            "failed": len(errors),
            "success_rate": processed / len(items) if items else 0,
        }


# Global error handler instance
_global_error_handler = ErrorHandler()


def get_error_handler() -> ErrorHandler:
    """Get the global error handler instance."""
    return _global_error_handler


# Convenience decorators for common retry patterns
def retry_api_calls(max_attempts: int = 3, delay: float = 1.0, jitter: bool = True):
    """Decorator for retrying API calls with intelligent backoff."""
    return ErrorHandler.retry_on_failure(
        max_attempts=max_attempts,
        delay=delay,
        exceptions=(APIError, RateLimitError, NetworkError),
        jitter=jitter,
    )


def retry_data_operations(max_attempts: int = 2, delay: float = 0.5):
    """Decorator for retrying data operations."""
    return ErrorHandler.retry_on_failure(
        max_attempts=max_attempts, delay=delay, exceptions=(DataError, ValidationError)
    )


def retry_github_operations(max_attempts: int = 3, delay: float = 2.0):
    """Decorator for retrying GitHub API operations."""
    return ErrorHandler.retry_on_failure(
        max_attempts=max_attempts,
        delay=delay,
        exceptions=(GitHubError, RateLimitError, NetworkError),
    )


def retry_pubmed_operations(max_attempts: int = 4, delay: float = 1.0):
    """Decorator for retrying PubMed API operations."""
    return ErrorHandler.retry_on_failure(
        max_attempts=max_attempts,
        delay=delay,
        exceptions=(PubMedError, RateLimitError, NetworkError),
    )


def handle_exceptions(
    default_return: Any = None, log_errors: bool = True, context: str = ""
):
    """Decorator to handle all exceptions and return a default value."""

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            return ErrorHandler.safe_execute(
                lambda: func(*args, **kwargs),
                default_return=default_return,
                log_errors=log_errors,
                context=context,
            )

        return wrapper

    return decorator


def circuit_breaker_api(failure_threshold: int = 5, recovery_timeout: int = 60):
    """Circuit breaker decorator for API operations."""
    return ErrorHandler.circuit_breaker(
        failure_threshold=failure_threshold,
        recovery_timeout=recovery_timeout,
        expected_exception=APIError,
    )


def validate_required_fields(required_fields: List[str]):
    """Decorator to validate required fields in data objects."""

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            # Check if first argument has required fields
            if args and hasattr(args[0], "__dict__"):
                obj = args[0]
                for field in required_fields:
                    if not hasattr(obj, field) or getattr(obj, field) is None:
                        raise ValidationError(
                            f"Required field '{field}' is missing or None", field=field
                        )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def log_execution_time(logger_instance: Optional[logging.Logger] = None):
    """Decorator to log function execution time."""

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            start_time = time.time()
            log = logger_instance or logger

            try:
                result = func(*args, **kwargs)
                execution_time = time.time() - start_time
                log.debug(f"{func.__name__} completed in {execution_time:.2f} seconds")
                return result
            except Exception as e:
                execution_time = time.time() - start_time
                log.warning(
                    f"{func.__name__} failed after {execution_time:.2f} seconds: {e}"
                )
                raise

        return wrapper

    return decorator
