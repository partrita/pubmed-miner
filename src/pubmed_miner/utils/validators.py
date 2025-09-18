"""
Configuration validation utilities and data validators.

This module provides comprehensive validation functions for configuration files,
data objects, and API responses used throughout the PubMed Miner system.

Requirements addressed:
- Data validation and normalization for all system components
- Configuration file validation and error reporting
"""
import re
import os
import json
from typing import Dict, Any, List, Optional, Union, Callable
from pathlib import Path
from datetime import datetime
from urllib.parse import urlparse

from .error_handler import ValidationError, ConfigurationError


class ConfigValidator:
    """Validates configuration files and settings."""
    
    @staticmethod
    def validate_topic_config(topic_data: Dict[str, Any]) -> List[str]:
        """Validate a single topic configuration.
        
        Args:
            topic_data: Dictionary containing topic configuration
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        
        # Required fields
        required_fields = ['name', 'query']
        for field in required_fields:
            if field not in topic_data:
                errors.append(f"Missing required field: {field}")
            elif not topic_data[field]:
                errors.append(f"Field '{field}' cannot be empty")
                
        # Validate name format (alphanumeric, hyphens, underscores only)
        if 'name' in topic_data:
            name = topic_data['name']
            if not re.match(r'^[a-zA-Z0-9_-]+$', name):
                errors.append(f"Topic name '{name}' contains invalid characters. Use only letters, numbers, hyphens, and underscores.")
                
        # Validate numeric fields
        numeric_fields = {
            'max_papers': (1, 10000),
            'essential_count': (1, 100)
        }
        
        for field, (min_val, max_val) in numeric_fields.items():
            if field in topic_data:
                value = topic_data[field]
                if not isinstance(value, int):
                    errors.append(f"Field '{field}' must be an integer")
                elif value < min_val or value > max_val:
                    errors.append(f"Field '{field}' must be between {min_val} and {max_val}")
                    
        # Validate essential_count <= max_papers
        if ('essential_count' in topic_data and 'max_papers' in topic_data and
            isinstance(topic_data['essential_count'], int) and 
            isinstance(topic_data['max_papers'], int)):
            if topic_data['essential_count'] > topic_data['max_papers']:
                errors.append("essential_count cannot be greater than max_papers")
                
        return errors
        
    @staticmethod
    def validate_github_config(github_data: Dict[str, Any]) -> List[str]:
        """Validate GitHub configuration.
        
        Args:
            github_data: Dictionary containing GitHub configuration
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        
        # Validate repository format
        if 'repository' not in github_data:
            errors.append("Missing required field: repository")
        else:
            repo = github_data['repository']
            if not repo:
                errors.append("Repository cannot be empty")
            elif '/' not in repo:
                errors.append("Repository must be in format 'owner/repo'")
            elif repo.count('/') != 1:
                errors.append("Repository format invalid. Use 'owner/repo'")
                
        # Validate issue labels
        if 'issue_labels' in github_data:
            labels = github_data['issue_labels']
            if not isinstance(labels, list):
                errors.append("issue_labels must be a list")
            elif not labels:
                errors.append("At least one issue label is required")
            else:
                for label in labels:
                    if not isinstance(label, str):
                        errors.append("All issue labels must be strings")
                    elif not label.strip():
                        errors.append("Issue labels cannot be empty")
                        
        return errors
        
    @staticmethod
    def validate_scoring_weights(weights_data: Dict[str, Any]) -> List[str]:
        """Validate scoring weights configuration.
        
        Args:
            weights_data: Dictionary containing scoring weights
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        
        required_weights = [
            'citation_weight',
            'impact_factor_weight', 
            'recency_weight',
            'relevance_weight'
        ]
        
        # Check all weights are present and numeric
        weights = {}
        for weight_name in required_weights:
            if weight_name not in weights_data:
                errors.append(f"Missing required weight: {weight_name}")
            else:
                value = weights_data[weight_name]
                if not isinstance(value, (int, float)):
                    errors.append(f"Weight '{weight_name}' must be a number")
                elif value < 0:
                    errors.append(f"Weight '{weight_name}' cannot be negative")
                else:
                    weights[weight_name] = float(value)
                    
        # Check weights sum to 1.0 (with small tolerance for floating point errors)
        if len(weights) == len(required_weights):
            total = sum(weights.values())
            if abs(total - 1.0) > 0.001:
                errors.append(f"Weights must sum to 1.0, got {total:.3f}")
                
        return errors
        
    @staticmethod
    def validate_file_paths(config_dir: str) -> List[str]:
        """Validate that configuration directory and files exist.
        
        Args:
            config_dir: Path to configuration directory
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        
        config_path = Path(config_dir)
        if not config_path.exists():
            errors.append(f"Configuration directory does not exist: {config_dir}")
            return errors
            
        if not config_path.is_dir():
            errors.append(f"Configuration path is not a directory: {config_dir}")
            return errors
            
        # Check for required files
        required_files = ['topics.yaml', 'settings.yaml']
        for filename in required_files:
            file_path = config_path / filename
            if not file_path.exists():
                errors.append(f"Required configuration file missing: {filename}")
            elif not file_path.is_file():
                errors.append(f"Configuration path is not a file: {filename}")
                
        return errors


class DataValidator:
    """Validates data objects and API responses."""
    
    @staticmethod
    def validate_pmid(pmid: Union[str, int]) -> bool:
        """Validate PubMed ID format.
        
        Args:
            pmid: PubMed ID to validate
            
        Returns:
            True if valid PMID format
        """
        if not pmid:
            return False
            
        pmid_str = str(pmid).strip()
        
        # PMID should be numeric and reasonable length
        if not pmid_str.isdigit():
            return False
            
        # PMIDs are typically 1-8 digits
        if len(pmid_str) < 1 or len(pmid_str) > 8:
            return False
            
        return True
        
    @staticmethod
    def validate_doi(doi: str) -> bool:
        """Validate DOI format.
        
        Args:
            doi: DOI to validate
            
        Returns:
            True if valid DOI format
        """
        if not doi:
            return False
            
        doi = doi.strip()
        
        # Basic DOI pattern: 10.xxxx/yyyy
        doi_pattern = r'^10\.\d{4,}/[^\s]+$'
        return bool(re.match(doi_pattern, doi))
        
    @staticmethod
    def validate_email(email: str) -> bool:
        """Validate email address format.
        
        Args:
            email: Email address to validate
            
        Returns:
            True if valid email format
        """
        if not email:
            return False
            
        email_pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
        return bool(re.match(email_pattern, email.strip()))
        
    @staticmethod
    def validate_url(url: str) -> bool:
        """Validate URL format.
        
        Args:
            url: URL to validate
            
        Returns:
            True if valid URL format
        """
        if not url:
            return False
            
        try:
            result = urlparse(url.strip())
            return all([result.scheme, result.netloc])
        except Exception:
            return False
            
    @staticmethod
    def validate_github_token(token: str) -> bool:
        """Validate GitHub token format.
        
        Args:
            token: GitHub token to validate
            
        Returns:
            True if valid token format
        """
        if not token:
            return False
            
        token = token.strip()
        
        # GitHub tokens are typically 40 characters (classic) or start with specific prefixes
        if len(token) == 40 and re.match(r'^[a-f0-9]{40}$', token):
            return True
            
        # New format tokens
        if token.startswith(('ghp_', 'gho_', 'ghu_', 'ghs_', 'ghr_')):
            return len(token) > 10
            
        return False
        
    @staticmethod
    def validate_paper_data(paper_data: Dict[str, Any]) -> List[str]:
        """Validate paper data structure.
        
        Args:
            paper_data: Dictionary containing paper data
            
        Returns:
            List of validation error messages
        """
        errors = []
        
        # Required fields
        required_fields = ['pmid', 'title', 'authors', 'journal']
        for field in required_fields:
            if field not in paper_data:
                errors.append(f"Missing required field: {field}")
            elif not paper_data[field]:
                errors.append(f"Field '{field}' cannot be empty")
                
        # Validate PMID
        if 'pmid' in paper_data:
            if not DataValidator.validate_pmid(paper_data['pmid']):
                errors.append(f"Invalid PMID format: {paper_data['pmid']}")
                
        # Validate DOI if present
        if 'doi' in paper_data and paper_data['doi']:
            if not DataValidator.validate_doi(paper_data['doi']):
                errors.append(f"Invalid DOI format: {paper_data['doi']}")
                
        # Validate publication date if present
        if 'publication_date' in paper_data and paper_data['publication_date']:
            if not isinstance(paper_data['publication_date'], datetime):
                try:
                    datetime.fromisoformat(str(paper_data['publication_date']))
                except ValueError:
                    errors.append(f"Invalid publication date format: {paper_data['publication_date']}")
                    
        # Validate authors list
        if 'authors' in paper_data:
            authors = paper_data['authors']
            if isinstance(authors, list):
                if not authors:
                    errors.append("Authors list cannot be empty")
                else:
                    for i, author in enumerate(authors):
                        if not isinstance(author, str) or not author.strip():
                            errors.append(f"Author at index {i} must be a non-empty string")
            elif isinstance(authors, str):
                if not authors.strip():
                    errors.append("Authors string cannot be empty")
            else:
                errors.append("Authors must be a list or string")
                
        return errors
        
    @staticmethod
    def validate_citation_count(count: Union[int, str]) -> bool:
        """Validate citation count value.
        
        Args:
            count: Citation count to validate
            
        Returns:
            True if valid citation count
        """
        try:
            count_int = int(count)
            return count_int >= 0
        except (ValueError, TypeError):
            return False
            
    @staticmethod
    def validate_impact_factor(factor: Union[float, str]) -> bool:
        """Validate impact factor value.
        
        Args:
            factor: Impact factor to validate
            
        Returns:
            True if valid impact factor
        """
        try:
            factor_float = float(factor)
            return 0.0 <= factor_float <= 100.0  # Reasonable range for impact factors
        except (ValueError, TypeError):
            return False
            
    @staticmethod
    def validate_score(score: Union[float, int]) -> bool:
        """Validate scoring value.
        
        Args:
            score: Score to validate
            
        Returns:
            True if valid score
        """
        try:
            score_float = float(score)
            return 0.0 <= score_float <= 100.0
        except (ValueError, TypeError):
            return False


class EnvironmentValidator:
    """Validates environment variables and system requirements."""
    
    @staticmethod
    def validate_required_env_vars(required_vars: List[str]) -> List[str]:
        """Validate that required environment variables are set.
        
        Args:
            required_vars: List of required environment variable names
            
        Returns:
            List of missing environment variables
        """
        missing_vars = []
        for var in required_vars:
            if not os.getenv(var):
                missing_vars.append(var)
        return missing_vars
        
    @staticmethod
    def validate_github_environment() -> List[str]:
        """Validate GitHub-related environment variables.
        
        Returns:
            List of validation error messages
        """
        errors = []
        
        # Check for GitHub token
        token = os.getenv('GITHUB_TOKEN')
        if not token:
            errors.append("GITHUB_TOKEN environment variable is required")
        elif not DataValidator.validate_github_token(token):
            errors.append("GITHUB_TOKEN format appears invalid")
            
        # Check for repository context in GitHub Actions
        if os.getenv('GITHUB_ACTIONS'):
            repo = os.getenv('GITHUB_REPOSITORY')
            if not repo:
                errors.append("GITHUB_REPOSITORY not set in GitHub Actions context")
            elif '/' not in repo:
                errors.append("GITHUB_REPOSITORY format invalid")
                
        return errors
        
    @staticmethod
    def validate_pubmed_environment() -> List[str]:
        """Validate PubMed-related environment variables.
        
        Returns:
            List of validation error messages
        """
        errors = []
        
        # Check for PubMed email (recommended by NCBI)
        email = os.getenv('PUBMED_EMAIL')
        if email and not DataValidator.validate_email(email):
            errors.append("PUBMED_EMAIL format appears invalid")
            
        return errors
        
    @staticmethod
    def validate_system_requirements() -> List[str]:
        """Validate system requirements and dependencies.
        
        Returns:
            List of validation error messages
        """
        errors = []
        
        # Check Python version
        import sys
        if sys.version_info < (3, 8):
            errors.append(f"Python 3.8+ required, found {sys.version}")
            
        # Check for required directories
        required_dirs = ['logs', 'cache']
        for dir_name in required_dirs:
            dir_path = Path(dir_name)
            if not dir_path.exists():
                try:
                    dir_path.mkdir(exist_ok=True)
                except Exception as e:
                    errors.append(f"Cannot create required directory '{dir_name}': {e}")
                    
        return errors


# Utility functions for validation
def validate_and_raise(
    data: Any,
    validator: Callable[[Any], bool],
    error_message: str,
    field_name: Optional[str] = None
) -> None:
    """Validate data and raise ValidationError if invalid.
    
    Args:
        data: Data to validate
        validator: Validation function
        error_message: Error message if validation fails
        field_name: Name of the field being validated
        
    Raises:
        ValidationError: If validation fails
    """
    if not validator(data):
        raise ValidationError(error_message, field=field_name, value=data)


def validate_json_structure(
    data: Dict[str, Any],
    required_fields: List[str],
    optional_fields: Optional[List[str]] = None
) -> List[str]:
    """Validate JSON structure against required and optional fields.
    
    Args:
        data: Dictionary to validate
        required_fields: List of required field names
        optional_fields: List of optional field names
        
    Returns:
        List of validation error messages
    """
    errors = []
    
    # Check required fields
    for field in required_fields:
        if field not in data:
            errors.append(f"Missing required field: {field}")
        elif data[field] is None:
            errors.append(f"Required field '{field}' cannot be null")
            
    # Check for unexpected fields
    if optional_fields is not None:
        allowed_fields = set(required_fields + optional_fields)
        for field in data.keys():
            if field not in allowed_fields:
                errors.append(f"Unexpected field: {field}")
                
    return errors


def sanitize_string(
    value: str,
    max_length: Optional[int] = None,
    allowed_chars: Optional[str] = None,
    strip_whitespace: bool = True
) -> str:
    """Sanitize and validate string input.
    
    Args:
        value: String to sanitize
        max_length: Maximum allowed length
        allowed_chars: Regex pattern for allowed characters
        strip_whitespace: Whether to strip leading/trailing whitespace
        
    Returns:
        Sanitized string
        
    Raises:
        ValidationError: If string fails validation
    """
    if not isinstance(value, str):
        raise ValidationError("Value must be a string", value=value)
        
    if strip_whitespace:
        value = value.strip()
        
    if max_length and len(value) > max_length:
        raise ValidationError(f"String too long (max {max_length} chars)", value=value)
        
    if allowed_chars and not re.match(allowed_chars, value):
        raise ValidationError(f"String contains invalid characters", value=value)
        
    return value


def normalize_journal_name(journal_name: str) -> str:
    """Normalize journal name for consistent matching.
    
    Args:
        journal_name: Raw journal name
        
    Returns:
        Normalized journal name
    """
    if not journal_name:
        return ""
        
    # Basic normalization
    normalized = journal_name.strip()
    
    # Remove common prefixes/suffixes
    prefixes_to_remove = ['The ', 'A ', 'An ']
    for prefix in prefixes_to_remove:
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix):]
            
    # Remove trailing periods and whitespace
    normalized = normalized.rstrip('. ')
    
    # Convert to title case for consistency
    normalized = normalized.title()
    
    return normalized


def validate_batch_data(
    items: List[Dict[str, Any]],
    validator_func: Callable[[Dict[str, Any]], List[str]],
    max_errors: int = 10
) -> Dict[str, Any]:
    """Validate a batch of data items.
    
    Args:
        items: List of data items to validate
        validator_func: Function to validate each item
        max_errors: Maximum number of errors to collect
        
    Returns:
        Dictionary with validation results
    """
    results = {
        'valid_items': [],
        'invalid_items': [],
        'total_items': len(items),
        'error_count': 0,
        'errors': []
    }
    
    for i, item in enumerate(items):
        try:
            errors = validator_func(item)
            if errors:
                results['invalid_items'].append({
                    'index': i,
                    'item': item,
                    'errors': errors
                })
                results['error_count'] += len(errors)
                results['errors'].extend([f"Item {i}: {error}" for error in errors])
                
                # Stop if too many errors
                if len(results['errors']) >= max_errors:
                    results['errors'].append(f"... and more (stopped at {max_errors} errors)")
                    break
            else:
                results['valid_items'].append(item)
                
        except Exception as e:
            results['invalid_items'].append({
                'index': i,
                'item': item,
                'errors': [f"Validation exception: {e}"]
            })
            results['error_count'] += 1
            
    results['success_rate'] = len(results['valid_items']) / results['total_items'] if results['total_items'] > 0 else 0
    
    return results