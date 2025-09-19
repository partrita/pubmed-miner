#!/usr/bin/env python3
"""
Configuration validation script for GitHub Actions
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    from pubmed_miner.utils.config_manager import ConfigurationManager
    
    print("Validating configuration...")
    config_manager = ConfigurationManager()
    config_manager.validate_config()
    print("✅ Configuration validation successful")
    
except Exception as e:
    print(f"❌ Configuration validation failed: {e}")
    sys.exit(1)