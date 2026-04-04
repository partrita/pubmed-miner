#!/usr/bin/env python3
"""
Setup script for Essential Papers Collection automation.

This script helps users configure the automation system by:
1. Validating configuration files
2. Testing GitHub API access
3. Providing setup instructions
"""

import os
import sys
from pathlib import Path

# Add src to path for imports
# Path(__file__).parent is scripts/
# Path(__file__).parent.parent is root/
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pubmed_miner.utils.config_manager import ConfigurationManager
from pubmed_miner.services.github_manager import GitHubIssuesManager


def check_environment():
    """Check if required environment variables are set."""
    print("🔍 Checking environment variables...")

    recommended_vars = {
        "GITHUB_TOKEN": "GitHub personal access token for API access (will use mock mode if not set)"
    }

    optional_vars = {
        "PUBMED_EMAIL": "Email address for PubMed API (defaults to example)"
    }

    for var, description in recommended_vars.items():
        if os.getenv(var):
            print(f"  ✅ {var}: Set")
        else:
            print(f"  ⚠️  {var}: Not set - {description}")

    for var, description in optional_vars.items():
        if os.getenv(var):
            print(f"  ✅ {var}: Set")
        else:
            print(f"  ⚠️  {var}: Not set - {description}")

    return (
        True  # Always return True since no vars are strictly required for local testing
    )


def validate_configuration():
    """Validate configuration files."""
    print("\n📋 Validating configuration files...")

    try:
        config_manager = ConfigurationManager()

        # Check if config files exist
        if not config_manager.topics_file.exists():
            print(f"  ❌ Topics file missing: {config_manager.topics_file}")
            return False

        if not config_manager.settings_file.exists():
            print(f"  ❌ Settings file missing: {config_manager.settings_file}")
            return False

        # Validate configuration
        config_manager.validate_config()
        print("  ✅ Configuration files are valid")

        # Show loaded topics
        topics = config_manager.load_topics()
        enabled_topics = [t for t in topics if t.enabled]
        print(f"  📊 Found {len(topics)} topics ({len(enabled_topics)} enabled)")

        for topic in enabled_topics:
            print(f"    - {topic.name}: {topic.query}")

        return True

    except Exception as e:
        print(f"  ❌ Configuration validation failed: {e}")
        return False


def test_github_access():
    """Test GitHub API access."""
    print("\n🐙 Testing GitHub API access...")

    if not os.getenv("GITHUB_TOKEN"):
        print("  ⚠️  GITHUB_TOKEN not set - will use mock mode for local testing")
        return True  # Return True for local testing

    try:
        config_manager = ConfigurationManager()
        github_config = config_manager.get_github_settings()

        # Test basic API access
        github_manager = GitHubIssuesManager(github_config)

        # Check if we're in mock mode
        if github_manager.mock_mode:
            print("  ⚠️  Running in mock mode - GitHub features will be simulated")
            return True

        # Try to get repository info (this tests authentication)
        import requests

        response = requests.get(
            f"https://api.github.com/repos/{github_config.repository}",
            headers=github_manager.headers,
            timeout=10,
        )

        if response.status_code == 200:
            repo_info = response.json()
            print(
                f"  ✅ Successfully connected to repository: {repo_info['full_name']}"
            )
            print(f"  📊 Repository has {repo_info['open_issues_count']} open issues")
            return True
        elif response.status_code == 404:
            print(f"  ❌ Repository not found: {github_config.repository}")
            print("     Make sure the repository exists and the token has access")
            return False
        else:
            print(f"  ❌ GitHub API error: {response.status_code} - {response.text}")
            return False

    except Exception as e:
        print(f"  ❌ GitHub API test failed: {e}")
        return False


def show_setup_instructions():
    """Show setup instructions for users."""
    print("\n📖 Setup Instructions:")
    print("=" * 50)

    print("\n1. 🔑 GitHub Token Setup:")
    print("   - Go to GitHub Settings > Developer settings > Personal access tokens")
    print("   - Create a new token with 'repo' and 'issues' permissions")
    print("   - Set the token as GITHUB_TOKEN environment variable or GitHub secret")

    print("\n2. 📧 PubMed Email (Optional but Recommended):")
    print("   - Method 1: Set PUBMED_EMAIL environment variable")
    print("   - Method 2: Edit config/settings.yaml and set pubmed.email")
    print("   - This is required by NCBI for API access")
    print("   - Environment variable takes priority over config file")

    print("\n3. ⚙️ Configuration Files:")
    print("   - Edit config/topics.yaml to add your research topics")
    print("   - Edit config/settings.yaml to configure GitHub repository")

    print("\n4. 🚀 GitHub Actions:")
    print(
        "   - The workflow is already configured in .github/workflows/collect-papers.yml"
    )
    print("   - It will run daily at 6:00 AM UTC")
    print("   - You can also trigger it manually from the Actions tab")

    print("\n5. 🧪 Testing:")
    print("   - Run this script again to validate your setup")
    print("   - Test manually with: python scripts/automated_collection.py")


def main():
    """Main setup function."""
    print("🔬 Essential Papers Collection - Setup Wizard")
    print("=" * 50)

    # Check environment
    env_ok = check_environment()

    # Validate configuration
    config_ok = validate_configuration()

    # Test GitHub access (only if environment is OK)
    github_ok = test_github_access() if env_ok else False

    # Show results
    print("\n📊 Setup Status:")
    print(f"  Environment: {'✅ OK' if env_ok else '❌ Issues'}")
    print(f"  Configuration: {'✅ OK' if config_ok else '❌ Issues'}")
    print(f"  GitHub Access: {'✅ OK' if github_ok else '❌ Issues'}")

    if env_ok and config_ok and github_ok:
        print("\n🎉 Setup complete! Your automation is ready to run.")
        print("   You can test it manually with: python scripts/automated_collection.py")
    else:
        print("\n⚠️  Setup incomplete. Please address the issues above.")
        show_setup_instructions()

    return env_ok and config_ok and github_ok


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
