import pytest
import tempfile
import csv
from pathlib import Path
from src.pubmed_miner.data.journal_database import JournalDatabase

def test_export_to_csv_sanitizes_csv_injection():
    """Verify that export_to_csv sanitizes payloads with leading whitespace."""
    db = JournalDatabase()

    # Add a journal with malicious data
    db.add_journal(
        name="  =cmd|' /C calc'!A0",  # CSV injection payload with leading space
        impact_factor=5.0,
        category="\t+Category",      # CSV injection payload with leading tab
        publisher=" \n -Malicious",   # CSV injection payload with leading newline
        issn="   @1234-5678"         # CSV injection payload
    )

    # Export to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        temp_path = f.name

    try:
        db.export_to_csv(temp_path)

        # Read the file and verify sanitization
        with open(temp_path, "r", newline="", encoding="utf-8") as csvfile:
            reader = csv.DictReader(csvfile)
            rows = list(reader)

            # Find our malicious row
            row = next(r for r in rows if "calc" in r["name"])
            assert row["name"] == "'=cmd|' /c calc'!a0"  # Normalized to lower case by add_journal and spaces stripped by normalization
            assert row["category"] == "'\t+Category"
            assert row["publisher"] == "' \n -Malicious"
            assert row["issn"] == "'   @1234-5678"
    finally:
        Path(temp_path).unlink(missing_ok=True)
