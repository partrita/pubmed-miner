import pytest
import pandas as pd
from src.pubmed_miner.core import make_table

# Expected columns for the DataFrame
EXPECTED_COLUMNS = ['mention', 'obj', 'prob', 'pmid']

def test_make_table_valid_input():
    """Test make_table with a valid list of JSON objects."""
    sample_json_list = [
        {
            "pmid": "12345",
            "annotations": [
                {"mention": "Aspirin", "obj": "drug", "prob": 0.9},
                {"mention": "Cancer", "obj": "disease", "prob": 0.85},
            ],
        },
        {
            "pmid": "67890",
            "annotations": [
                {"mention": "Diabetes", "obj": "disease", "prob": 0.92},
            ],
        },
    ]
    df = make_table(sample_json_list)

    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    assert len(df) == 3, "DataFrame should have 3 rows (2 from first PMID, 1 from second)"

    # Check data integrity for the first annotation of the first PMID
    assert df.iloc[0]['mention'] == "Aspirin"
    assert df.iloc[0]['obj'] == "drug"
    assert df.iloc[0]['prob'] == 0.9
    assert df.iloc[0]['pmid'] == 12345 # pmid should be int after astype

    # Check data integrity for the second annotation of the first PMID
    assert df.iloc[1]['mention'] == "Cancer"
    assert df.iloc[1]['obj'] == "disease"
    assert df.iloc[1]['prob'] == 0.85
    assert df.iloc[1]['pmid'] == 12345

    # Check data integrity for the annotation of the second PMID
    assert df.iloc[2]['mention'] == "Diabetes"
    assert df.iloc[2]['obj'] == "disease"
    assert df.iloc[2]['prob'] == 0.92
    assert df.iloc[2]['pmid'] == 67890

def test_make_table_empty_input_list():
    """Test make_table with an empty list."""
    df = make_table([])
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    assert len(df) == 0, "DataFrame should be empty (0 rows)"

def test_make_table_input_none():
    """Test make_table with None as input."""
    df = make_table(None)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    assert len(df) == 0, "DataFrame should be empty (0 rows)"

def test_make_table_malformed_input_missing_annotations():
    """Test make_table with input missing the 'annotations' key."""
    malformed_json_list = [
        {"pmid": "11111"}, # Missing 'annotations'
        {"pmid": "22222", "annotations": []}, # Empty 'annotations'
    ]
    df = make_table(malformed_json_list)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"

    # Case 1: 'annotations' key is missing - this item will be skipped by the `isinstance(annotations, list)` check
    # or by the initial `if not isinstance(item, dict) or 'pmid' not in item:` check if pmid is also missing.
    # If pmid is present but annotations is missing, `item.get("annotations")` is None, so it's skipped.
    # Case 2: 'annotations' is an empty list - this will create a row with NAs for annotation fields.
    # The test input has one of each.
    # The first entry {"pmid": "11111"} will be skipped (annotations is None).
    # The second entry {"pmid": "22222", "annotations": []} will produce one row with NAs.
    assert len(df) == 1, "DataFrame should have one row for the entry with empty annotations"
    assert df.iloc[0]['pmid'] == 22222
    assert pd.isna(df.iloc[0]['mention'])
    assert pd.isna(df.iloc[0]['obj'])
    assert pd.isna(df.iloc[0]['prob'])

def test_make_table_malformed_input_no_pmid_in_meta():
    """Test make_table with input where pmid is not in the main dict (item skipped)."""
    # This tests robustness if the structure isn't what's expected by pd.json_normalize's meta.
    malformed_json_list = [
        {
            "annotations": [
                {"mention": "Paracetamol", "obj": "drug", "prob": 0.95, "pmid": "33333"},
            ]
        }
    ]
    df = make_table(malformed_json_list)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    # The current make_table skips items if 'pmid' is not a top-level key in the item dictionary.
    df = make_table(malformed_json_list)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    assert len(df) == 0, "DataFrame should be empty as the item is skipped due to missing top-level pmid"

def test_make_table_annotations_not_a_list():
    """Test make_table with 'annotations' being a dictionary instead of a list."""
    malformed_json_list = [
        {
            "pmid": "44444",
            "annotations": {"mention": "Flu", "obj": "disease", "prob": 0.7} # Should be a list of dicts
        }
    ]
    df = make_table(malformed_json_list)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    # pd.json_normalize might raise an error or return an empty df for this record_path
    # The error handling in make_table should catch this and return an empty df with expected columns
    assert len(df) == 0, "DataFrame should be empty if 'annotations' is not a list"

def test_make_table_mixed_valid_and_malformed_input():
    """Test make_table with a mix of valid and malformed entries."""
    mixed_json_list = [
        { # Valid entry
            "pmid": "12345",
            "annotations": [{"mention": "Aspirin", "obj": "drug", "prob": 0.9}],
        },
        {"pmid": "55555"}, # Malformed: Missing 'annotations'
        { # Valid entry
            "pmid": "67890",
            "annotations": [{"mention": "Diabetes", "obj": "disease", "prob": 0.92}],
        },
        { # Malformed: 'annotations' is not a list
            "pmid": "44444",
            "annotations": {"mention": "Flu", "obj": "disease", "prob": 0.7}
        }
    ]
    df = make_table(mixed_json_list)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    # Expect 2 rows from the valid entries. Malformed entries should be skipped or result in no rows.
    assert len(df) == 2, "DataFrame should contain rows from valid entries only"
    assert df.iloc[0]['pmid'] == 12345
    assert df.iloc[1]['pmid'] == 67890
    assert df.iloc[0]['mention'] == "Aspirin"
    assert df.iloc[1]['mention'] == "Diabetes"

def test_make_table_annotations_with_missing_fields():
    """Test make_table where annotations have missing optional fields (e.g. prob)."""
    json_list_missing_fields = [
        {
            "pmid": "77777",
            "annotations": [
                {"mention": "Headache", "obj": "symptom"}, # Missing 'prob'
            ],
        }
    ]
    df = make_table(json_list_missing_fields)
    assert isinstance(df, pd.DataFrame), "Should return a pandas DataFrame"
    assert list(df.columns) == EXPECTED_COLUMNS, "DataFrame should have the expected columns"
    assert len(df) == 1, "DataFrame should have 1 row"
    assert df.iloc[0]['mention'] == "Headache"  # Will be 'Headache' string
    assert df.iloc[0]['obj'] == "symptom"    # Will be 'symptom' string
    assert pd.isna(df.iloc[0]['prob']), "Missing 'prob' should result in NA"
    assert df.iloc[0]['pmid'] == 77777 # Compare with int
