from Bio import Entrez
import requests
import json, time
import pandas as pd
from tqdm import tqdm
from typing import List, Dict, Optional, Union
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Constants
ENTREZ_EMAIL = "youremail@dot.com"  # Always tell NCBI who you are
BERN2_API_URL = "http://bern2.korea.ac.kr/pubmed"  # Remember 100 requests for 100 seconds
Entrez.email = ENTREZ_EMAIL


def get_pmid(query: str, num_records: str) -> List[str]:
    """
    Retrieves a list of PMIDs for a given query.
    num_records should be below 10000.
    """
    try:
        handle = Entrez.esearch(
            db="pubmed", term=query, retmax=num_records, sort="pub+date", retmode="xml"
        )
        records: Dict = Entrez.read(handle)
        handle.close()
        return records.get("IdList", [])
    except Exception as e:
        logging.error(f"Error in get_pmid for query '{query}': {e}")
        return []


def query_pmid(pmids: List[str], url: str = BERN2_API_URL) -> Optional[List[Dict]]:
    # headers = {"Content-Type": "application/json", "charset": "UTF-8", "Accept": "*/*"}
    if not pmids:
        return None
    try:
        response = requests.get(url + "/" + ",".join(pmids))
        response.raise_for_status()  # Raises HTTPError for bad responses (4XX or 5XX)
        return response.json()
    except requests.exceptions.RequestException as e:
        logging.error(f"Error querying PMIDs {pmids}: {e}")
        return None
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON response for PMIDs {pmids}: {e}")
        return None


def make_table(json_list: Optional[List[Dict]]) -> pd.DataFrame:
    """
    Converts a list of JSON objects (from BERN2 API) into a pandas DataFrame.
    """
    expected_columns = ['mention', 'obj', 'prob', 'pmid']
    all_records_df = pd.DataFrame(columns=expected_columns)

    if not json_list:
        return all_records_df

    processed_dfs = []
    for item in json_list:
        try:
            # Ensure item is a dictionary and 'pmid' is present for meta
            if not isinstance(item, dict) or 'pmid' not in item:
                # logging.warning(f"Skipping item due to missing 'pmid' or not a dict: {item}")
                continue

            # Attempt to normalize this single item
            # 'annotations' must be a list for record_path to work as expected.
            # If 'annotations' is not a list or missing, json_normalize with errors='ignore'
            # on a single item might return an empty df or df without 'annotations' fields.
            annotations = item.get("annotations")
            if not isinstance(annotations, list):
                # logging.warning(f"Skipping item {item.get('pmid')} due to 'annotations' not being a list: {annotations}")
                # If we want to create rows with NAs for mention, obj, prob for this pmid:
                # temp_df = pd.DataFrame([{'pmid': item.get('pmid')}]) 
                # for col in ['mention', 'obj', 'prob']: temp_df[col] = pd.NA
                # processed_dfs.append(temp_df[expected_columns]) # ensure order and all columns
                continue # Or just skip this item

            if not annotations: # Empty annotations list
                 # Create a record with NA for annotation fields but with the pmid
                temp_df = pd.DataFrame([{'pmid': item.get('pmid')}])
                for col in ['mention', 'obj', 'prob']:
                    temp_df[col] = pd.NA
                # Ensure all expected_columns are present, especially if 'pmid' was the only one.
                for col in expected_columns:
                    if col not in temp_df.columns:
                         temp_df[col] = pd.NA
                processed_dfs.append(temp_df[expected_columns])
                continue
            
            df_item = pd.json_normalize(item, record_path=["annotations"], meta=["pmid"], errors='raise') # Raise error to catch it below
            
            # Ensure all expected columns are present after normalization, fill with NA if not
            for col in expected_columns:
                if col not in df_item.columns:
                    df_item[col] = pd.NA
            
            processed_dfs.append(df_item[expected_columns])

        except (AttributeError, KeyError, TypeError, ValueError, RuntimeError) as e: # Added RuntimeError for potential pandas issues
            logging.warning(f"Skipping item {item.get('pmid', 'Unknown PMID')} due to error: {e}. Item: {item}")
            # Optionally, create a row with NA for this pmid if it's important to represent it
            # temp_df = pd.DataFrame([{'pmid': item.get('pmid')}])
            # for col in ['mention', 'obj', 'prob']: temp_df[col] = pd.NA
            # processed_dfs.append(temp_df[expected_columns])
            continue # Skip this item on error

    if not processed_dfs:
        return all_records_df # Returns empty DataFrame with expected_columns

    all_records_df = pd.concat(processed_dfs, ignore_index=True)

    # Consolidate type conversion after concatenation
    # Ensure 'pmid' column exists before trying to convert it
    if 'pmid' not in all_records_df.columns:
        all_records_df['pmid'] = pd.NA # Add pmid column if somehow missing
        
    # Convert types, be robust to existing NAs
    try:
        # Use pd.to_numeric for 'prob' for better NA handling before converting to float16
        if 'prob' in all_records_df.columns:
            all_records_df['prob'] = pd.to_numeric(all_records_df['prob'], errors='coerce')
        
        all_records_df = all_records_df.astype(
            {
                "prob": "float16", 
                "pmid": "Int32",  # Use nullable Int32 type to support NA alongside integers
                "mention": "string", # Use pandas string type to keep pd.NA as NA
                "obj": "string"      # Use pandas string type to keep pd.NA as NA
            }
        )

    except Exception as e:
        logging.error(f"Error during type conversion in make_table: {e}")
        # Return dataframe with unconverted types but correct columns if conversion fails
        return all_records_df[expected_columns]

    return all_records_df[expected_columns] # Ensure column order


def get_bern2(pmids: List[str]) -> pd.DataFrame:
    """
    Queries BERN2 API for a list of PMIDs and returns a concatenated DataFrame of annotations.
    """
    all_annotations_df = pd.DataFrame(columns=['mention', 'obj', 'prob', 'pmid'])
    if not pmids:
        return all_annotations_df

    for pmid_chunk in tqdm([pmids[i:i + 10] for i in range(0, len(pmids), 10)], desc="Querying BERN2 and processing results"): # Process in chunks for API efficiency
        queried_data = query_pmid(pmid_chunk) # query_pmid can now handle a list
        if queried_data:
            new_df = make_table(queried_data)
            if not new_df.empty:
                try:
                    all_annotations_df = pd.concat([all_annotations_df, new_df], ignore_index=True)
                except Exception as e: # Catch error during concat
                    logging.error(f"Error concatenating DataFrame for PMIDs {pmid_chunk}: {e}")
        time.sleep(1)  # delay for BERN2 API - consider if this is needed per chunk or per request
    return all_annotations_df


def main():
    """
    Main function to run the PubMed mining script.
    """
    entrez_query: str = input("What do you want to search in PubMed? ")
    max_records: str = "10000" # Define as a string as per get_pmid's original expectation

    pmids_list: List[str] = get_pmid(entrez_query, max_records)
    num_pmids: int = len(pmids_list)

    if not pmids_list:
        logging.info(f"No PMIDs found for query: {entrez_query}")
        return

    logging.info(f"Number of search results is: {num_pmids}!")

    results_df: pd.DataFrame = get_bern2(pmids_list)

    if results_df.empty:
        logging.info(f"No annotations found for the PMIDs from query: {entrez_query}")
    else:
        logging.info(f"Shape of result table is: {results_df.shape}.")
        
        output_dir = Path("./output")
        output_dir.mkdir(parents=True, exist_ok=True) # Ensure output directory exists
        
        # Sanitize entrez_query to be a valid filename
        safe_query_filename = "".join(c if c.isalnum() else "_" for c in entrez_query)
        if not safe_query_filename: # Handle empty or all-non-alphanumeric queries
            safe_query_filename = "pubmed_results"
        file_path = output_dir / f"{safe_query_filename}.csv"
        
        try:
            results_df.to_csv(file_path, index=False)
            logging.info(f"'{file_path}' file saved in output directory.")
        except Exception as e:
            logging.error(f"Error saving CSV file to '{file_path}': {e}")

    logging.info("Everything is completed! It is time for analysis.")


if __name__ == "__main__":
    main()
