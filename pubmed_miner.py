from Bio import Entrez
import requests
import json, time
import pandas as pd
from tqdm import tqdm

Entrez.email = "youremail@dot.com"  # Always tell NCBI who you are
url = "http://bern2.korea.ac.kr/pubmed"  # Remember 100 requests for 100 seconds


def get_pmid(query, num_list):
    """
    num_list should be below 10000
    """
    handle = Entrez.esearch(
        db="pubmed", term=query, retmax=num_list, sort="pub+date", retmode="xml"
    )
    records = Entrez.read(handle)
    return records["IdList"]


def query_pmid(pmids, url="http://bern2.korea.ac.kr/pubmed"):
    # headers = {"Content-Type": "application/json", "charset": "UTF-8", "Accept": "*/*"}
    try:
        return requests.get(url + "/" + ",".join(pmids)).json()
    except:
        pass


def make_table(json_list):
    try:
        df = pd.json_normalize(json_list, record_path=["annotations"], meta=["pmid"])
        refine_df = df[["mention", "obj", "prob", "pmid"]]
        # use pd.astype() function for save memory
        refine_df = refine_df.astype(
            {
                "prob": "float16",
                "pmid": "int32",
            }
        )
        return refine_df
    except:
        return []


def get_bern2(pmids):
    temp_df = pd.DataFrame()
    for i in tqdm(pmids, unit="pmid"):
        # print(query_pmid(i))
        new_df = make_table(query_pmid([i]))  # only list works!
        time.sleep(1)  # delay for BERN2 API
        try:
            temp_df = pd.concat([temp_df, new_df], ignore_index=True)
        except:
            pass
    return temp_df


if __name__ == "__main__":
    entrez_query = input("What do want to search in PubMed? ")
    pmids = get_pmid(entrez_query, "10000")  # string should be used for number
    num_pmids = len(pmids)
    print(f"Number of search results is: {num_pmids}!")
    temp_df = get_bern2(pmids)
    print(f"Shape of result table is: {temp_df.shape}.")
    temp_df.to_csv(f"./output/{entrez_query}.csv")
    print(f"{entrez_query}.csv fsile saved in output directory.")
    print(f"Everything is completed! It is time for analysis.")
