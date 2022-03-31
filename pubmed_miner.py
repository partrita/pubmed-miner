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
        return refine_df
    except:
        return []


if __name__ == "__main__":
    entrez_query = input("What do want to search in PubMed? ")
    pmids = get_pmid(entrez_query, "100")  # string should be used for numberv
    result_df = pd.DataFrame()
    print(f"PubMed search for {entrez_query} is Done!")
    for i in tqdm(pmids, unit=" pmid"):
        # print(query_pmid(i))
        new_df = make_table(query_pmid([i]))  # only list works!
        time.sleep(1)
        try:
            # print(result_df)
            # print(new_df)
            result_df = pd.concat([result_df, new_df], ignore_index=True)
        except:
            pass

    print("Concet the dataframe is Done")
    result_df.to_csv(f"{entrez_query}.csv")
    print("Results saved")
