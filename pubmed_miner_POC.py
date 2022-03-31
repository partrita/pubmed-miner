from Bio import Entrez
import requests
import json
import pandas as pd

Entrez.email = "youremail@dot.com"  # Always tell NCBI who you are
url = "http://bern2.korea.ac.kr/pubmed"  # Remember 100 requests for 100 seconds


def get_pmid(query, num_list):
    handle = Entrez.esearch(
        db="pubmed", term=query, retmax=num_list, sort="pub+date", retmode="xml"
    )
    records = Entrez.read(handle)
    return records["IdList"]


def query_pmid(pmids):
    return requests.get(url + "/" + ",".join(pmids)).json()


def make_table(json_list):
    meta = ["pmid"]
    df = pd.json_normalize(json_list, record_path=["annotations"], meta=meta)
    refine_df = df[["mention", "obj", "prob", "pmid"]]
    return refine_df


if __name__ == "__main__":
    entrez_query = input("What do want to search in PubMed? ")
    pmids = get_pmid(entrez_query, "5")
    print(pmids)
    df = make_table(query_pmid(pmids))
    print(df)
    df.to_csv(f"{entrez_query}.csv")
    print("Results saved")
