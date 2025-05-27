# PubMed Miner

`pubmed-miner` is a Python tool that retrieves PubMed article PMIDs based on a search query, analyzes these articles using the BERN2 API for biomedical named entity recognition and normalization, and then saves the structured annotations into a CSV file.

The core logic resides in `src/pubmed_miner/core.py`.

## Features

*   Fetches PMIDs from PubMed using Biopython's Entrez module.
*   Processes PMIDs through the BERN2 API to extract biomedical annotations.
*   Normalizes and structures the BERN2 output into a pandas DataFrame.
*   Saves the resulting table to a CSV file in the `output/` directory.

## Tools & APIs Used

*   **Biopython:** For interacting with the NCBI PubMed database.
    ```
    @article{cock2009biopython,
     title={Biopython: freely available Python tools for computational molecular biology and bioinformatics},
     author={Cock, Peter JA and Antao, Tiago and Chang, Jeffrey T and Chapman, Brad A and Cox, Cymon J and Dalke, Andrew and Friedberg, Iddo and Hamelryck, Thomas and Kauff, Frank and Wilczynski, Bartek and others},
     journal={Bioinformatics},
     volume={25},
     number={11},
     pages={1422--1423},
     year={2009},
     publisher={Oxford University Press}
    }
    ```
*   **BERN2 API:** For biomedical named entity recognition and normalization.
    ```
    @article{sung2022bern2,
     title={BERN2: an advanced neural biomedical namedentity recognition and normalization tool},
     author={Sung, Mujeen and Jeong, Minbyul and Choi, Yonghwa and Kim, Donghyeon and Lee, Jinhyuk and Kang, Jaewoo},
     year={2022},
     eprint={2201.02080},
     archivePrefix={arXiv},
     primaryClass={cs.CL}
    }
    ```
*   **Pandas:** For data manipulation and CSV output.
*   **uv:** For environment and package management.

## Setup and Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/pubmed-miner.git # Replace with actual URL
    cd pubmed-miner
    ```

2.  **Create a virtual environment and install dependencies using `uv`:**
    Ensure `uv` is installed (see [uv installation guide](https://astral.sh/uv/install.sh)).
    ```bash
    uv venv
    source .venv/bin/activate  # On Windows use: .venv\Scripts\activate
    uv pip install -e .[test]
    ```
    This installs the package in editable mode (`-e .`) along with its test dependencies (`[test]`).

## Usage

To run the script:

```bash
python src/pubmed_miner/core.py
```

The script will then prompt you to enter your PubMed search query. For example:
```
What do you want to search in PubMed? Aspirin and cancer
```

The script will:
1.  Fetch PMIDs related to "Aspirin and cancer".
2.  Query the BERN2 API with these PMIDs.
3.  Process the results and save them to a CSV file in the `output/` directory (e.g., `output/Aspirin_and_cancer.csv`).

**Note:** Remember to set your email address in `src/pubmed_miner/core.py` for the `Entrez.email` variable, as NCBI requires this for Entrez API access.
```python
# In src/pubmed_miner/core.py
ENTREZ_EMAIL = "your_actual_email@example.com"
Entrez.email = ENTREZ_EMAIL
```

## Running Tests

Tests are located in the `tests/` directory and can be run using `pytest`:

```bash
pytest tests/
```
This command will discover and run all tests in the specified directory.

## How it Works

1.  The user is prompted for a PubMed query via the command line.
2.  The `get_pmid` function (using Biopython's `Entrez` module) retrieves a list of PMIDs matching the query.
3.  The `get_bern2` function iterates through these PMIDs (in chunks), sending them to the BERN2 API via the `query_pmid` function.
4.  `query_pmid` makes a GET request to the BERN2 API.
5.  The JSON response from BERN2 is processed by `make_table` into a structured pandas DataFrame, focusing on annotations like mentions, object types (e.g., drug, disease), and confidence probabilities.
6.  The main script then saves this DataFrame to a CSV file in the `output/` directory, named after the original query.
```