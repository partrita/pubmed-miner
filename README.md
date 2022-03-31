# pubmed-miner

Python code for make table data of pubmed article analysis. What i do is just connect two tool.


# Tools that I use.

## Biopython

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

## BERN2 

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

## PDM for package management

more detail at https://pdm.fming.dev/

# Mode of action

```bash
python pubmed_miner.py
```

1. Ask user input for Pubmed query.
2. Biopython `entrez` module return bunch of article's pmid.
3. After that BERN2 used biological text analysis and return json file.
4. `pandas` used for trim the data and save `csv` file.