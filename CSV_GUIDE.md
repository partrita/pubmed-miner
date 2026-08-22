# CSV Collections 가이드

`CSVManager`를 사용하여 수집한 논문 정보를 CSV 파일로 저장하고 관리하는 방법을 설명합니다.

## 📚 개요

`CSVManager`는 PubMed에서 수집한 논문 데이터를 CSV 형식으로 저장, 로드, 추가하는 기능을 제공합니다.

### 주요 기능
- ✅ 논문 정보를 CSV로 저장
- ✅ 기존 CSV에 논문 추가 (append)
- ✅ CSV에서 논문 데이터 로드
- ✅ 점수 정보를 포함한 CSV 저장
- ✅ 자동 디렉토리 생성

---

## 🚀 기본 사용법

### 1. 기본 논문 저장하기

```python
from src.pubmed_miner.utils import CSVManager
from src.pubmed_miner.models import Paper
from datetime import datetime

# 논문 객체 생성
papers = [
    Paper(
        pmid="12345678",
        title="Machine Learning in Healthcare",
        authors=["John Doe", "Jane Smith"],
        journal="Nature Machine Intelligence",
        publication_date=datetime(2024, 1, 15),
        doi="10.1234/sample.2024.001",
        abstract="This paper discusses ML applications...",
        topic="ai-healthcare",  # Topic 정보 추가
    )
]

# CSV로 저장
CSVManager.save_papers(papers, "data/collections.csv")
```

### 2. 점수 정보와 함께 저장하기

```python
from src.pubmed_miner.models import ScoredPaper

scored_papers = [
    ScoredPaper(
        pmid="87654321",
        title="Deep Learning Research",
        authors=["Alice Johnson"],
        journal="Science",
        publication_date=datetime(2024, 1, 10),
        citation_count=150,
        impact_factor=42.5,
        score=95.0,
        rank=1,
        doi="10.1234/score.2024.001",
        topic="cancer-immunotherapy",  # Topic 정보 추가
    )
]

# 점수 정보 포함하여 저장
CSVManager.update_collection(scored_papers, "data/collections.csv")
```

### 3. 기존 CSV에 논문 추가하기

```python
# 새로운 논문 추가
new_papers = [
    Paper(
        pmid="99999999",
        title="New Research",
        authors=["New Author"],
        journal="New Journal",
        publication_date=datetime(2024, 1, 19),
        doi="10.1234/new.2024.001",
        topic="bioinformatics",  # Topic 정보 추가
    )
]

# append=True로 추가
CSVManager.append_papers(new_papers, "data/collections.csv")
```

### 4. CSV에서 논문 로드하기

```python
# CSV 파일 읽기
papers = CSVManager.load_papers("data/collections.csv")

# 데이터 확인
for paper in papers:
    print(f"PMID: {paper['pmid']}")
    print(f"Title: {paper['title']}")
    print(f"Authors: {paper['authors']}")
```

---

## 🔄 완전한 워크플로우 예제

PubMed에서 논문을 검색하여 CSV로 저장하는 전체 과정:

```python
from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.utils import CSVManager

# 1. 서비스 초기화
service = PaperCollectionService(email="your.email@example.com")

# 2. PubMed에서 논문 검색
query = "machine learning AND healthcare"
topic = "ai-healthcare"  # 주제 정보
pmids = service.search_papers(query, max_results=100)

# 3. 논문 상세 정보 조회 (topic 정보 포함)
papers = service.get_paper_details(pmids, topic=topic)

# 4. CSV로 저장
CSVManager.save_papers(papers, "data/collections.csv")

print(f"✓ {len(papers)}개의 논문을 저장했습니다.")
```

---

## 📋 CSV 파일 구조

### 기본 논문 (HEADERS)

```
pmid,title,authors,journal,publication_date,doi,abstract,topic
```

**컬럼 설명:**
- `pmid`: PubMed ID (필수)
- `title`: 논문 제목 (필수)
- `authors`: 저자 (세미콜론으로 구분, 예: "John Doe; Jane Smith")
- `journal`: 저널명 (필수)
- `publication_date`: 발행 날짜 (ISO format: YYYY-MM-DD)
- `doi`: DOI (선택)
- `abstract`: 초록 (선택)
- `topic`: 주제/토픽 (선택, 예: "ai-drug-discovery", "cancer-immunotherapy")

### 점수 포함 논문 (SCORED_HEADERS)

```
pmid,title,authors,journal,publication_date,doi,abstract,topic,citation_count,impact_factor,score,rank
```

**추가 컬럼:**
- `citation_count`: 인용 수
- `impact_factor`: 임팩트 팩터
- `score`: 점수
- `rank`: 순위

---

## 🛠️ API 참조

### `CSVManager.save_papers(papers, filepath, include_scoring=False, append=False)`

논문들을 CSV 파일로 저장합니다.

**파라미터:**
- `papers` (List[Paper]): 저장할 논문 목록
- `filepath` (str): CSV 파일 경로
- `include_scoring` (bool): 점수 정보 포함 여부 (기본값: False)
- `append` (bool): 기존 파일에 추가할지 여부 (기본값: False)

**예외:**
- `ValueError`: 빈 논문 목록
- `IOError`: 파일 쓰기 실패

### `CSVManager.append_papers(papers, filepath)`

기존 CSV 파일에 논문을 추가합니다.

**파라미터:**
- `papers` (List[Paper]): 추가할 논문 목록
- `filepath` (str): CSV 파일 경로

### `CSVManager.load_papers(filepath)`

CSV 파일에서 논문 데이터를 읽습니다.

**반환값:** List[Dict] - 논문 정보 딕셔너리 목록

**예외:**
- `FileNotFoundError`: 파일 미존재
- `IOError`: 파일 읽기 실패

### `CSVManager.update_collection(papers, filepath)`

점수 정보를 포함하여 CSV를 업데이트합니다.

**파라미터:**
- `papers` (List[ScoredPaper]): 점수가 포함된 논문 목록
- `filepath` (str): CSV 파일 경로

---

## 💡 실용적인 예제들

### 예제 1: 여러 검색어로 수집 및 저장 (토픽 포함)

```python
from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.utils import CSVManager

service = PaperCollectionService(email="your.email@example.com")

# 여러 주제로 검색
topics = {
    "ai-drug-discovery": "machine learning AND drug discovery",
    "cancer-immunotherapy": "cancer AND immunotherapy",
    "protein-design": "de novo design AND artificial intelligence AND protein design",
}

all_papers = []

for topic_name, query in topics.items():
    print(f"\n🔍 Collecting papers for {topic_name}...")
    pmids = service.search_papers(query, max_results=50)
    papers = service.get_paper_details(pmids, topic=topic_name)
    all_papers.extend(papers)

# 한번에 저장
CSVManager.save_papers(all_papers, "data/collections.csv")
print(f"✓ 총 {len(all_papers)}개 논문 저장 완료")
```

### 예제 2: 정기적으로 새 논문 추가

```python
from datetime import datetime

# 기존 데이터 로드
existing_papers = CSVManager.load_papers("data/collections.csv")
print(f"기존 논문: {len(existing_papers)}개")

# 새로운 논문 수집 (새로운 토픽)
service = PaperCollectionService(email="your.email@example.com")
new_topic = "bioinformatics"
pmids = service.search_papers("bioinformatics", max_results=20)
new_papers = service.get_paper_details(pmids, topic=new_topic)

# 추가 저장
if new_papers:
    CSVManager.append_papers(new_papers, "data/collections.csv")
    print(f"✓ {len(new_papers)}개 논문 추가됨 (topic: {new_topic})")
    print(f"총 {len(existing_papers) + len(new_papers)}개 논문")
```

### 예제 3: 토픽별 논문 분석

```python
from collections import Counter

# CSV 로드
papers = CSVManager.load_papers("data/collections.csv")

# 토픽별 논문 수
topics = [p["topic"] for p in papers if p["topic"]]
topic_counts = Counter(topics)

print("📊 토픽별 논문 수:")
for topic, count in topic_counts.most_common():
    print(f"  {topic}: {count}개")

# 토픽별 저널 분석
for topic in set(topics):
    topic_papers = [p for p in papers if p["topic"] == topic]
    journals = [p["journal"] for p in topic_papers]
    journal_counts = Counter(journals)
    print(f"\n📰 {topic} - 상위 저널:")
    for journal, count in journal_counts.most_common(3):
        print(f"  {journal}: {count}개")
```

### 예제 4: CSV 데이터 분석

```python
from collections import Counter

# CSV 로드
papers = CSVManager.load_papers("data/collections.csv")

# 저널별 논문 수
journals = [p["journal"] for p in papers]
journal_counts = Counter(journals)

print("📊 저널별 논문 수:")
for journal, count in journal_counts.most_common(10):
    print(f"  {journal}: {count}개")

# 저자 분석
all_authors = []
for paper in papers:
    authors = [a.strip() for a in paper["authors"].split(";")]
    all_authors.extend(authors)

author_counts = Counter(all_authors)
```

# 한번에 저장
CSVManager.save_papers(all_papers, "data/collections.csv")
print(f"✓ 총 {len(all_papers)}개 논문 저장 완료")
```

---

## ⚙️ 고급 옵션

### 데이터 형식 확인

CSV 헤더를 확인할 수 있습니다:

```python
# CSVManager.HEADERS 확인
print(CSVManager.HEADERS)
# ['pmid', 'title', 'authors', 'journal', 'publication_date', 'doi', 'abstract', 'topic']

print(CSVManager.SCORED_HEADERS)
# ['pmid', 'title', 'authors', 'journal', 'publication_date', 'doi', 'abstract', 'topic',
#  'citation_count', 'impact_factor', 'score', 'rank']
```

### 경로 처리

```python
from pathlib import Path

# 경로 자동 생성
filepath = Path("data/collections.csv")
filepath.parent.mkdir(parents=True, exist_ok=True)

# 또는 CSVManager가 자동으로 생성
CSVManager.save_papers(papers, "data/collections.csv")
```

---

## ✅ 테스트

```bash
# 모든 테스트 실행
pytest tests/unit/test_csv_manager.py -v

# 특정 테스트만 실행
pytest tests/unit/test_csv_manager.py::TestCSVManager::test_topic_field_included -v
```

---

## 📝 주의사항

1. **이메일 설정**: PubMed API 사용 시 유효한 이메일 주소 필수
2. **Rate Limiting**: PubMed는 요청 제한이 있으므로 대량 수집 시 시간 고려
3. **문자 인코딩**: CSV 파일은 UTF-8 인코딩으로 저장됨
4. **author 형식**: 여러 저자는 세미콜론(`;`)으로 구분됨
5. **날짜 형식**: ISO 8601 형식 (YYYY-MM-DD)
6. **topic 필드**: 선택 사항이지만 논문 분류를 위해 권장됨

---

## 🔗 관련 리소스

- [PubMed API 문서](https://www.ncbi.nlm.nih.gov/home/develop/api/)
- [Paper 모델](../../models/paper.py)
- [PaperCollectionService](../../services/paper_collection.py)
- [Topics 설정](../../config/topics.yaml)

---

## 📞 문제 해결

### CSV 저장 실패
- 디렉토리 권한 확인
- 경로 유효성 확인
- 디스크 공간 확인

### 로드 시 파일 미존재 오류
- 파일 경로 절대 경로 사용
- 파일 생성 후 다시 시도

### 토픽 정보 없음
- Paper/ScoredPaper 생성 시 topic 파라미터 전달 확인
- get_paper_details 호출 시 topic 파라미터 전달 확인

---

최종 업데이트: 2024년 1월
