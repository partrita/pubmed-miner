# CSV Collections 저장 기능 구현 완료 (토픽 정보 추가됨)

## 📋 구현 내용

수집한 논문의 정보를 `data/collections.csv`에 저장하고, **토픽(topic) 정보**도 함께 관리하는 완전한 CSV 관리 시스템을 구현했습니다.

---

## 📁 추가된 파일들

### 1. 핵심 모듈
- **[src/pubmed_miner/utils/csv_manager.py](src/pubmed_miner/utils/csv_manager.py)**
  - `CSVManager` 클래스: CSV 파일 저장/로드/추가 기능
  - 기본 논문과 점수 정보 포함 논문 모두 지원
  - **토픽(topic) 필드 포함**
  - 자동 디렉토리 생성
  - 에러 핸들링

### 2. 테스트
- **[tests/unit/test_csv_manager.py](tests/unit/test_csv_manager.py)**
  - **14개의 포괄적인 테스트** (토픽 필드 테스트 포함)
  - 모든 테스트 통과 ✓

### 3. 문서 및 예제
- **[CSV_GUIDE.md](CSV_GUIDE.md)** - 상세 사용 가이드 (토픽 예제 포함)
- **[example_csv_usage.py](example_csv_usage.py)** - 기본 사용 예제 (토픽 포함)
- **[integration_example.py](integration_example.py)** - PubMed 통합 예제 (토픽 지원)

---

## 🚀 주요 기능

### 1. **토픽 정보와 함께 논문 저장**
```python
from src.pubmed_miner.utils import CSVManager
from src.pubmed_miner.models import Paper

papers = [
    Paper(
        pmid="12345678",
        title="Machine Learning in Healthcare",
        authors=["John Doe"],
        journal="Nature",
        publication_date=datetime(2024, 1, 15),
        topic="ai-healthcare",  # 토픽 정보
    )
]

CSVManager.save_papers(papers, "data/collections.csv")
```

### 2. **토픽과 함께 논문 수집**
```python
from src.pubmed_miner.services.paper_collection import PaperCollectionService

service = PaperCollectionService(email="your.email@example.com")
pmids = service.search_papers("machine learning AND healthcare", max_results=100)

# 토픽 정보와 함께 논문 상세 정보 조회
papers = service.get_paper_details(pmids, topic="ai-healthcare")

CSVManager.save_papers(papers, "data/collections.csv")
```

### 3. **여러 토픽으로 논문 수집**
```python
topics = {
    "ai-drug-discovery": "machine learning AND drug discovery",
    "cancer-immunotherapy": "cancer AND immunotherapy",
}

for topic_name, query in topics.items():
    pmids = service.search_papers(query, max_results=50)
    papers = service.get_paper_details(pmids, topic=topic_name)
    # 기존 파일에 추가
    CSVManager.append_papers(papers, "data/collections.csv")
```

### 4. **토픽별 논문 분석**
```python
papers = CSVManager.load_papers("data/collections.csv")

# 토픽별 논문 수 집계
from collections import Counter

topics = [p["topic"] for p in papers if p["topic"]]
topic_counts = Counter(topics)

print("📊 토픽별 논문 수:")
for topic, count in topic_counts.most_common():
    print(f"  {topic}: {count}개")
```

---

## 📊 CSV 파일 구조

### 기본 컬럼 (Paper 객체)
```
pmid | title | authors | journal | publication_date | doi | abstract | topic
```

### 확장 컬럼 (ScoredPaper 객체)
```
pmid | title | authors | journal | publication_date | doi | abstract | topic | 
citation_count | impact_factor | score | rank
```

**예제 데이터:**
```csv
pmid,title,authors,journal,publication_date,doi,abstract,topic,citation_count,impact_factor,score,rank
11111111,"Highly Cited Research","Dr. Einstein","Science","2023-06-01T00:00:00","10.1234/score.2023.001","",cancer-immunotherapy,150,42.5,95.5,1
22222222,"Important Study","Prof. Newton","Nature","2023-07-15T00:00:00","10.1234/score.2023.002","",de-novo-protein-design,87,39.8,78.3,2
99999999,"Additional Paper","New Author","Journal of Examples","2024-01-19T00:00:00","10.1234/additional.2024.001","",bioinformatics,,,,
```

---

## ✅ 테스트 결과

```
tests/unit/test_csv_manager.py::TestCSVManager
✓ test_save_papers_creates_file
✓ test_save_papers_with_correct_headers
✓ test_save_papers_contains_data
✓ test_save_scored_papers_with_scoring_info
✓ test_append_papers
✓ test_load_papers
✓ test_save_papers_empty_list_raises_error
✓ test_append_papers_empty_list_raises_error
✓ test_load_nonexistent_file_raises_error
✓ test_authors_joined_with_semicolon
✓ test_publication_date_as_isoformat
✓ test_update_collection
✓ test_topic_field_included
✓ test_topic_value_saved_correctly

======================== 14 passed in 0.08s ========================
```

---

## 🔄 완전한 워크플로우 예제

```python
from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.utils import CSVManager

# 1. PubMed 검색
service = PaperCollectionService(email="your.email@example.com")
pmids = service.search_papers("machine learning AND healthcare", max_results=100)

# 2. 논문 정보 조회
papers = service.get_paper_details(pmids)

# 3. CSV 저장
CSVManager.save_papers(papers, "data/collections.csv")

print(f"✓ {len(papers)}개의 논문을 저장했습니다.")
```

---

## 💡 사용 시나리오

### 시나리오 1: 초기 수집 및 저장
```python
# 새로운 CSV 파일 생성
CSVManager.save_papers(papers, "data/collections.csv")
```

### 시나리오 2: 정기적인 추가 수집
```python
# 기존 CSV에 새 논문 추가
CSVManager.append_papers(new_papers, "data/collections.csv")
```

### 시나리오 3: 점수 계산 후 업데이트
```python
# 스코어가 포함된 논문으로 업데이트
CSVManager.update_collection(scored_papers, "data/collections.csv")
```

### 시나리오 4: 데이터 분석
```python
# CSV에서 데이터 로드하여 분석
papers = CSVManager.load_papers("data/collections.csv")
# 추가 분석 수행...
```

---

## 🛠️ API 참조

### `CSVManager.save_papers()`
- 새 CSV 파일 생성 또는 덮어쓰기
- 점수 정보 선택적 포함 가능

### `CSVManager.append_papers()`
- 기존 CSV에 데이터 추가
- 헤더 자동 관리

### `CSVManager.load_papers()`
- CSV 파일에서 데이터 로드
- 딕셔너리 리스트 반환

### `CSVManager.update_collection()`
- ScoredPaper 정보 저장
- 점수/순위 정보 포함

---

## 📦 의존성

### 기존 의존성 활용
- `csv` (Python 표준 라이브러리)
- `pathlib` (Python 표준 라이브러리)
- `datetime` (Python 표준 라이브러리)
- `Paper`, `ScoredPaper` 모델

### 추가 요구사항
- 없음 (모두 표준 라이브러리 사용)

---

## 📋 파일 구조

```
/workspaces/pubmed-miner/
├── src/pubmed_miner/
│   ├── utils/
│   │   ├── __init__.py (업데이트됨)
│   │   └── csv_manager.py ✨ (새로 추가)
│   ├── models/
│   │   └── paper.py (기존)
│   └── services/
│       └── paper_collection.py (기존)
├── tests/unit/
│   └── test_csv_manager.py ✨ (새로 추가)
├── data/
│   └── collections.csv ✨ (생성됨)
├── example_csv_usage.py ✨ (새로 추가)
├── integration_example.py ✨ (새로 추가)
├── CSV_GUIDE.md ✨ (새로 추가)
└── IMPLEMENTATION_SUMMARY.md ✨ (이 파일)
```

---

## ⚡ 빠른 시작

1. **기본 사용**
```bash
python example_csv_usage.py
```

2. **통합 예제 (PubMed 필요)**
```bash
# Entrez 이메일 설정 필수
python integration_example.py
```

3. **테스트**
```bash
pytest tests/unit/test_csv_manager.py -v
```

---

## 🎯 특징

✨ **사용 용이성**
- 직관적인 API
- 자동 디렉토리 생성
- 명확한 에러 메시지

🔒 **안정성**
- 입력값 검증
- 예외 처리
- 파일 무결성 보장

📊 **유연성**
- 기본 논문 및 점수 정보 포함 논문 지원
- 저장/로드/추가 기능
- 헤더 자동 관리

📝 **확장성**
- 추가 컬럼 확장 가능
- 커스텀 필터링 가능

---

## 📚 문서

자세한 사용 방법은 [CSV_GUIDE.md](CSV_GUIDE.md)를 참고하세요.

주요 내용:
- 기본 사용법
- 완전한 워크플로우
- API 참조
- 실용적인 예제들
- 문제 해결

---

## 🔔 다음 단계

### 선택적 개선 사항
1. 대용량 파일을 위한 배치 처리 최적화
2. CSV 필터링 및 검색 기능
3. 데이터 검증 및 정제 기능
4. 다른 포맷 지원 (JSON, Excel 등)
5. 데이터베이스 백업 기능

---

최종 업데이트: 2024년 1월 19일
구현 상태: ✅ 완료 및 테스트됨
