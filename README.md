# PubMed Miner - 필수 논문 추천 도구

`pubmed-miner`는 PubMed에서 필수 연구 논문을 자동으로 발견, 분석, 추천하는 지능형 Python 도구입니다. 고급 논문 수집, 인용 분석, 자동화된 GitHub 통합을 결합하여 연구자들이 해당 분야의 가장 중요한 출판물을 최신 상태로 유지할 수 있도록 도와줍니다.

## 🚀 주요 기능

### 핵심 기능
- 자동화된 PubMed 검색: 구성 가능한 검색 쿼리로 BioPython의 Entrez API를 사용하여 논문 검색
- 지능형 논문 점수 매기기: 인용, 저널 임팩트 팩터, 최신성, 관련성을 고려한 다중 요소 점수 알고리즘
- 필수 논문 선택: 각 연구 주제에 대해 가장 중요한 논문을 자동으로 식별
- GitHub 통합: 큐레이션된 논문 목록으로 GitHub 이슈 생성 및 관리
- 변경 추적: 새로운 논문과 순위 변경을 시간에 따라 모니터링

### 고급 기능
- 인용 분석: 여러 소스(PMC, Crossref)에서 인용 수 수집
- 저널 임팩트 팩터: 더 나은 논문 평가를 위한 저널 품질 지표 통합
- 배치 처리: 속도 제한을 통한 대용량 논문 데이터셋의 효율적 처리
- 캐싱 시스템: API 호출을 최소화하고 성능을 향상시키는 스마트 캐싱
- 오류 처리: 강력한 오류 복구 및 재시도 메커니즘

### 자동화 기능
- GitHub Actions 통합: 완전 자동화된 일일 논문 수집
- 구성 관리: YAML 기반 주제 및 시스템 구성
- 로깅 및 모니터링: 문제 해결 및 모니터링을 위한 포괄적 로깅
- 유연한 스케줄링: 사용자 정의 가능한 자동화 일정

## 🛠️ 기술 스택

### 핵심 의존성
- BioPython - PubMed/NCBI API 통합
- PyGithub - 이슈 관리를 위한 GitHub API 통합
- PyYAML - 구성 파일 관리
- Requests - 외부 API용 HTTP 클라이언트
- Pandas - 데이터 조작 및 분석
- Python-dateutil - 고급 날짜/시간 처리

### 외부 API 및 서비스
- PubMed/Entrez API - 논문 메타데이터 및 검색
- PMC LinkOut - 인용 수 검색
- Crossref API - 대체 인용 소스
- GitHub API - 이슈 생성 및 관리
- Journal Impact Factor 데이터베이스 - 저널 품질 지표

### 개발 도구
- Pytest - 커버리지를 포함한 테스트 프레임워크
- Black - 코드 포맷팅
- MyPy - 정적 타입 검사
- Pre-commit - 코드 품질을 위한 Git 훅

## 📦 설치

### 사전 요구사항
- Python 3.8 이상
- Git
- GitHub 계정 (자동화 기능용)

### 빠른 설치

1. 저장소 복제:
   ```bash
   git clone https://github.com/yourusername/pubmed-miner.git
   cd pubmed-miner
   ```

2. uv로 설치:
   ```bash
   uv sync
   ```

## 🚀 빠른 시작

### 1. 기본 설정

먼저 구성 파일을 설정합니다:

```bash
uv run python setup_automation.py
```

이 명령은 `config/` 디렉토리에 필요한 구성 파일을 생성합니다.

### 2. 연구 주제 구성

`config/topics.yaml` 편집:
```yaml
topics:
  - name: "machine-learning-healthcare"
    query: "machine learning AND healthcare"
    max_papers: 1000
    essential_count: 15
    enabled: true
  
  - name: "cancer-immunotherapy"
    query: "cancer AND immunotherapy"
    max_papers: 800
    essential_count: 20
    enabled: true
```

### 3. GitHub 통합 구성

`config/settings.yaml` 편집:
```yaml
github:
  repository: "your-username/your-repo"
  issue_labels:
    - "essential-papers"
    - "automated"

scoring_weights:
  citation_weight: 0.4
  impact_factor_weight: 0.3
  recency_weight: 0.2
  relevance_weight: 0.1
```

### 4. 환경 변수 설정

PubMed 이메일은 다음 두 가지 방법으로 설정할 수 있습니다:

**방법 1: 환경 변수 사용**
```bash
export GITHUB_TOKEN="your_github_token_here"
export PUBMED_EMAIL="your_email@example.com"  # 선택사항이지만 권장
```

**방법 2: 설정 파일 사용**
`config/settings.yaml` 파일에서 직접 설정:
```yaml
pubmed:
  email: "your_email@example.com"  # 여기에 실제 이메일 주소 입력
```

> **참고**: 환경 변수가 설정 파일보다 우선순위가 높습니다.

### 5. 수동 수집 실행

시스템을 수동으로 테스트:
```bash
uv run python automated_collection.py
```

이 명령은 다음을 수행합니다:
1. 주제와 일치하는 논문을 PubMed에서 검색
2. 인용 및 저널 임팩트 데이터 수집
3. 중요도에 따라 논문 점수 매기기 및 순위 지정
4. 필수 논문 목록으로 GitHub 이슈 생성/업데이트

## 🤖 GitHub Actions 자동화

### 자동화된 수집 설정

1. 이 저장소를 GitHub 계정으로 포크/복제

2. GitHub Secrets 설정:
   - 저장소 Settings → Secrets and variables → Actions로 이동
   - `GITHUB_TOKEN` 추가 (GitHub에서 자동 제공)
   - 선택적으로 이메일과 함께 `PUBMED_EMAIL` 추가

3. 워크플로우 구성:
   워크플로우 파일 `.github/workflows/collect-papers.yml`은 이미 다음과 같이 실행되도록 구성되어 있습니다:
   - 매일 UTC 오전 6시
   - main 브랜치로 푸시할 때 (테스트용)
   - Actions 탭에서 수동 트리거

4. 일정 사용자 정의 (선택사항):
   ```yaml
   schedule:
     - cron: '0 6 * * *'  # 매일 UTC 오전 6시
     - cron: '0 18 * * 0' # 매주 일요일 오후 6시 UTC
   ```

### 자동화가 수행하는 작업

1. 구성된 주제와 일치하는 논문을 PubMed에서 검색
2. 인용, 저널 임팩트 팩터를 포함한 메타데이터 수집
3. 가중 알고리즘을 사용하여 논문 점수 매기기
4. 필수 논문 선택 (각 주제별 상위 N개)
5. 형식화된 논문 목록으로 GitHub 이슈 생성/업데이트
6. 변경 사항 추적 및 새 논문이나 순위 변경에 대한 댓글 추가

## 📊 점수 알고리즘 작동 방식

시스템은 정교한 다중 요소 점수 알고리즘을 사용합니다:

### 점수 요소

1. 인용 수 (기본 가중치 40%)
   - PMC LinkOut 또는 Crossref의 원시 인용 수
   - 출판 연령으로 정규화

2. 저널 임팩트 팩터 (기본 가중치 30%)
   - Journal Citation Reports (JCR) 데이터
   - JCR 데이터가 없는 저널의 추정 임팩트 팩터

3. 최신성 (기본 가중치 20%)
   - 지수 감쇠를 적용한 출판 날짜
   - 최근 발견과 기존 지식의 균형

4. 관련성 (기본 가중치 10%)
   - 제목과 초록 기반 쿼리 일치 점수
   - 키워드 밀도 및 의미적 매칭

### 점수 공식

```
최종 점수 = (인용 점수 × 0.4) + 
           (임팩트 팩터 점수 × 0.3) + 
           (최신성 점수 × 0.2) + 
           (관련성 점수 × 0.1)
```

모든 점수는 가중치 적용 전에 0-1 범위로 정규화됩니다.

## 🧪 테스트

### 로컬 테스트 (GitHub 토큰 불필요)

GitHub 자격 증명 없이 로컬 개발 및 테스트를 위해:

```bash
# 로컬 테스트 스크립트 실행
uv run python test_local_without_token.py
```

이 스크립트는 다음을 보여줍니다:
- 모의 모드 기능
- 논문 점수 매기기 및 순위 지정
- 이슈 형식화
- GitHub 액세스 없이 모든 핵심 기능

### 모든 테스트 실행
```bash
uv run pytest
```

### 특정 테스트 카테고리 실행
```bash
# 단위 테스트만
uv run pytest tests/unit/

# 통합 테스트만
uv run pytest tests/integration/

# 커버리지 리포트와 함께
uv run pytest --cov=src/pubmed_miner --cov-report=html
```

### 테스트 구성
```bash
# 자세한 출력으로 테스트 실행
uv run pytest -v

# 테스트 실행 후 첫 번째 실패에서 중지
uv run pytest -x

# 빠른 테스트만 실행 (느린 통합 테스트 건너뛰기)
uv run pytest -m "not slow"
```

### 개발용 모의 모드

시스템은 GitHub 토큰이 없을 때 자동으로 감지하여 모의 모드로 전환합니다:

```python
# 코드나 테스트에서
from pubmed_miner.models import GitHubConfig
from pubmed_miner.services.github_manager import GitHubIssuesManager

# 자동으로 모의 모드를 사용합니다
config = GitHubConfig(
    token="mock_token_for_local_testing",  # 또는 빈 문자열
    repository="test/repo",
    issue_labels=["test"]
)

manager = GitHubIssuesManager(config)
# manager.mock_mode는 True가 됩니다
```

## 📁 프로젝트 구조

```
pubmed-miner/
├── book_src/                   # MdBook 소스 (문서)
├── src/pubmed_miner/           # 메인 패키지
│   ├── models/                 # 데이터 모델 및 구성
│   ├── services/               # 핵심 비즈니스 로직 서비스
│   ├── scoring/                # 논문 점수 알고리즘
│   ├── utils/                  # 유틸리티 함수 및 헬퍼
│   └── data/                   # 정적 데이터 및 데이터베이스
├── tests/                      # 테스트 스위트
│   ├── unit/                   # 단위 테스트
│   └── integration/            # 통합 테스트
├── config/                     # 구성 파일
│   ├── topics.yaml            # 연구 주제 구성
│   └── settings.yaml          # 시스템 설정
├── .github/workflows/          # GitHub Actions 워크플로우
├── automated_collection.py     # 메인 자동화 스크립트
├── setup_automation.py        # 설정 및 검증 스크립트
└── pyproject.toml             # 프로젝트 구성 및 의존성
```
```
## 🔧 구성 참조

### 주제 구성 (`config/topics.yaml`)

```yaml
topics:
  - name: "topic-identifier"           # 주제의 고유 식별자
    query: "search terms AND keywords" # PubMed 검색 쿼리
    max_papers: 1000                   # 검색할 최대 논문 수
    essential_count: 15                # 선택할 필수 논문 수
    enabled: true                      # 이 주제 활성화/비활성화
    
  # 고급 구성
  - name: "advanced-topic"
    query: "complex query[MeSH] AND (term1 OR term2)"
    max_papers: 500
    essential_count: 10
    enabled: true
    filters:                           # 선택적 필터
      min_citation_count: 5
      min_impact_factor: 2.0
      date_range: "2020:2024"
```

### 시스템 설정 (`config/settings.yaml`)

```yaml
# GitHub 통합
github:
  repository: "username/repository-name"
  issue_labels:
    - "essential-papers"
    - "automated"
    - "research"

# 점수 알고리즘 가중치
scoring_weights:
  citation_weight: 0.4      # 인용 수 중요도 (0.0-1.0)
  impact_factor_weight: 0.3 # 저널 임팩트 팩터 중요도
  recency_weight: 0.2       # 출판 최신성 중요도
  relevance_weight: 0.1     # 쿼리 관련성 중요도

# API 속도 제한
rate_limits:
  pubmed_delay: 0.5         # PubMed API 호출 간 초
  citation_delay: 1.0       # 인용 API 호출 간 초
  github_delay: 0.5         # GitHub API 호출 간 초

# 캐싱 설정
cache_settings:
  citation_cache_days: 30   # 인용 데이터 캐시 일수
  impact_factor_cache_days: 365  # 임팩트 팩터 데이터 캐시 일수
  paper_metadata_cache_days: 7   # 논문 메타데이터 캐시 일수

# 로깅 구성
logging:
  level: "INFO"             # DEBUG, INFO, WARNING, ERROR
  file: "logs/pubmed_miner.log"
  max_file_size: "10MB"
  backup_count: 5
```

## 🚨 문제 해결

### 일반적인 문제 및 해결책

#### 1. GitHub 토큰 문제

문제: `GitHub token not found` 또는 `Authentication failed`

해결책:
- `GITHUB_TOKEN` 환경 변수가 설정되었는지 확인
- 토큰 권한 확인 (`repo` 및 `issues` 범위 필요)
- 토큰이 만료되지 않았는지 확인
- GitHub Actions의 경우 시크릿이 올바르게 구성되었는지 확인

```bash
# 로컬에서 토큰 테스트
export GITHUB_TOKEN="your_token_here"
uv run python -c "
import os
from github import Github
g = Github(os.environ['GITHUB_TOKEN'])
print(f'Authenticated as: {g.get_user().login}')
"
```

#### 2. PubMed API 문제

문제: `Rate limit exceeded` 또는 `API timeout`

해결책:
- `config/settings.yaml`에서 지연 시간 증가
- `PUBMED_EMAIL` 환경 변수 설정
- 인터넷 연결 확인
- PubMed 서비스 상태 확인

```yaml
# config/settings.yaml에서
rate_limits:
  pubmed_delay: 1.0  # 기본값 0.5에서 증가
```

#### 3. 논문을 찾을 수 없음

문제: `No papers found for topic` 또는 빈 결과

해결책:
- [PubMed Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/)를 사용하여 PubMed 쿼리 구문 확인
- `config/topics.yaml`에서 주제 구성 확인
- 주제가 활성화되었는지 확인 (`enabled: true`)
- 더 넓은 검색어 시도

#### 4. 인용 데이터 누락

문제: 논문의 인용 수가 0이거나 임팩트 팩터가 누락됨

해결책:
- 매우 최근 논문의 경우 정상적인 현상
- 논문이 인용 데이터베이스에 색인되었는지 확인
- 저널 이름이 올바르게 매칭되었는지 확인
- 인용 중요도를 줄이기 위해 점수 가중치 조정 고려

#### 5. GitHub 이슈가 생성되지 않음

문제: 자동화가 실행되지만 이슈가 나타나지 않음

해결책:
- 저장소 권한 확인 (이슈가 활성화되어야 함)
- `config/settings.yaml`에서 저장소 이름 확인
- 오류에 대한 GitHub Actions 로그 확인
- 필수 논문이 발견되었는지 확인 (로그 확인)

### 디버그 모드

자세한 문제 해결을 위해 디버그 로깅 활성화:

```yaml
# config/settings.yaml에서
logging:
  level: "DEBUG"
```

또는 디버그 플래그로 실행:
```bash
uv run python automated_collection.py --debug
```

### 로그 파일

자세한 오류 정보는 다음 로그 파일을 확인하세요:
- `logs/pubmed_miner.log` - 메인 애플리케이션 로그
- GitHub Actions 로그 - 저장소의 Actions 탭에서 확인 가능

## ❓ 자주 묻는 질문

### 일반 사용법

Q: 자동화를 얼마나 자주 실행해야 하나요?
A: 대부분의 연구 분야에서는 일일 실행이 잘 작동합니다. 빠르게 발전하는 분야의 경우 하루 두 번, 느린 분야의 경우 주간 실행으로도 충분할 수 있습니다.

Q: 비의학 연구에도 사용할 수 있나요?
A: 시스템은 PubMed/의학 연구에 최적화되어 있지만, PubMed에서 사용 가능한 다른 도메인에 대해 쿼리를 조정할 수 있습니다.

Q: 몇 개의 주제를 구성할 수 있나요?
A: 엄격한 제한은 없지만 API 속도 제한을 고려하세요. 5-10개 주제로 시작하여 성능에 따라 확장하세요.

### 기술적 질문

Q: 점수 알고리즘은 얼마나 정확한가요?
A: 알고리즘은 여러 품질 지표를 결합하며 가중치를 통해 조정 가능합니다. 높은 영향력의 논문을 표면화하도록 설계되었지만 전문가 판단을 보완하는 것이지 대체하는 것은 아닙니다.

Q: GitHub 이슈 형식을 사용자 정의할 수 있나요?
A: 네, `src/pubmed_miner/services/github_manager.py`의 `GitHubIssuesManager` 클래스를 수정하세요.

Q: 구성을 어떻게 백업하나요?
A: `config/` 디렉토리에 모든 설정이 포함되어 있습니다. 버전 관리에 커밋하거나 정기적으로 백업하세요.

Q: 여러 저장소에서 실행할 수 있나요?
A: 네, 별도의 인스턴스를 실행하거나 구성에서 여러 저장소를 지원하도록 코드를 수정할 수 있습니다.

### 성능 질문

Q: 일반적인 실행 시간은 얼마나 걸리나요?
A: 주제와 논문 수에 따라 다릅니다. API 지연을 포함하여 1000개 논문이 있는 주제당 2-5분을 예상하세요.

Q: 프로세스를 어떻게 가속화할 수 있나요?
A: 주제당 `max_papers`를 줄이거나, 캐시 보존 기간을 늘리거나, 주제를 병렬로 실행하세요 (코드 수정 필요).

Q: API 속도 제한은 무엇인가요?
A: PubMed: ~초당 3개 요청, GitHub: 시간당 5000개 요청. 시스템에는 내장된 속도 제한이 포함되어 있습니다.

## 🤝 기여하기

기여를 환영합니다! 시작하는 방법은 다음과 같습니다:

### 개발 환경 설정

1. 저장소를 포크하고 복제
2. 개발 환경 설정:
   ```bash
   uv sync --extra dev --extra test
   uv run pre-commit install
   ```
3. 모든 것이 작동하는지 테스트 실행:
   ```bash
   uv run pytest
   ```

### GitHub 토큰 없이 로컬 테스트

시스템은 이제 GitHub 토큰 없이 로컬 테스트를 지원합니다. `GITHUB_TOKEN`이 제공되지 않을 때:

- GitHub 기능이 **모의 모드**로 실행됩니다
- 이슈가 생성되는 대신 시뮬레이션되고 로그됩니다
- 다른 모든 기능 (PubMed 검색, 점수 매기기 등)은 정상적으로 작동합니다

로컬에서 테스트하려면:

```bash
# GitHub 토큰 없이 실행 (모의 모드 사용)
uv run python automated_collection.py

# 설정 검증 실행
uv run python setup_automation.py

# 특정 테스트 실행
uv run pytest tests/unit/
uv run pytest tests/integration/ --skip-external
```

모의 모드는 GitHub에서 수행되었을 작업을 로그합니다:
```
[MOCK MODE] Would create/update issue for topic: machine-learning
[MOCK MODE] Issue would contain 15 papers
[MOCK MODE] Paper 1: Deep Learning in Medical Imaging (Score: 95.2)
```

### 코드 스타일

코드 품질을 유지하기 위해 여러 도구를 사용합니다:
- Black - 코드 포맷팅
- isort - 임포트 정렬
- mypy - 타입 검사
- pytest - 테스트

모든 검사 실행:
```bash
# 코드 포맷팅
uv run black src/ tests/

# 임포트 정렬
uv run isort src/ tests/

# 타입 검사
uv run mypy src/

# 테스트 실행
uv run pytest
```

### 변경 사항 제출

1. 기능 브랜치 생성: `git checkout -b feature-name`
2. 변경 사항을 만들고 테스트 추가
3. 모든 검사가 통과하는지 확인: `uv run pytest && uv run black --check src/ && uv run mypy src/`
4. 명확한 메시지로 커밋
5. 푸시하고 풀 리퀘스트 생성

### 기여 영역

- 새로운 점수 알고리즘 - 대체 논문 순위 방법
- 추가 데이터 소스 - 더 많은 인용 또는 임팩트 팩터 소스
- UI 개선 - 더 나은 GitHub 이슈 포맷팅
- 성능 최적화 - 더 빠른 데이터 수집 및 처리
- 문서화 - 예제, 튜토리얼, 가이드

## 📄 라이선스

이 프로젝트는 MIT 라이선스 하에 라이선스됩니다 - 자세한 내용은 [LICENSE](LICENSE) 파일을 참조하세요.

## 🙏 감사의 말

### API 및 서비스
- NCBI/PubMed - 의학 문헌에 대한 무료 액세스 제공
- GitHub - 호스팅 및 자동화 인프라 제공
- Crossref - 인용 데이터 제공
- Journal Citation Reports - 임팩트 팩터 데이터 제공

### 오픈 소스 라이브러리
- BioPython - NCBI API 통합
- PyGithub - GitHub API 클라이언트
- Requests - HTTP 클라이언트 라이브러리
- Pandas - 데이터 조작
- PyYAML - 구성 관리

### 인용

연구에서 이 도구를 사용하는 경우 다음과 같이 인용해 주세요:

```bibtex
@software{pubmed_miner,
  title={PubMed Miner: Essential Papers Recommender},
  author={Kim, Taeyoon},
  year={2024},
  url={https://github.com/example/pubmed-miner}
}
```

---

연구 커뮤니티를 위해 ❤️로 제작되었습니다

질문, 문제 또는 제안 사항이 있으시면 GitHub에서 [이슈를 열어주세요](https://github.com/example/pubmed-miner/issues).