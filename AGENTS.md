# Project Instructions (for AI Agents)

## Language Rule — Always Use ASD-STE100
- Always write all agent-facing text in ASD-STE100 Simplified Technical English.
- Use short sentences. Use approved words. Avoid ambiguity.
- Do not use idioms. Do not use long nouns. Do not use passive voice unless needed.
- Apply this rule to: comments, commit messages, pull request descriptions, documentation, and chat replies.


## Jules Learnings

Consolidated from the former memory files (`sentinel.md`, `bolt.md`, `palette.md`; directories removed). These are hard-won lessons from past sessions and must be respected in future changes.

### Security (`sentinel.md`)

1. **Sensitive Data Exposure in `GitHubConfig`** (2026-03-05)
   - Default dataclass `repr()` includes all fields, so secrets (e.g., GitHub PATs) leak into logs even when masked in settings files.
   - Always use `field(repr=False)` for any sensitive field in dataclasses.

2. **Missing Timeouts on External API Calls** (2026-03-29)
   - `requests` calls without a `timeout` can hang indefinitely, exhausting thread pools (DoS risk).
   - Always set an explicit timeout (e.g., `timeout=10`) on every external API call (`requests.get/post/patch`, etc.).

3. **SQL Injection Risk in Database Export** (2026-03-30)
   - Table names interpolated into SQL (where parameterization is impossible) must be validated against a strict allowlist, e.g. `["citations", "impact_factors", "paper_metadata"]`.

4. **CSV Injection (Formula Injection) Bypass** (2026-04-08)
   - Formula-trigger checks (`=`, `+`, `-`, `@`) can be bypassed by leading whitespace, which Excel still executes.
   - Always `.lstrip()` string values before checking the first character against formula triggers.

5. **XSS Vulnerability in `MdBookManager._format_page_content`** (2026-04-21)
   - Untrusted PubMed fields (`title`, `journal`, `abstract`, `authors`) were concatenated directly into HTML.
   - Always apply `html.escape()` (or equivalent context-aware escaping) to external text before embedding it in HTML.

6. **XSS Vulnerability in `MdBookManager.update_monthly_page`** (2026-04-24)
   - Same class of bug for table/link generation (`title`, `journal`, `pmid`, `doi` inside `<a href="...">`).
   - Same rule: escape all external strings embedded in HTML.

7. **Weak Hashing Algorithm in Cache Generation** (2026-05-11)
   - `hashlib.md5()` is flagged by Bandit (B324); use `hashlib.sha256()` instead, even for non-cryptographic caching.
   - Suppress false-positive SQL injection lints only after rigorous allowlist validation, using `# nosec B608`.

### Performance (`bolt.md`)

1. **Parallelize citation fetching** (2025-05-22)
   - Sequential IO-bound fetching of citations from external APIs is a major bottleneck; use `ThreadPoolExecutor`.
   - Always look for batch operations involving external API calls that can be parallelized.

2. **Thread-safe rate limiting** (2025-05-22)
   - `time.sleep()`-based rate limiting breaks under concurrency (threads pass the same check simultaneously).
   - Use a lock to reserve time slots per request so limits are strictly enforced across all threads.

### Accessibility (`palette.md`)

1. **Accessible Markdown Tables** (2025-05-15)
   - Standard markdown tables cannot carry ARIA labels or semantic hints.
   - Inject well-formed HTML tags (e.g., `<a aria-label="...">`, `<span title="...">`) directly into markdown table rows to improve screen-reader support without breaking table structure; prefer HTML injection for critical interactive elements in generated markdown content to ensure WCAG compliance.

---

## CSV Collections & Architecture (Consolidated from Copilot Contributions)

Commit `4d8d27c` ("update by copilot", 2026-01-19) implemented the CSV-based paper collection feature and management system:

### 1. Key Components & Implementation
- **Core Utility (`src/pubmed_miner/utils/csv_manager.py`)**:
  - `CSVManager.save_papers(papers, filepath, include_scoring=False, append=False)`: Save papers to CSV (creates parent directory automatically).
  - `CSVManager.append_papers(papers, filepath)`: Append new papers without duplicating headers.
  - `CSVManager.load_papers(filepath)`: Read CSV into list of dictionaries.
  - `CSVManager.update_collection(papers, filepath)`: Save or update with `ScoredPaper` objects.
- **Model Integrations**:
  - `Paper` & `ScoredPaper` (`src/pubmed_miner/models/paper.py`): Supports `topic` field for domain-specific categorization.
  - Integration with `PaperCollectionService` and automated pipeline (`scripts/automated_collection.py`).
- **Tests**: `tests/unit/test_csv_manager.py`
- **Examples**: `examples/csv_usage.py`, `examples/integration.py`

### 2. CSV Schema Reference
- **Basic Headers (`CSVManager.HEADERS`)**:
  `pmid,title,authors,journal,publication_date,doi,abstract,topic`
- **Scored Headers (`CSVManager.SCORED_HEADERS`)**:
  `pmid,title,authors,journal,publication_date,doi,abstract,topic,citation_count,impact_factor,score,rank`
- **Formatting Conventions**:
  - `authors`: Joined with semicolon (`;`), e.g., `"John Doe; Jane Smith"`.
  - `publication_date`: ISO 8601 formatted string (`YYYY-MM-DD` or `YYYY-MM-DDTHH:MM:SS`).
  - Encoding: Always UTF-8.

---

## Completion Status Markers [STRICT - EVERY TEXT RESPONSE]

At the end of **every** text-only response, include exactly one marker on its own line:

```
TASK_STATUS: COMPLETE
TASK_STATUS: BLOCKED
TASK_STATUS: INCOMPLETE
```

**Mapping to Turn Continuity rules:**
- `COMPLETE` → used when task is fully done and verified (aligns with Turn Continuity reason #1)
- `BLOCKED` → used when hard blocker or unrecoverable error (aligns with Turn Continuity reasons #2 and #3)
- `INCOMPLETE` → **platform-forced only** when OpenCode terminates mid-task due to limits; never choose deliberately

**Placement rules:**
- The marker must be the last line of the response
- No additional text after the marker
- Only ONE marker per response
