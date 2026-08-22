# AGENTS.md

Guidelines and consolidated history for AI coding agents working on this repository.

This repository has been developed with the help of multiple AI agents:

- **Jules** (Google): security, performance, and accessibility improvements — session learnings recorded in `.Jules/`
- **GitHub Copilot**: implementation of the CSV paper-collection feature (commit `4d8d27c`, "update by copilot")

---

## Jules Learnings

Consolidated from the `.Jules/` memory files (`sentinel.md`, `bolt.md`, `palette.md`). These are hard-won lessons from past sessions and must be respected in future changes.

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

## GitHub Copilot Contributions

Commit `4d8d27c` ("update by copilot", 2026-01-19) implemented the CSV-based paper collection feature:

- Added reference documentation: `CSV_GUIDE.md`, `IMPLEMENTATION_SUMMARY.md`
- Added core module: `src/pubmed_miner/utils/csv_manager.py`
- Added tests: `tests/unit/test_csv_manager.py`
- Integrated CSV handling into:
  - `src/pubmed_miner/services/paper_collection.py`
  - `src/pubmed_miner/models/paper.py`
  - `src/pubmed_miner/utils/__init__.py`
- Updated configuration/data: `config/topics.yaml`, `data/collections.csv`
- Usage examples (now under `examples/`): `examples/csv_usage.py`, `examples/integration.py`

When modifying anything related to CSV collection/export, consult `CSV_GUIDE.md` and keep `IMPLEMENTATION_SUMMARY.md` up to date.

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
