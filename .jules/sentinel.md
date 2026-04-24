## 2026-03-05 - [Sensitive Data Exposure in GitHubConfig]
**Vulnerability:** The GitHub personal access token was being included in the default `repr()` of the `GitHubConfig` and `SystemConfig` dataclasses. This led to the secret being leaked in logs when these configuration objects were logged for debugging or during function call logging.
**Learning:** Default dataclass `repr()` includes all fields. Even when secrets are "masked" in settings files, once loaded into memory as dataclasses, they can be accidentally exposed through automated logging or error reporting that captures the object's state.
**Prevention:** Always use `field(repr=False)` for any sensitive data in dataclasses to ensure they are excluded from the default string representation.

## 2026-03-29 - [Missing Timeouts on External API Calls]
**Vulnerability:** External HTTP requests via the `requests` library lacked a `timeout` parameter in `GitHubIssuesManager`. This creates a Denial of Service (DoS) risk where threads can hang indefinitely if the API server is slow or non-responsive.
**Learning:** The `requests` library in Python does not enforce a timeout by default. Indefinite hangs can lead to thread pool exhaustion and application crashes.
**Prevention:** Always explicitly set a `timeout` parameter (e.g., `timeout=10`) for every external API call made using `requests.get`, `requests.post`, `requests.patch`, etc.

## 2026-03-30 - [SQL Injection Risk in Database Export]
**Vulnerability:** The `CacheManager.export_cache_data` function allowed SQL injection by directly formatting user-provided table names (`table_name`) into a raw SQL query string (`f"SELECT * FROM {table}"`) without validation.
**Learning:** Even utility or internal-facing functions can be exposed to untrusted input. When parameterized queries cannot be used for table names, strict validation against a predefined allowlist is necessary.
**Prevention:** Validate input strings against an explicit list of allowed names (e.g., `["citations", "impact_factors", "paper_metadata"]`) before interpolating them into SQL queries.

## 2026-04-08 - [CSV Injection (Formula Injection) Bypass]
**Vulnerability:** The existing CSV injection protection checked if the first character of a string field was a formula trigger (`=`, `+`, `-`, `@`). However, this check could be bypassed if the malicious payload started with whitespace characters (e.g., space or tab) because Excel still executes formulas with leading spaces.
**Learning:** Basic character matching for CSV injection is insufficient if it doesn't account for how spreadsheet applications parse the input. Leading whitespace must be stripped before evaluating whether a payload starts with an executable formula character.
**Prevention:** Always use `.lstrip()` on string values before checking if the first character is a known trigger character (e.g., `stripped_value and stripped_value[0] in ('=', '+', '-', '@')`).
## 2026-04-21 - [XSS Vulnerability in MdBookManager]
**Vulnerability:** The `MdBookManager._format_page_content` method concatenated external, potentially untrusted strings from PubMed (such as `paper.title`, `paper.journal`, `paper.abstract` and `paper.authors`) directly into HTML markup (`<div class="paper-title">...</div>`) without escaping them. This exposed the application to Cross-Site Scripting (XSS).
**Learning:** Any dynamic generation of HTML or Markdown that embeds external string data must consider XSS attacks, even if the data originates from a "trusted" external API like PubMed, because the content can still contain HTML-like syntax or malicious payloads.
**Prevention:** Always use `html.escape()` or an equivalent context-aware escaping mechanism on external text fields before embedding them into HTML content.

## 2026-04-24 - [XSS Vulnerability in MdBookManager table generation]
**Vulnerability:** The `MdBookManager.update_monthly_page` method concatenated external, potentially untrusted strings from PubMed (such as `paper.title`, `paper.journal`, `paper.pmid` and `paper.doi`) directly into HTML markup (`<a href="...">...</a>`) without escaping them. This exposed the application to Cross-Site Scripting (XSS).
**Learning:** Any dynamic generation of HTML or Markdown that embeds external string data must consider XSS attacks, even if the data originates from a "trusted" external API like PubMed, because the content can still contain HTML-like syntax or malicious payloads.
**Prevention:** Always use `html.escape()` or an equivalent context-aware escaping mechanism on external text fields before embedding them into HTML content.
