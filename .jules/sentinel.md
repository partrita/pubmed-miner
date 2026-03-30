## 2026-03-05 - [Sensitive Data Exposure in GitHubConfig]
**Vulnerability:** The GitHub personal access token was being included in the default `repr()` of the `GitHubConfig` and `SystemConfig` dataclasses. This led to the secret being leaked in logs when these configuration objects were logged for debugging or during function call logging.
**Learning:** Default dataclass `repr()` includes all fields. Even when secrets are "masked" in settings files, once loaded into memory as dataclasses, they can be accidentally exposed through automated logging or error reporting that captures the object's state.
**Prevention:** Always use `field(repr=False)` for any sensitive data in dataclasses to ensure they are excluded from the default string representation.

## 2026-03-29 - [Missing Timeouts on External API Calls]
**Vulnerability:** External HTTP requests via the `requests` library lacked a `timeout` parameter in `GitHubIssuesManager`. This creates a Denial of Service (DoS) risk where threads can hang indefinitely if the API server is slow or non-responsive.
**Learning:** The `requests` library in Python does not enforce a timeout by default. Indefinite hangs can lead to thread pool exhaustion and application crashes.
**Prevention:** Always explicitly set a `timeout` parameter (e.g., `timeout=10`) for every external API call made using `requests.get`, `requests.post`, `requests.patch`, etc.
