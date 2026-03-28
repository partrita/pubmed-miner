## 2026-03-05 - [Sensitive Data Exposure in GitHubConfig]
**Vulnerability:** The GitHub personal access token was being included in the default `repr()` of the `GitHubConfig` and `SystemConfig` dataclasses. This led to the secret being leaked in logs when these configuration objects were logged for debugging or during function call logging.
**Learning:** Default dataclass `repr()` includes all fields. Even when secrets are "masked" in settings files, once loaded into memory as dataclasses, they can be accidentally exposed through automated logging or error reporting that captures the object's state.
**Prevention:** Always use `field(repr=False)` for any sensitive data in dataclasses to ensure they are excluded from the default string representation.
## 2024-05-24 - Missing Timeout Parameter in External HTTP Requests
**Vulnerability:** External HTTP requests via the `requests` library were missing a `timeout` parameter, allowing potential Denial of Service (DoS) attacks or process hangs due to resource exhaustion when interacting with external APIs (like GitHub's).
**Learning:** The default behavior of `requests` is to wait indefinitely for a response if no timeout is explicitly specified. This behavior can cause indefinite hangs and exhaust thread pools or processes if the downstream service is unresponsive.
**Prevention:** Always explicitly set a `timeout` parameter (e.g., `timeout=10`) for all external HTTP calls in `requests` to ensure the application fails fast and gracefully handles slow network responses.
