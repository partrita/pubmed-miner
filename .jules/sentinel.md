## 2026-03-05 - [Sensitive Data Exposure in GitHubConfig]
**Vulnerability:** The GitHub personal access token was being included in the default `repr()` of the `GitHubConfig` and `SystemConfig` dataclasses. This led to the secret being leaked in logs when these configuration objects were logged for debugging or during function call logging.
**Learning:** Default dataclass `repr()` includes all fields. Even when secrets are "masked" in settings files, once loaded into memory as dataclasses, they can be accidentally exposed through automated logging or error reporting that captures the object's state.
**Prevention:** Always use `field(repr=False)` for any sensitive data in dataclasses to ensure they are excluded from the default string representation.
