## 2025-05-15 - [Accessible Markdown Tables]
**Learning:** When generating markdown tables for documentation (like mdBook), standard markdown syntax often lacks the ability to add ARIA labels or semantic hints. Injecting well-formed HTML tags (e.g., `<a>` with `aria-label`, `<span>` with `title`) directly into the markdown table rows significantly improves screen reader support and provides contextual information without breaking the table structure.
**Action:** Always prefer HTML injection for critical interactive elements in generated markdown content to ensure WCAG compliance.
