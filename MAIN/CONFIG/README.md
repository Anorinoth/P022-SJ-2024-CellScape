# CONFIG Directory

This directory contains project configuration files and standards.

## Files

### CODING_STANDARDS.md
**Primary project coding standards and conventions.**

Key principle: **Fail-Fast Approach**
- No error suppression (no tryCatch, suppressWarnings, etc.)
- Operations fail immediately with full error messages
- Fix root causes instead of masking problems

### Usage

Before writing any code in this project:
1. Read `CODING_STANDARDS.md`
2. Follow the fail-fast principle
3. Let errors surface naturally

## Cursor AI Integration

This project is configured to enforce these standards through:

- **`.cursorrules`** (project root): Configures Cursor AI to follow these standards
- **`CODING_STANDARDS.md`** (this directory): Full documentation
- **`SpicyFlow.Rmd`** header: Quick reference reminder

Cursor AI will automatically:
- Avoid suggesting error suppression code
- Recommend fail-fast approaches
- Prioritize transparency over convenience

## Quick Reference

❌ **Never Use:**
- `tryCatch()`
- `try()`
- `suppressWarnings()`
- `suppressMessages()`

✅ **Always:**
- Let operations fail with full errors
- Read and fix error messages
- Document complex operations
- Use explicit error messages

---

See `CODING_STANDARDS.md` for complete details.

