# Project Coding Standards - SJ_04.2

## Error Handling Philosophy

### ❌ DO NOT USE: tryCatch, try, suppressWarnings, suppressMessages

**Rule: Let operations fail loudly and immediately.**

### Rationale

1. **Fail-Fast Approach**: Errors should halt execution immediately so problems are caught early
2. **Transparent Debugging**: Full error messages and stack traces help identify root causes
3. **Data Integrity**: Silent failures can corrupt results without detection
4. **Reproducibility**: Hidden errors make workflows unreproducible

### What This Means

#### ❌ AVOID:
```r
# DON'T hide errors
tryCatch({
  result <- risky_operation()
}, error = function(e) {
  message("Something went wrong")
  return(NULL)
})

# DON'T suppress warnings
suppressWarnings(problematic_function())

# DON'T suppress messages
suppressMessages(library(package))
```

#### ✅ INSTEAD:
```r
# Let it fail if there's a problem
result <- risky_operation()

# Show all warnings
problematic_function()

# Show all messages
library(package)
```

### Exceptions

None. If an operation truly needs error handling, fix the underlying issue instead.

### When Something Fails

1. **Read the error message carefully** - it tells you what's wrong
2. **Check prerequisites** - missing packages, corrupt files, wrong paths
3. **Fix the root cause** - don't mask it with error handling
4. **Verify the fix** - ensure the operation succeeds cleanly

---

## Additional Standards

### Code Organization

- Use meaningful variable names
- Comment complex operations
- Keep chunks focused on single tasks
- Use consistent naming conventions (snake_case for variables)

### Documentation

- Document all custom functions
- Explain non-obvious operations
- Link to relevant documentation/papers

### Reproducibility

- Set seeds before random operations: `set.seed(51773)`
- Use explicit package versions via renv
- Document system requirements

---

## Integration with Cursor AI

This project includes a `.cursorrules` file at the project root that configures Cursor AI to follow these standards automatically. The AI assistant will:

- Avoid suggesting tryCatch or error suppression
- Recommend fail-fast approaches
- Suggest fixing root causes over workarounds
- Follow project naming and organization conventions

If you notice the AI suggesting code that violates these standards, remind it to check the `.cursorrules` file or this document.

---

**Last Updated:** 2025-10-22  
**Applies to:** All R and Rmd files in this project  
**Enforced by:** `.cursorrules` (for Cursor AI) and code review

