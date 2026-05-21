# Security Policy

## Reporting a Vulnerability

Please do not open a public issue for a potential security vulnerability. Contact the maintainer privately, or use GitHub private vulnerability reporting if it is enabled for the repository.

Include:

- A concise description of the issue.
- Steps to reproduce it.
- Affected versions or commits, if known.
- Any relevant input files, with sensitive data removed.

## Supported Versions

ChromoSort is under active pre-1.0 development. Security fixes are expected to land on `main` first and then be included in the next tagged release.

## Scope

ChromoSort is a local command-line toolkit. The main security concerns are untrusted input files, generated HTML review artifacts, and dependency vulnerabilities. Treat HTML files generated from untrusted inputs as files to review locally, and avoid serving them on shared infrastructure without normal web-security precautions.
