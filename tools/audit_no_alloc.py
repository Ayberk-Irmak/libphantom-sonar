#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Source-level check that the library allocates nothing.

Comments and string literals are stripped before scanning, so prose like
"free to get right" does not trip the audit -- an audit that cries wolf gets
disabled, which is worse than no audit.

This is the cheap first pass. The authoritative check is `nm` over the built
archive (see the link-audit step in .github/workflows/ci.yml): source greps can
be fooled, symbol tables cannot.

    python3 tools/audit_no_alloc.py [paths...]      default: include/ src/
"""

import re
import sys
from pathlib import Path

BANNED = [
    (r"\bnew\b(?!\s*line)", "operator new"),
    (r"\bdelete\b", "operator delete"),
    (r"\bmalloc\b", "malloc"),
    (r"\bcalloc\b", "calloc"),
    (r"\brealloc\b", "realloc"),
    (r"\bfree\s*\(", "free()"),
    (r"std::vector\b", "std::vector"),
    (r"std::string\b", "std::string"),
    (r"std::(unordered_)?(map|set)\b", "std::map / std::set"),
    (r"std::(shared_ptr|unique_ptr)\b", "smart pointer"),
    (r"std::function\b", "std::function"),
]

BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)
LINE_COMMENT = re.compile(r"//[^\n]*")
STRING_LIT = re.compile(r'"(?:[^"\\\n]|\\.)*"')
CHAR_LIT = re.compile(r"'(?:[^'\\\n]|\\.)*'")


def strip_noncode(text: str) -> str:
    """Blank out comments and literals, preserving line numbering."""

    def blank(match):
        return re.sub(r"[^\n]", " ", match.group(0))

    for pattern in (BLOCK_COMMENT, LINE_COMMENT, STRING_LIT, CHAR_LIT):
        text = pattern.sub(blank, text)
    return text


def audit(paths):
    findings = []
    scanned = 0
    for root in paths:
        root = Path(root)
        files = [root] if root.is_file() else sorted(root.rglob("*.[ch]pp"))
        for path in files:
            scanned += 1
            code = strip_noncode(path.read_text(encoding="utf-8"))
            for lineno, line in enumerate(code.splitlines(), 1):
                for pattern, label in BANNED:
                    if re.search(pattern, line):
                        findings.append((path, lineno, label, line.strip()))
    return scanned, findings


def main(argv):
    paths = argv[1:] or ["include", "src"]
    scanned, findings = audit(paths)

    for path, lineno, label, line in findings:
        print(f"{path}:{lineno}: {label}\n    {line}", file=sys.stderr)

    if findings:
        print(
            f"\nFAIL: {len(findings)} allocating construct(s) in {scanned} file(s)",
            file=sys.stderr,
        )
        return 1

    print(f"clean: no allocating constructs in {scanned} library source file(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
