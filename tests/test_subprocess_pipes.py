#!/usr/bin/env python3
"""Guard against the check_call(..., stdout=PIPE) deadlock.

subprocess.check_call() and .call() wait for the child without ever reading
the pipes they created, so the child blocks forever once it writes more than
the 64 KB pipe buffer. Python's own docs say not to combine them with PIPE.

This bit dante_tir for real: `mmseqs easy-cluster` (a chain of chatty
sub-commands) hung in pipe_write while dante_tir.py sat in wait(), on inputs
big enough to cross the buffer -- i.e. exactly the large genomes that matter.

Use subprocess.run(), which drains the pipes while waiting, or DEVNULL when the
output is not wanted. This test walks the AST of every module in the repo so a
future call site cannot reintroduce the pattern.
"""

import ast
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BLOCKING = {"check_call", "call"}


def is_pipe(node):
    """True for `subprocess.PIPE` / a bare `PIPE` name."""
    if isinstance(node, ast.Attribute) and node.attr == "PIPE":
        return True
    return isinstance(node, ast.Name) and node.id == "PIPE"


def offending_calls(path):
    with open(path) as fh:
        tree = ast.parse(fh.read(), filename=path)
    bad = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if name not in BLOCKING:
            continue
        piped = [kw.arg for kw in node.keywords
                 if kw.arg in ("stdout", "stderr") and is_pipe(kw.value)]
        if piped:
            bad.append((node.lineno, name, ",".join(piped)))
    return bad


def main():
    failures = []
    checked = 0
    for entry in sorted(os.listdir(ROOT)):
        if not entry.endswith(".py"):
            continue
        path = os.path.join(ROOT, entry)
        checked += 1
        for lineno, name, streams in offending_calls(path):
            failures.append("%s:%d subprocess.%s() with %s=PIPE"
                            % (entry, lineno, name, streams))

    if failures:
        print("FAIL: undrained pipes can deadlock the child:")
        for f in failures:
            print("  " + f)
        print("  use subprocess.run(...) (it drains while waiting) "
              "or DEVNULL if the output is unwanted")
        return 1
    print("  no blocking subprocess call pipes its child's output (%d modules)"
          % checked)
    print("test_subprocess_pipes: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
