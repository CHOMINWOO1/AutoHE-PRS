"""Run the OpenFHE-only integration tests and emit a machine-readable summary."""

from __future__ import annotations

import json
import sys
import unittest


def main() -> int:
    suite = unittest.defaultTestLoader.discover(
        "tests",
        pattern="test_autohe_openfhe_integration.py",
    )
    result = unittest.TextTestRunner(stream=sys.stderr, verbosity=2).run(suite)
    payload = {
        "status": "ok" if result.wasSuccessful() else "failed",
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if result.wasSuccessful() else 2


if __name__ == "__main__":
    raise SystemExit(main())
