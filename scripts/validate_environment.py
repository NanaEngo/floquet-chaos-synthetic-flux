#!/usr/bin/env python3
"""Validate and record the local FL_QOM reconstruction environment."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REQUIRED = ("numpy", "scipy", "sympy", "matplotlib", "h5py", "joblib", "qutip", "quimb")

def git_revision(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None

def git_status(path: Path) -> list[str]:
    try:
        output = subprocess.check_output(
            ["git", "-C", str(path), "status", "--short"], text=True
        )
        return [line for line in output.splitlines() if line]
    except (OSError, subprocess.CalledProcessError):
        return []

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fl-qom-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.fl_qom_root.resolve()
    if not root.is_dir():
        print(f"ERROR: FL_QOM root does not exist: {root}", file=sys.stderr)
        return 2

    modules = {}
    missing = []
    for name in REQUIRED:
        try:
            module = importlib.import_module(name)
            modules[name] = getattr(module, "__version__", "installed")
        except Exception as exc:  # pragma: no cover - environment dependent
            modules[name] = f"ERROR: {exc}"
            missing.append(name)

    record = {
        "status": "PASS" if not missing else "FAIL",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "executable": sys.executable,
        "fl_qom_root": str(root),
        "fl_qom_revision": git_revision(root),
        "fl_qom_worktree_status": git_status(root),
        "compatibility_patch_files": {
            "src/fl_qom/analysis/non_gaussian.py": str(root / "src/fl_qom/analysis/non_gaussian.py"),
            "src/fl_qom/analysis/quantum/master_equation_solver.py": str(root / "src/fl_qom/analysis/quantum/master_equation_solver.py"),
            "tests/test_trackb_benchmark_consistency.py": str(root / "tests/test_trackb_benchmark_consistency.py"),
        },
        "modules": modules,
        "missing_or_failed_modules": missing,
    }
    pyproject = root / "pyproject.toml"
    if pyproject.is_file():
        record["fl_qom_pyproject_sha256"] = sha256(pyproject)
    patch_files = [
        root / "src/fl_qom/analysis/non_gaussian.py",
        root / "src/fl_qom/analysis/quantum/master_equation_solver.py",
        root / "tests/test_trackb_benchmark_consistency.py",
    ]
    record["compatibility_patch_sha256"] = {
        str(path.relative_to(root)): sha256(path)
        for path in patch_files if path.is_file()
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2, sort_keys=True))
    return 0 if not missing else 1

if __name__ == "__main__":
    raise SystemExit(main())
