#!/usr/bin/env python3
"""Write MANIFEST.md: every published file, its size and its sha256.

Also enforces the size budget. The repo is consumed over
raw.githubusercontent.com, one file at a time, and Nextflow re-stages remote
files on every run with no shared cache -- so every megabyte here is paid again
on every job of every CI matrix. Keep it small.
"""

import hashlib
import subprocess
import sys
from pathlib import Path

import config

TOTAL_BUDGET_MB = 100
# The largest single file is the E. coli map-ont minimap2 index at 13.3 MB.
FILE_BUDGET_MB = 15


def generator_version() -> tuple[str, bool]:
    """The scripts-branch commit this dataset was built from.

    Data and generators live on different branches, so a published dataset
    cannot point at its own source by commit. Record it here instead.
    """
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"],
                                    capture_output=True, text=True,
                                    check=True).stdout.strip())
        return sha, dirty
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown", True


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    root = Path(config.DATA_DIR)
    files = sorted(p for p in root.rglob("*") if p.is_file() and not p.name.startswith("."))

    rows, total, oversized = [], 0, []
    for path in files:
        size = path.stat().st_size
        total += size
        if size > FILE_BUDGET_MB * 1024 * 1024:
            oversized.append((path, size))
        rows.append((str(path), size, sha256(path)))

    groups: dict[str, list[int]] = {}
    for path, size, _ in rows:
        key = "/".join(Path(path).parts[:3])
        groups.setdefault(key, []).append(size)

    sha, dirty = generator_version()
    lines = [
        "# Manifest",
        "",
        f"{len(rows)} files, {total / 1048576:.1f} MB total "
        f"(budget {TOTAL_BUDGET_MB} MB).",
        "",
        "| | |",
        "| --- | --- |",
        f"| generator commit | `{sha}`{' **+ uncommitted changes**' if dirty else ''} |",
        f"| seed | `{config.SEED}` |",
        f"| published base URL | `{config.BASE_URL}` |",
        "",
        "| area | files | MB |",
        "| --- | ---: | ---: |",
    ]
    for key in sorted(groups):
        sizes = groups[key]
        lines.append(f"| `{key}` | {len(sizes)} | {sum(sizes) / 1048576:.1f} |")

    lines += ["", "## Files", "", "| file | bytes | sha256 |", "| --- | ---: | --- |"]
    lines += [f"| `{path}` | {size} | `{digest[:16]}` |" for path, size, digest in rows]

    Path("MANIFEST.md").write_text("\n".join(lines) + "\n")
    print(f"MANIFEST.md: {len(rows)} files, {total / 1048576:.1f} MB")

    failed = False
    for path, size in oversized:
        print(f"OVER FILE BUDGET: {path} is {size / 1048576:.1f} MB "
              f"(limit {FILE_BUDGET_MB} MB)")
        failed = True
    if total > TOTAL_BUDGET_MB * 1024 * 1024:
        print(f"OVER TOTAL BUDGET: {total / 1048576:.1f} MB "
              f"(limit {TOTAL_BUDGET_MB} MB)")
        failed = True
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
