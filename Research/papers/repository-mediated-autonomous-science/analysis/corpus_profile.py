#!/usr/bin/env python3
"""Profile claim-promotion and claim-correction language in a Git corpus.

This is a deliberately shallow descriptive measure. It classifies commit
subjects, not scientific validity. The dictionaries are printed with the
result so later readers can reproduce or challenge the coding.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path


PROMOTION_PATTERN = re.compile(
    r"\b(discover(?:y|ed)?|confirm(?:ed|s)?|validat(?:e|ed|ion)|complete|"
    r"proof|law|universal|breakthrough|success|supported|resolved)\b",
    re.IGNORECASE,
)

CORRECTION_PATTERN = re.compile(
    r"\b(falsif(?:y|ied|ication)|fail(?:s|ed|ure)?|overturn(?:ed|s)?|"
    r"correct(?:ed|ion)?|refut(?:e|ed|ation)|artifact|illusion|not a|no |"
    r"collapse|dissolve(?:s|d)?|reanalysis|cannot|never)\b",
    re.IGNORECASE,
)


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args],
        text=True,
        stderr=subprocess.DEVNULL,
    ).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def profile(
    raw_log: str,
    *,
    repository_head: str,
    tracked_file_count: int,
    log_source: str,
    log_sha256: str | None = None,
) -> dict[str, object]:
    records = []
    classes = Counter()
    months = Counter()

    for line in raw_log.splitlines():
        commit, authored_at, subject = line.split("\t", 2)
        promotion = bool(PROMOTION_PATTERN.search(subject))
        correction = bool(CORRECTION_PATTERN.search(subject))
        label = (
            "both"
            if promotion and correction
            else "promotion_only"
            if promotion
            else "correction_only"
            if correction
            else "neither"
        )
        classes[label] += 1
        months[authored_at[:7]] += 1
        records.append(
            {
                "commit": commit,
                "authored_at": authored_at,
                "subject": subject,
                "promotion_marker": promotion,
                "correction_marker": correction,
                "class": label,
            }
        )

    ordered = sorted(records, key=lambda record: record["authored_at"])
    first = ordered[0] if ordered else None
    last = ordered[-1] if ordered else None

    return {
        "repository_head": repository_head,
        "commit_count": len(records),
        "tracked_file_count": tracked_file_count,
        "first_commit": first,
        "last_commit": last,
        "commit_subject_classes": dict(classes),
        "commits_by_month": dict(sorted(months.items())),
        "log_source": log_source,
        "log_sha256": log_sha256,
        "promotion_regex": PROMOTION_PATTERN.pattern,
        "correction_regex": CORRECTION_PATTERN.pattern,
        "interpretation_warning": (
            "Commit-subject discourse only. Counts are not discovery, error, "
            "or research-quality rates; named lineages require manual audit."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--repo",
        type=Path,
        help="Path to the gnosis-research Git repository.",
    )
    source.add_argument(
        "--log-file",
        type=Path,
        help="Frozen TSV from: git log --reverse --format=%%H%%x09%%aI%%x09%%s",
    )
    parser.add_argument(
        "--repository-head",
        help="Required with --log-file: full commit hash represented by the snapshot.",
    )
    parser.add_argument(
        "--tracked-file-count",
        type=int,
        help="Required with --log-file: git ls-files count at the snapshot.",
    )
    args = parser.parse_args()

    if args.repo:
        repo = args.repo.resolve()
        raw_log = git(repo, "log", "--reverse", "--format=%H%x09%aI%x09%s")
        result = profile(
            raw_log,
            repository_head=git(repo, "rev-parse", "HEAD"),
            tracked_file_count=len(git(repo, "ls-files").splitlines()),
            log_source=str(repo),
        )
    else:
        if args.repository_head is None or args.tracked_file_count is None:
            parser.error(
                "--log-file requires --repository-head and --tracked-file-count"
            )
        log_file = args.log_file.resolve()
        result = profile(
            log_file.read_text(encoding="utf-8"),
            repository_head=args.repository_head,
            tracked_file_count=args.tracked_file_count,
            log_source=str(log_file),
            log_sha256=sha256(log_file),
        )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
