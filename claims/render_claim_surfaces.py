#!/usr/bin/env python3
"""Validate the claim ledger and deterministically render v2 public surfaces."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "claims" / "claim_ledger.json"
FREEZE_PATH = ROOT / "claims" / "v1-freeze-manifest.json"
SNAPSHOT_DIR_RELATIVE = "claims/v1-snapshot"
README_OUTPUT = ROOT / "claims" / "generated" / "README-claims.md"
V2_OUTPUT = ROOT / "docs" / "whitepaper-v2" / "Synchronism_Whitepaper_v2.md"

REQUIRED_CLAIM_FIELDS = {
    "id",
    "title",
    "domain",
    "claim_level",
    "status",
    "novelty",
    "statement",
    "public_summary",
    "evidence",
    "source_refs",
    "updated",
}
OPTIONAL_CLAIM_FIELDS = {"falsifier_or_next_test", "supersedes"}
STATUSES = {
    "confirmed",
    "audited",
    "established",
    "supported",
    "open",
    "underdetermined",
    "null",
    "refuted",
    "erroneous",
    "historical",
}
DOMAINS = {"program", "ontology", "physics", "applied", "methodology"}
CLAIM_LEVELS = {"program", "umbrella", "realization", "method", "application"}
NOVELTY_STATES = {"novel", "not_novel", "unresolved", "not_applicable"}
CLAIM_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9._-]+$")
DATE_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}$")
STATUS_LABELS = {
    "confirmed": "Confirmed",
    "audited": "Audited",
    "established": "Established",
    "supported": "Supported (bounded)",
    "open": "Open",
    "underdetermined": "Underdetermined",
    "null": "Null result",
    "refuted": "Refuted realization",
    "erroneous": "Erroneous",
    "historical": "Historical",
}


def load_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def validate_ledger(ledger: dict[str, object]) -> list[dict[str, object]]:
    if ledger.get("schema_version") != 1:
        raise ValueError("unsupported claim-ledger schema_version")
    if set(ledger.get("status_definitions", {})) != STATUSES:
        raise ValueError("status_definitions must exactly cover the allowed statuses")
    if not DATE_PATTERN.fullmatch(str(ledger.get("as_of", ""))):
        raise ValueError("as_of must use YYYY-MM-DD")
    claims = ledger.get("claims")
    if not isinstance(claims, list) or not claims:
        raise ValueError("claims must be a non-empty list")

    ids: set[str] = set()
    for index, claim in enumerate(claims):
        if not isinstance(claim, dict):
            raise ValueError(f"claim {index} is not an object")
        missing = REQUIRED_CLAIM_FIELDS - claim.keys()
        if missing:
            raise ValueError(f"claim {index} missing fields: {sorted(missing)}")
        unknown = claim.keys() - REQUIRED_CLAIM_FIELDS - OPTIONAL_CLAIM_FIELDS
        if unknown:
            raise ValueError(f"claim {index} has unknown fields: {sorted(unknown)}")
        claim_id = claim["id"]
        if (
            not isinstance(claim_id, str)
            or not CLAIM_ID_PATTERN.fullmatch(claim_id)
            or claim_id in ids
        ):
            raise ValueError(f"claim id is missing or duplicated: {claim_id!r}")
        ids.add(claim_id)
        if claim["status"] not in STATUSES:
            raise ValueError(f"{claim_id}: unknown status {claim['status']!r}")
        if claim["domain"] not in DOMAINS:
            raise ValueError(f"{claim_id}: unknown domain {claim['domain']!r}")
        if claim["claim_level"] not in CLAIM_LEVELS:
            raise ValueError(
                f"{claim_id}: unknown claim level {claim['claim_level']!r}"
            )
        if claim["novelty"] not in NOVELTY_STATES:
            raise ValueError(f"{claim_id}: unknown novelty {claim['novelty']!r}")
        if not DATE_PATTERN.fullmatch(str(claim["updated"])):
            raise ValueError(f"{claim_id}: updated must use YYYY-MM-DD")
        for list_field in ("evidence", "source_refs"):
            value = claim[list_field]
            if not isinstance(value, list) or not value:
                raise ValueError(f"{claim_id}: {list_field} must be non-empty")
    return claims


def verify_v1_freeze() -> None:
    manifest = load_json(FREEZE_PATH)
    frozen_at_commit = manifest.get("frozen_at_commit")
    if not isinstance(frozen_at_commit, str):
        raise ValueError("v1 freeze manifest needs frozen_at_commit")

    # A hash freeze can only be applied to a path nothing else writes. Twice now
    # the manifest froze a *living* path and the gate went red on a legitimate
    # update: 2026-07-24 (the live PREDICTIONS.md, written by the research lane)
    # and 2026-07-25 (docs/whitepaper/Synchronism_Whitepaper.pdf, rewritten by
    # build_whitepaper.yml, which does `git add docs/ whitepaper/build/` on every
    # push to main -- that rebuild changed 14 bytes of PDF metadata with the
    # source .md byte-identical, so the gate failed on a content-free rebuild).
    # Both were the same class of error, so the rule is now enforced in code
    # rather than remembered: every frozen artifact lives under the snapshot
    # directory, which no workflow and no research lane writes.
    stray = sorted(
        path
        for path in manifest["files"]
        if not path.startswith(f"{SNAPSHOT_DIR_RELATIVE}/")
    )
    if stray:
        raise ValueError(
            "v1 freeze manifest may only freeze snapshot copies under "
            f"{SNAPSHOT_DIR_RELATIVE}/; these are live paths that other writers "
            "own, so freezing them makes the gate fail on legitimate updates: "
            + ", ".join(stray)
            + ". Fix: `git show <frozen_at_commit>:<path> > "
            f"{SNAPSHOT_DIR_RELATIVE}/<basename>`, then point the manifest at "
            "the snapshot (the sha256 is unchanged, so the evidence guarantee "
            "is unchanged)."
        )

    failures = []
    for relative, expected in manifest["files"].items():
        try:
            blob = subprocess.check_output(
                ["git", "-C", str(ROOT), "show", f"{frozen_at_commit}:{relative}"],
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as error:
            detail = error.stderr.decode("utf-8", errors="replace").strip()
            failures.append(
                f"{relative}: unavailable at {frozen_at_commit} ({detail})"
            )
            continue

        actual = hashlib.sha256(blob).hexdigest()
        if actual != expected:
            failures.append(
                f"{relative}@{frozen_at_commit}: expected {expected}, found {actual}"
            )
            continue

        # The manifest freezes committed Git objects, not checkout-transformed
        # bytes. Separately require the current index/worktree to remain
        # equivalent under Git's configured clean filters.
        unchanged = subprocess.run(
            [
                "git",
                "-C",
                str(ROOT),
                "diff",
                "--quiet",
                frozen_at_commit,
                "--",
                relative,
            ],
            check=False,
        )
        if unchanged.returncode == 1:
            failures.append(f"{relative}: differs from frozen commit")
        elif unchanged.returncode != 0:
            failures.append(
                f"{relative}: git diff failed with {unchanged.returncode}"
            )
    if failures:
        raise ValueError("v1 freeze violation:\n" + "\n".join(failures))


def render_readme(ledger: dict[str, object], claims: list[dict[str, object]]) -> str:
    rows = [
        "# Synchronism claim status",
        "",
        "<!-- GENERATED by claims/render_claim_surfaces.py; do not hand-edit. -->",
        "",
        f"**Ledger date:** {ledger['as_of']}",
        "**Scope:** v2 public claims currently migrated to the machine-readable ledger.",
        "",
        "This surface keeps two propositions separate: several concrete realizations",
        "are refuted or erroneous, while the umbrella ontology remains underdetermined.",
        "Applied usefulness is evidence of generative design value, not confirmation of",
        "fundamental physics.",
        "",
        "| Claim | Status | Public statement |",
        "|---|---|---|",
    ]
    for claim in claims:
        rows.append(
            f"| `{claim['id']}` | **{STATUS_LABELS[claim['status']]}** | "
            f"{claim['public_summary']} |"
        )
    rows.extend(
        [
            "",
            "## Confirmation policy",
            "",
            "A preregistered prediction can receive confirmatory evidence from archival",
            "data when the hypothesis, parameters, analysis, and success criterion were",
            "fixed before exposure to a genuine holdout. Instrument ownership is not the",
            "boundary between confirmation and refutation.",
            "",
            "## Provenance",
            "",
            "The source of this page is [`claims/claim_ledger.json`](../claim_ledger.json).",
            "The previous README, prediction ledger, and complete whitepaper are preserved",
            "as historical v1 under [`claims/v1-freeze-manifest.json`](../v1-freeze-manifest.json).",
            "",
        ]
    )
    return "\n".join(rows)


def render_v2(ledger: dict[str, object], claims: list[dict[str, object]]) -> str:
    groups = [
        ("Ontology and commitments", {"ontology"}),
        ("Empirical results and corrections", {"physics", "program"}),
        ("Applied systems work", {"applied"}),
        ("Research methodology", {"methodology"}),
    ]
    rows = [
        "# Synchronism: epistemic-status edition (v2)",
        "",
        "<!-- GENERATED by claims/render_claim_surfaces.py; do not hand-edit. -->",
        "",
        f"**Generated from claim ledger dated {ledger['as_of']}.**",
        "",
        "## Reader's contract",
        "",
        "This is a concise status edition, not a replacement history. The historical",
        "v1 manuscript is intentionally preserved because its layers of assertion,",
        "testing, correction, and recorrection are part of the research record.",
        "",
        "Synchronism currently has zero confirmed novel physics predictions. Several",
        "specific realizations have been refuted or found erroneous. That does not, by",
        "itself, refute every realization compatible with the umbrella ontology; the",
        "umbrella remains underdetermined rather than confirmed. Separately, the",
        "vocabulary has demonstrated applied value in software and governance design.",
        "",
    ]
    for heading, domains in groups:
        rows.extend([f"## {heading}", ""])
        for claim in claims:
            if claim["domain"] not in domains:
                continue
            rows.extend(
                [
                    f"### {claim['title']}",
                    "",
                    f"**Status:** {STATUS_LABELS[claim['status']]}",
                    f"**Novelty:** {claim['novelty'].replace('_', ' ')}",
                    "",
                    str(claim["statement"]),
                    "",
                    "**Evidence recorded in the ledger:**",
                    "",
                ]
            )
            rows.extend(f"- {item}" for item in claim["evidence"])
            if claim.get("falsifier_or_next_test"):
                rows.extend(
                    [
                        "",
                        f"**Next discriminating test:** {claim['falsifier_or_next_test']}",
                    ]
                )
            rows.extend(["", "**Sources:** " + "; ".join(claim["source_refs"]), ""])
    rows.extend(
        [
            "## Open-bet discipline",
            "",
            "A result is confirmatory when a specific, temporally prior claim meets a",
            "prospectively fixed criterion on genuine holdout evidence. Public archival",
            "data can satisfy that standard. Retrospective parameter choice, model-family",
            "selection, cuts, or success criteria make the result exploratory instead.",
            "",
            "A useful next registration is the proposed SPARC-by-Cassini squeeze: define",
            "the transition-function family and parameter interval, freeze the joint",
            "analysis and every outcome branch, then evaluate the existing independent",
            "datasets. The result may confirm the registered compatibility region, refute",
            "it, or leave the family underpowered; all three branches must be stated first.",
            "",
            "## Versioning and lineage",
            "",
            "This document is generated from `claims/claim_ledger.json`. It may be updated",
            "only by changing a claim record and regenerating. The v1 preservation hashes",
            "are in `claims/v1-freeze-manifest.json`; a failed hash check is a review",
            "blocker, not an invitation to refresh the baseline silently.",
            "",
        ]
    )
    return "\n".join(rows)


def write_or_check(path: Path, content: str, check: bool) -> bool:
    content = content.rstrip() + "\n"
    if check:
        if not path.exists() or path.read_text(encoding="utf-8") != content:
            print(f"stale generated surface: {path.relative_to(ROOT)}", file=sys.stderr)
            return False
        return True
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate and compare generated outputs without writing.",
    )
    args = parser.parse_args()

    ledger = load_json(LEDGER_PATH)
    claims = validate_ledger(ledger)

    # The freeze check is an evidence-integrity ALARM, not a precondition for
    # rendering: the generated surfaces are a deterministic function of the
    # ledger and have no dependency on the v1 artifacts. Letting the alarm abort
    # the renderer is what turned one legitimate prediction update into a
    # two-day outage of the README and whitepaper-v2 surfaces (2026-07-24/25).
    # So: raise in --check (CI must fail on a violation), but in the write path
    # report loudly, render anyway, and exit non-zero at the end.
    freeze_ok = True
    try:
        verify_v1_freeze()
    except ValueError as error:
        if args.check:
            raise
        freeze_ok = False
        print(f"WARNING: {error}", file=sys.stderr)
        print(
            "WARNING: rendering anyway -- the generated surfaces do not depend "
            "on the v1 artifacts. Exit status will be non-zero.",
            file=sys.stderr,
        )

    results = [
        write_or_check(README_OUTPUT, render_readme(ledger, claims), args.check),
        write_or_check(V2_OUTPUT, render_v2(ledger, claims), args.check),
    ]
    if not all(results):
        return 1
    action = "checked" if args.check else "rendered"
    freeze_note = "v1 freeze verified" if freeze_ok else "V1 FREEZE VIOLATED (see above)"
    print(f"{action} {len(claims)} claims; {freeze_note}")
    return 0 if freeze_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
