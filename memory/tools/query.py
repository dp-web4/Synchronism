#!/usr/bin/env python3
"""
Synchronism Memory Query Tool

Search and retrieve findings from the research memory database.
Supports text search, domain filtering, and salience-based ranking.

Usage:
    python query.py "search term"
    python query.py --domain galaxy_physics
    python query.py --type derivation
    python query.py --high-salience
    python query.py --parameters
    python query.py --predictions --untested
"""

import sqlite3
import argparse
import json
from pathlib import Path
from datetime import datetime

DB_PATH = Path(__file__).parent.parent / "knowledge.db"


def get_connection():
    """Get database connection."""
    if not DB_PATH.exists():
        print(f"Database not found at {DB_PATH}")
        print("Run: python init_db.py to create it")
        return None
    return sqlite3.connect(DB_PATH)


def search_findings(query: str = None, domain: str = None, finding_type: str = None,
                   high_salience: bool = False, limit: int = 20):
    """Search findings with various filters."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    if query:
        # Full-text search
        cursor.execute("""
            SELECT f.id, f.title, f.domain, f.finding_type, f.summary,
                   f.salience, f.validation_level, f.session_id
            FROM findings f
            JOIN findings_fts fts ON f.rowid = fts.rowid
            WHERE findings_fts MATCH ?
            ORDER BY f.salience DESC
            LIMIT ?
        """, (query, limit))
    else:
        # Filtered search
        conditions = ["1=1"]
        params = []

        if domain:
            conditions.append("domain = ?")
            params.append(domain)

        if finding_type:
            conditions.append("finding_type = ?")
            params.append(finding_type)

        if high_salience:
            conditions.append("salience > 0.5")

        params.append(limit)

        cursor.execute(f"""
            SELECT id, title, domain, finding_type, summary,
                   salience, validation_level, session_id
            FROM findings
            WHERE {' AND '.join(conditions)}
            ORDER BY salience DESC
            LIMIT ?
        """, params)

    results = cursor.fetchall()
    conn.close()

    if not results:
        print("No findings found.")
        return

    print(f"\n{'='*80}")
    print(f"Found {len(results)} findings")
    print(f"{'='*80}\n")

    for row in results:
        id_, title, domain, ftype, summary, salience, validation, session = row
        print(f"[{id_}] (Session #{session})")
        print(f"  Title: {title}")
        print(f"  Domain: {domain} | Type: {ftype} | Validation: {validation}")
        print(f"  Salience: {salience:.3f}")
        print(f"  Summary: {summary[:200]}..." if len(summary) > 200 else f"  Summary: {summary}")
        print()


def show_parameters():
    """Show all tracked parameters."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()
    cursor.execute("""
        SELECT name, symbol, value_derived, value_empirical, value_current,
               derivation_status, derivation_summary
        FROM parameters
        ORDER BY name
    """)

    results = cursor.fetchall()
    conn.close()

    if not results:
        print("No parameters tracked yet.")
        return

    print(f"\n{'='*80}")
    print("SYNCHRONISM PARAMETERS")
    print(f"{'='*80}\n")

    for row in results:
        name, symbol, derived, empirical, current, status, summary = row
        print(f"{symbol or name} ({name})")
        print(f"  Derived: {derived} | Empirical: {empirical} | Current: {current}")
        print(f"  Status: {status}")
        if summary:
            print(f"  Derivation: {summary}")
        print()


def show_predictions(untested_only: bool = False):
    """Show predictions and their status."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    if untested_only:
        cursor.execute("""
            SELECT id, title, quantitative_claim, discriminates_from,
                   status, priority, test_method
            FROM predictions
            WHERE status = 'untested'
            ORDER BY
                CASE priority
                    WHEN 'critical' THEN 1
                    WHEN 'high' THEN 2
                    WHEN 'medium' THEN 3
                    ELSE 4
                END
        """)
    else:
        cursor.execute("""
            SELECT id, title, quantitative_claim, discriminates_from,
                   status, priority, test_method
            FROM predictions
            ORDER BY
                CASE priority
                    WHEN 'critical' THEN 1
                    WHEN 'high' THEN 2
                    WHEN 'medium' THEN 3
                    ELSE 4
                END
        """)

    results = cursor.fetchall()
    conn.close()

    if not results:
        print("No predictions tracked yet.")
        return

    print(f"\n{'='*80}")
    print("SYNCHRONISM PREDICTIONS")
    print(f"{'='*80}\n")

    for row in results:
        id_, title, claim, discriminates, status, priority, method = row
        status_emoji = {
            'untested': '⏳',
            'confirmed': '✅',
            'falsified': '❌',
            'partially_tested': '🔄',
            'inconclusive': '❓'
        }.get(status, '?')

        print(f"{status_emoji} [{priority.upper()}] {title}")
        if claim:
            print(f"   Claim: {claim}")
        if discriminates:
            print(f"   Discriminates from: {discriminates}")
        if method:
            print(f"   Test: {method[:80]}...")
        print()


def show_validations():
    """Show validation test results."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()
    cursor.execute("""
        SELECT test_name, dataset, success_rate, median_chi_squared,
               sample_size, parameters_json, session_id, created
        FROM validations
        ORDER BY created DESC
        LIMIT 20
    """)

    results = cursor.fetchall()
    conn.close()

    if not results:
        print("No validations recorded yet.")
        return

    print(f"\n{'='*80}")
    print("VALIDATION RESULTS")
    print(f"{'='*80}\n")

    for row in results:
        name, dataset, success, chi2, n, params, session, created = row
        print(f"[Session #{session}] {name}")
        print(f"  Dataset: {dataset} (n={n})")
        print(f"  Success: {success*100:.1f}% | Median χ²: {chi2:.2f}")
        if params:
            print(f"  Parameters: {params}")
        print()


def show_stats():
    """Show database statistics."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    # Count findings by domain
    cursor.execute("SELECT domain, COUNT(*) FROM findings GROUP BY domain")
    domain_counts = cursor.fetchall()

    # Count findings by type
    cursor.execute("SELECT finding_type, COUNT(*) FROM findings GROUP BY finding_type")
    type_counts = cursor.fetchall()

    # Count predictions by status
    cursor.execute("SELECT status, COUNT(*) FROM predictions GROUP BY status")
    pred_counts = cursor.fetchall()

    # Total counts
    cursor.execute("SELECT COUNT(*) FROM findings")
    total_findings = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM parameters")
    total_params = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM predictions")
    total_preds = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM validations")
    total_vals = cursor.fetchone()[0]

    conn.close()

    print(f"\n{'='*80}")
    print("SYNCHRONISM MEMORY STATISTICS")
    print(f"{'='*80}\n")

    print(f"Total: {total_findings} findings, {total_params} parameters, "
          f"{total_preds} predictions, {total_vals} validations\n")

    if domain_counts:
        print("By Domain:")
        for domain, count in domain_counts:
            print(f"  {domain}: {count}")
        print()

    if type_counts:
        print("By Type:")
        for ftype, count in type_counts:
            print(f"  {ftype}: {count}")
        print()

    if pred_counts:
        print("Predictions:")
        for status, count in pred_counts:
            print(f"  {status}: {count}")


def main():
    parser = argparse.ArgumentParser(description="Query Synchronism research memory")
    parser.add_argument("query", nargs="?", help="Search term")
    parser.add_argument("--domain", "-d", help="Filter by domain")
    parser.add_argument("--type", "-t", dest="finding_type", help="Filter by finding type")
    parser.add_argument("--high-salience", "-s", action="store_true", help="Only high salience findings")
    parser.add_argument("--parameters", "-p", action="store_true", help="Show parameters")
    parser.add_argument("--predictions", action="store_true", help="Show predictions")
    parser.add_argument("--untested", "-u", action="store_true", help="Only untested predictions")
    parser.add_argument("--validations", "-v", action="store_true", help="Show validation results")
    parser.add_argument("--stats", action="store_true", help="Show statistics")
    parser.add_argument("--limit", "-l", type=int, default=20, help="Max results")

    args = parser.parse_args()

    if args.stats:
        show_stats()
    elif args.parameters:
        show_parameters()
    elif args.predictions:
        show_predictions(args.untested)
    elif args.validations:
        show_validations()
    else:
        search_findings(
            query=args.query,
            domain=args.domain,
            finding_type=args.finding_type,
            high_salience=args.high_salience,
            limit=args.limit
        )


if __name__ == "__main__":
    main()
