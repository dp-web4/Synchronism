#!/usr/bin/env python3
"""
Synchronism Memory Add Tool

Add findings, parameters, predictions, and validations to the research memory.

Usage:
    python add.py finding --interactive
    python add.py parameter --name gamma --symbol γ --value 2.0 --status first_principles
    python add.py prediction --title "Void BTFR offset" --claim "0.28 dex"
    python add.py validation --test "SPARC rotation curves" --success 0.52
"""

import sqlite3
import argparse
import json
from pathlib import Path
from datetime import datetime
import uuid

DB_PATH = Path(__file__).parent.parent / "knowledge.db"


def get_connection():
    """Get database connection."""
    if not DB_PATH.exists():
        print(f"Database not found at {DB_PATH}")
        print("Run: python init_db.py to create it")
        return None
    return sqlite3.connect(DB_PATH)


def generate_id(prefix: str, title: str) -> str:
    """Generate a readable ID from title."""
    slug = title.lower()
    slug = ''.join(c if c.isalnum() else '-' for c in slug)
    slug = '-'.join(filter(None, slug.split('-')))[:40]
    return f"{prefix}-{slug}"


def add_finding_interactive():
    """Add a finding interactively."""
    print("\n=== Add Finding ===\n")

    # Required fields
    title = input("Title: ").strip()
    if not title:
        print("Title is required.")
        return

    summary = input("Summary (1-2 sentences): ").strip()
    if not summary:
        print("Summary is required.")
        return

    print("\nDomains: cosmology, galaxy_physics, coherence_math, implementation, validation, methodology")
    domain = input("Domain: ").strip()

    print("\nTypes: derivation, validation, prediction, connection, methodology, failure, insight, question")
    finding_type = input("Type: ").strip()

    session_id = input("Session ID (e.g., 78): ").strip()

    # Optional fields
    details = input("Details (optional, press Enter to skip): ").strip() or None
    subdomain = input("Subdomain (optional): ").strip() or None

    print("\nValidation: proven, validated, theoretical, speculative, falsified")
    validation_level = input("Validation level [theoretical]: ").strip() or "theoretical"

    # SNARC scores
    print("\nSNARC scores (0.0-1.0, press Enter for default 0.5)")
    try:
        surprise = float(input("Surprise [0.5]: ").strip() or 0.5)
        novelty = float(input("Novelty [0.5]: ").strip() or 0.5)
        arousal = float(input("Arousal [0.5]: ").strip() or 0.5)
        reward = float(input("Reward [0.5]: ").strip() or 0.5)
        conflict = float(input("Conflict [0.0]: ").strip() or 0.0)
    except ValueError:
        print("Invalid score, using defaults.")
        surprise = novelty = arousal = reward = 0.5
        conflict = 0.0

    tags = input("Tags (comma-separated): ").strip()
    tags_json = json.dumps([t.strip() for t in tags.split(',')]) if tags else None

    # Generate ID
    finding_id = generate_id("f", title)

    # Insert
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()
    try:
        cursor.execute("""
            INSERT INTO findings (
                id, title, summary, details, domain, subdomain, finding_type,
                session_id, validation_level, surprise, novelty, arousal, reward, conflict, tags
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (finding_id, title, summary, details, domain, subdomain, finding_type,
              session_id, validation_level, surprise, novelty, arousal, reward, conflict, tags_json))
        conn.commit()
        print(f"\n✅ Added finding: {finding_id}")
    except sqlite3.IntegrityError as e:
        print(f"\n❌ Error: {e}")
    finally:
        conn.close()


def add_finding(title: str, summary: str, domain: str, finding_type: str,
                session_id: str = None, details: str = None, subdomain: str = None,
                validation_level: str = "theoretical", surprise: float = 0.5,
                novelty: float = 0.5, arousal: float = 0.5, reward: float = 0.5,
                conflict: float = 0.0, tags: list = None):
    """Add a finding programmatically."""
    conn = get_connection()
    if not conn:
        return None

    finding_id = generate_id("f", title)
    tags_json = json.dumps(tags) if tags else None

    cursor = conn.cursor()
    try:
        cursor.execute("""
            INSERT INTO findings (
                id, title, summary, details, domain, subdomain, finding_type,
                session_id, validation_level, surprise, novelty, arousal, reward, conflict, tags
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (finding_id, title, summary, details, domain, subdomain, finding_type,
              session_id, validation_level, surprise, novelty, arousal, reward, conflict, tags_json))
        conn.commit()
        return finding_id
    except sqlite3.IntegrityError as e:
        print(f"Error: {e}")
        return None
    finally:
        conn.close()


def add_parameter(name: str, symbol: str = None, value_derived: float = None,
                  value_empirical: float = None, value_current: float = None,
                  derivation_status: str = "unknown", derivation_session: str = None,
                  derivation_summary: str = None, physical_meaning: str = None,
                  units: str = None):
    """Add or update a parameter."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()
    try:
        cursor.execute("""
            INSERT INTO parameters (
                name, symbol, value_derived, value_empirical, value_current,
                derivation_status, derivation_session, derivation_summary,
                physical_meaning, units
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(name) DO UPDATE SET
                symbol = excluded.symbol,
                value_derived = COALESCE(excluded.value_derived, value_derived),
                value_empirical = COALESCE(excluded.value_empirical, value_empirical),
                value_current = COALESCE(excluded.value_current, value_current),
                derivation_status = excluded.derivation_status,
                derivation_session = COALESCE(excluded.derivation_session, derivation_session),
                derivation_summary = COALESCE(excluded.derivation_summary, derivation_summary),
                physical_meaning = COALESCE(excluded.physical_meaning, physical_meaning),
                units = COALESCE(excluded.units, units),
                updated = CURRENT_TIMESTAMP
        """, (name, symbol, value_derived, value_empirical, value_current,
              derivation_status, derivation_session, derivation_summary,
              physical_meaning, units))
        conn.commit()
        print(f"✅ Parameter '{name}' saved.")
    except sqlite3.Error as e:
        print(f"❌ Error: {e}")
    finally:
        conn.close()


def add_prediction(title: str, description: str, quantitative_claim: str = None,
                   test_method: str = None, required_data: str = None,
                   discriminates_from: str = None, priority: str = "medium",
                   source_session: str = None):
    """Add a prediction."""
    conn = get_connection()
    if not conn:
        return

    pred_id = generate_id("p", title)

    cursor = conn.cursor()
    try:
        cursor.execute("""
            INSERT INTO predictions (
                id, title, description, quantitative_claim, test_method,
                required_data, discriminates_from, priority, source_session
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (pred_id, title, description, quantitative_claim, test_method,
              required_data, discriminates_from, priority, source_session))
        conn.commit()
        print(f"✅ Prediction '{pred_id}' added.")
    except sqlite3.IntegrityError as e:
        print(f"❌ Error: {e}")
    finally:
        conn.close()


def add_validation(test_name: str, dataset: str, success_rate: float,
                   median_chi_squared: float = None, sample_size: int = None,
                   parameters: dict = None, comparison_model: str = None,
                   comparison_result: str = None, session_id: str = None,
                   notes: str = None):
    """Add a validation result."""
    conn = get_connection()
    if not conn:
        return

    val_id = generate_id("v", f"{test_name}-{dataset}")

    cursor = conn.cursor()
    try:
        cursor.execute("""
            INSERT INTO validations (
                id, test_name, dataset, success_rate, median_chi_squared,
                sample_size, parameters_json, comparison_model, comparison_result,
                session_id, notes
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (val_id, test_name, dataset, success_rate, median_chi_squared,
              sample_size, json.dumps(parameters) if parameters else None,
              comparison_model, comparison_result, session_id, notes))
        conn.commit()
        print(f"✅ Validation '{val_id}' recorded.")
    except sqlite3.IntegrityError as e:
        print(f"❌ Error: {e}")
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description="Add to Synchronism research memory")
    subparsers = parser.add_subparsers(dest="command", help="What to add")

    # Finding subcommand
    finding_parser = subparsers.add_parser("finding", help="Add a finding")
    finding_parser.add_argument("--interactive", "-i", action="store_true", help="Interactive mode")
    finding_parser.add_argument("--title", help="Finding title")
    finding_parser.add_argument("--summary", help="Brief summary")
    finding_parser.add_argument("--domain", help="Domain")
    finding_parser.add_argument("--type", dest="finding_type", help="Finding type")
    finding_parser.add_argument("--session", help="Session ID")

    # Parameter subcommand
    param_parser = subparsers.add_parser("parameter", help="Add a parameter")
    param_parser.add_argument("--name", required=True, help="Parameter name")
    param_parser.add_argument("--symbol", help="Symbol (e.g., γ)")
    param_parser.add_argument("--derived", type=float, help="Derived value")
    param_parser.add_argument("--empirical", type=float, help="Empirical value")
    param_parser.add_argument("--current", type=float, help="Current value")
    param_parser.add_argument("--status", help="Derivation status")
    param_parser.add_argument("--session", help="Derivation session")
    param_parser.add_argument("--meaning", help="Physical meaning")

    # Prediction subcommand
    pred_parser = subparsers.add_parser("prediction", help="Add a prediction")
    pred_parser.add_argument("--title", required=True, help="Prediction title")
    pred_parser.add_argument("--description", required=True, help="Description")
    pred_parser.add_argument("--claim", help="Quantitative claim")
    pred_parser.add_argument("--test", help="Test method")
    pred_parser.add_argument("--data", help="Required data")
    pred_parser.add_argument("--discriminates", help="What it discriminates from")
    pred_parser.add_argument("--priority", default="medium", help="Priority level")
    pred_parser.add_argument("--session", help="Source session")

    # Validation subcommand
    val_parser = subparsers.add_parser("validation", help="Add a validation result")
    val_parser.add_argument("--test", required=True, help="Test name")
    val_parser.add_argument("--dataset", required=True, help="Dataset used")
    val_parser.add_argument("--success", type=float, required=True, help="Success rate (0-1)")
    val_parser.add_argument("--chi2", type=float, help="Median chi-squared")
    val_parser.add_argument("--n", type=int, help="Sample size")
    val_parser.add_argument("--session", help="Session ID")

    args = parser.parse_args()

    if args.command == "finding":
        if args.interactive:
            add_finding_interactive()
        elif args.title and args.summary and args.domain and args.finding_type:
            add_finding(args.title, args.summary, args.domain, args.finding_type,
                       session_id=args.session)
        else:
            print("Use --interactive or provide --title, --summary, --domain, --type")

    elif args.command == "parameter":
        add_parameter(
            name=args.name,
            symbol=args.symbol,
            value_derived=args.derived,
            value_empirical=args.empirical,
            value_current=args.current,
            derivation_status=args.status or "unknown",
            derivation_session=args.session,
            physical_meaning=args.meaning
        )

    elif args.command == "prediction":
        add_prediction(
            title=args.title,
            description=args.description,
            quantitative_claim=args.claim,
            test_method=args.test,
            required_data=args.data,
            discriminates_from=args.discriminates,
            priority=args.priority,
            source_session=args.session
        )

    elif args.command == "validation":
        add_validation(
            test_name=args.test,
            dataset=args.dataset,
            success_rate=args.success,
            median_chi_squared=args.chi2,
            sample_size=args.n,
            session_id=args.session
        )

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
