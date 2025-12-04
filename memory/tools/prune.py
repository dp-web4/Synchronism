#!/usr/bin/env python3
"""
Synchronism Memory Prune Tool

SNARC-based salience pruning for the research memory database.
Removes low-salience findings to keep memory focused and relevant.

Usage:
    python prune.py --dry-run              # Show what would be pruned
    python prune.py --threshold 0.3        # Prune below salience 0.3
    python prune.py --keep 100             # Keep top 100 by salience
    python prune.py --decay                # Apply time-based decay
    python prune.py --superseded           # Remove superseded findings
"""

import sqlite3
import argparse
from pathlib import Path
from datetime import datetime, timedelta

DB_PATH = Path(__file__).parent.parent / "knowledge.db"


def get_connection():
    """Get database connection."""
    if not DB_PATH.exists():
        print(f"Database not found at {DB_PATH}")
        return None
    return sqlite3.connect(DB_PATH)


def show_salience_distribution():
    """Show distribution of salience scores."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    # Get distribution
    cursor.execute("""
        SELECT
            CASE
                WHEN salience < 0.2 THEN '0.0-0.2 (very low)'
                WHEN salience < 0.4 THEN '0.2-0.4 (low)'
                WHEN salience < 0.6 THEN '0.4-0.6 (medium)'
                WHEN salience < 0.8 THEN '0.6-0.8 (high)'
                ELSE '0.8-1.0 (very high)'
            END as bucket,
            COUNT(*) as count
        FROM findings
        GROUP BY bucket
        ORDER BY bucket
    """)

    results = cursor.fetchall()
    conn.close()

    print("\n=== Salience Distribution ===\n")
    for bucket, count in results:
        bar = '█' * min(count, 50)
        print(f"{bucket}: {count:3d} {bar}")


def apply_time_decay(decay_per_day: float = 0.01, dry_run: bool = True):
    """Apply time-based decay to salience scores."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    # Get findings with their age
    cursor.execute("""
        SELECT id, title, salience, created,
               julianday('now') - julianday(created) as age_days
        FROM findings
        WHERE status = 'active'
        ORDER BY salience
    """)

    results = cursor.fetchall()

    print("\n=== Time Decay Analysis ===\n")
    print(f"Decay rate: {decay_per_day} per day")
    print(f"{'Dry run' if dry_run else 'APPLYING DECAY'}\n")

    updates = []
    for id_, title, salience, created, age_days in results:
        decay = decay_per_day * age_days
        new_salience = max(0.1, salience - decay)  # Floor at 0.1

        if new_salience != salience:
            updates.append((id_, title, salience, new_salience, age_days))
            if dry_run:
                print(f"  {id_}: {salience:.3f} → {new_salience:.3f} ({age_days:.0f} days old)")

    if not dry_run and updates:
        for id_, title, old_sal, new_sal, _ in updates:
            cursor.execute("UPDATE findings SET salience = ? WHERE id = ?", (new_sal, id_))
        conn.commit()
        print(f"Updated {len(updates)} findings.")

    conn.close()


def prune_by_threshold(threshold: float, dry_run: bool = True):
    """Prune findings below salience threshold."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    cursor.execute("""
        SELECT id, title, salience, domain, finding_type, session_id
        FROM findings
        WHERE salience < ? AND status = 'active'
        ORDER BY salience
    """, (threshold,))

    results = cursor.fetchall()

    if not results:
        print(f"No findings below threshold {threshold}")
        return

    print(f"\n=== Prune by Threshold ({threshold}) ===\n")
    print(f"{'DRY RUN - ' if dry_run else ''}Found {len(results)} findings to prune:\n")

    for id_, title, salience, domain, ftype, session in results:
        print(f"  [{salience:.3f}] {title[:50]} (Session #{session})")

    if not dry_run:
        cursor.execute("""
            UPDATE findings SET status = 'deprecated'
            WHERE salience < ? AND status = 'active'
        """, (threshold,))
        conn.commit()
        print(f"\nDeprecated {len(results)} findings.")
    else:
        print(f"\nRun without --dry-run to prune.")

    conn.close()


def prune_keep_top(keep_count: int, dry_run: bool = True):
    """Keep only top N findings by salience."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    # Get count
    cursor.execute("SELECT COUNT(*) FROM findings WHERE status = 'active'")
    total = cursor.fetchone()[0]

    if total <= keep_count:
        print(f"Only {total} active findings, nothing to prune (keeping {keep_count})")
        return

    prune_count = total - keep_count

    # Get findings to prune
    cursor.execute("""
        SELECT id, title, salience
        FROM findings
        WHERE status = 'active'
        ORDER BY salience ASC
        LIMIT ?
    """, (prune_count,))

    results = cursor.fetchall()

    print(f"\n=== Prune to Keep Top {keep_count} ===\n")
    print(f"{'DRY RUN - ' if dry_run else ''}Pruning {prune_count} of {total} findings:\n")

    for id_, title, salience in results:
        print(f"  [{salience:.3f}] {title[:60]}")

    if not dry_run:
        ids_to_prune = [r[0] for r in results]
        placeholders = ','.join('?' * len(ids_to_prune))
        cursor.execute(f"""
            UPDATE findings SET status = 'deprecated'
            WHERE id IN ({placeholders})
        """, ids_to_prune)
        conn.commit()
        print(f"\nDeprecated {prune_count} findings.")
    else:
        print(f"\nRun without --dry-run to prune.")

    conn.close()


def prune_superseded(dry_run: bool = True):
    """Remove findings marked as superseded."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()

    cursor.execute("""
        SELECT id, title, salience
        FROM findings
        WHERE status = 'superseded'
    """)

    results = cursor.fetchall()

    if not results:
        print("No superseded findings to prune.")
        return

    print(f"\n=== Prune Superseded ===\n")
    print(f"{'DRY RUN - ' if dry_run else ''}Found {len(results)} superseded findings:\n")

    for id_, title, salience in results:
        print(f"  {title[:60]}")

    if not dry_run:
        cursor.execute("DELETE FROM findings WHERE status = 'superseded'")
        conn.commit()
        print(f"\nDeleted {len(results)} superseded findings.")
    else:
        print(f"\nRun without --dry-run to delete.")

    conn.close()


def boost_salience(finding_id: str, boost: float = 0.1):
    """Boost salience of a specific finding (accessed = more important)."""
    conn = get_connection()
    if not conn:
        return

    cursor = conn.cursor()
    cursor.execute("""
        UPDATE findings
        SET salience = MIN(1.0, salience + ?),
            last_accessed = CURRENT_TIMESTAMP,
            access_count = access_count + 1
        WHERE id = ?
    """, (boost, finding_id))

    if cursor.rowcount > 0:
        conn.commit()
        print(f"Boosted salience of {finding_id} by {boost}")
    else:
        print(f"Finding {finding_id} not found")

    conn.close()


def main():
    parser = argparse.ArgumentParser(description="Prune Synchronism research memory")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Show what would be pruned")
    parser.add_argument("--threshold", "-t", type=float, help="Prune below this salience")
    parser.add_argument("--keep", "-k", type=int, help="Keep only top N by salience")
    parser.add_argument("--decay", "-d", action="store_true", help="Apply time-based decay")
    parser.add_argument("--decay-rate", type=float, default=0.01, help="Decay per day")
    parser.add_argument("--superseded", "-s", action="store_true", help="Prune superseded findings")
    parser.add_argument("--distribution", action="store_true", help="Show salience distribution")
    parser.add_argument("--boost", help="Boost salience of specific finding ID")
    parser.add_argument("--boost-amount", type=float, default=0.1, help="Boost amount")

    args = parser.parse_args()

    # Default to dry run for safety
    dry_run = args.dry_run or not any([args.threshold, args.keep, args.decay, args.superseded])

    if args.distribution:
        show_salience_distribution()
    elif args.boost:
        boost_salience(args.boost, args.boost_amount)
    elif args.threshold:
        prune_by_threshold(args.threshold, dry_run)
    elif args.keep:
        prune_keep_top(args.keep, dry_run)
    elif args.decay:
        apply_time_decay(args.decay_rate, dry_run)
    elif args.superseded:
        prune_superseded(dry_run)
    else:
        show_salience_distribution()
        print("\nUse --threshold, --keep, --decay, or --superseded to prune.")
        print("Add --dry-run to preview without changes.")


if __name__ == "__main__":
    main()
