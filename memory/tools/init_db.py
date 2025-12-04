#!/usr/bin/env python3
"""
Initialize the Synchronism research memory database.

Usage:
    python init_db.py           # Create fresh database
    python init_db.py --reset   # Reset existing database
"""

import sqlite3
import argparse
from pathlib import Path

SCHEMA_PATH = Path(__file__).parent.parent / "schema.sql"
DB_PATH = Path(__file__).parent.parent / "knowledge.db"


def init_database(reset: bool = False):
    """Initialize the database from schema."""

    if DB_PATH.exists():
        if reset:
            print(f"Resetting database at {DB_PATH}")
            DB_PATH.unlink()
        else:
            print(f"Database already exists at {DB_PATH}")
            print("Use --reset to recreate it.")
            return

    if not SCHEMA_PATH.exists():
        print(f"Schema not found at {SCHEMA_PATH}")
        return

    # Read schema
    with open(SCHEMA_PATH, 'r') as f:
        schema = f.read()

    # Create database
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    try:
        cursor.executescript(schema)
        conn.commit()
        print(f"✅ Database created at {DB_PATH}")

        # Verify tables
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]
        print(f"   Tables: {', '.join(tables)}")

    except sqlite3.Error as e:
        print(f"❌ Error creating database: {e}")
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description="Initialize Synchronism memory database")
    parser.add_argument("--reset", "-r", action="store_true", help="Reset existing database")
    args = parser.parse_args()

    init_database(reset=args.reset)


if __name__ == "__main__":
    main()
