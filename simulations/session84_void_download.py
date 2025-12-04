#!/usr/bin/env python3
"""
Session #84: Void Catalog Download

Downloads the Douglass+ 2023 void catalog from VizieR for BTFR analysis.

Data source: J/ApJS/265/7
- 776,500 galaxies with void membership
- Multiple tables for different cosmologies and pruning methods

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
Session: #84
"""

import numpy as np
from pathlib import Path
import json
from datetime import datetime

try:
    from astroquery.vizier import Vizier
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False


def download_void_catalog():
    """
    Download Douglass+ 2023 void catalog from VizieR.
    """
    print("=" * 70)
    print("SESSION #84: VOID CATALOG DOWNLOAD (Douglass+ 2023)")
    print("=" * 70)

    if not ASTROQUERY_AVAILABLE:
        print("\nERROR: astroquery not available")
        return None

    # Configure Vizier
    print("\nConfiguring VizieR query...")
    v = Vizier(
        columns=['*'],  # Get all columns
        row_limit=-1
    )

    print("Downloading void catalog (J/ApJS/265/7)...")
    print("This catalog has multiple tables - getting main galaxy table...")

    try:
        # Query all tables in the catalog
        catalogs = v.get_catalogs('J/ApJS/265/7')

        if catalogs is None or len(catalogs) == 0:
            print("ERROR: No data returned from VizieR")
            return None

        print(f"\nFound {len(catalogs)} tables in catalog")
        for i, cat in enumerate(catalogs):
            print(f"  Table {i}: {len(cat)} rows, columns: {cat.colnames[:5]}...")

        return catalogs

    except Exception as e:
        print(f"ERROR downloading catalog: {e}")
        return None


def analyze_void_tables(catalogs):
    """
    Analyze the void catalog tables.
    """
    print("\n" + "=" * 70)
    print("VOID CATALOG ANALYSIS")
    print("=" * 70)

    tables_info = []

    for i, cat in enumerate(catalogs):
        info = {
            'index': i,
            'n_rows': len(cat),
            'columns': list(cat.colnames)
        }
        tables_info.append(info)

        print(f"\n--- Table {i} ---")
        print(f"Rows: {len(cat)}")
        print(f"Columns: {cat.colnames}")

        # Sample first few rows
        if len(cat) > 0:
            print("Sample data (first 3 rows):")
            for j in range(min(3, len(cat))):
                print(f"  Row {j}: {dict(zip(cat.colnames[:3], [cat[col][j] for col in cat.colnames[:3]]))}")

    return tables_info


def find_galaxy_table(catalogs):
    """
    Find the main galaxy table with void membership.
    """
    print("\n" + "=" * 70)
    print("FINDING GALAXY TABLE")
    print("=" * 70)

    # The galaxy table should be the largest one
    largest_idx = 0
    largest_size = 0

    for i, cat in enumerate(catalogs):
        if len(cat) > largest_size:
            largest_size = len(cat)
            largest_idx = i

    print(f"\nLargest table: {largest_idx} with {largest_size} rows")

    # Check if it looks like galaxy data (should have coordinates)
    galaxy_table = catalogs[largest_idx]
    has_coords = any('RA' in col.upper() or 'DEC' in col.upper() or col in ['x', 'y', 'z']
                     for col in galaxy_table.colnames)

    if has_coords:
        print("Table appears to contain galaxy coordinates ✓")
    else:
        print("WARNING: Table may not be galaxy data")
        print(f"Columns: {galaxy_table.colnames}")

    return galaxy_table, largest_idx


def save_void_catalog(galaxy_table, tables_info):
    """
    Save the void galaxy catalog to local files.
    """
    output_dir = Path(__file__).parent / 'void_data'
    output_dir.mkdir(exist_ok=True)

    # Save main galaxy table
    galaxy_path = output_dir / 'void_galaxies.csv'

    data = {}
    for col in galaxy_table.colnames:
        try:
            data[col] = np.array(galaxy_table[col])
        except:
            pass

    with open(galaxy_path, 'w') as f:
        cols = list(data.keys())
        f.write(','.join(cols) + '\n')
        for i in range(len(galaxy_table)):
            row = []
            for col in cols:
                val = data[col][i]
                try:
                    if np.isnan(val):
                        row.append('')
                    else:
                        row.append(str(val))
                except (TypeError, ValueError):
                    row.append(str(val) if val else '')
            f.write(','.join(row) + '\n')

    print(f"\nVoid galaxy table saved to: {galaxy_path}")

    # Save table info
    info_path = output_dir / 'void_tables_info.json'
    with open(info_path, 'w') as f:
        json.dump(tables_info, f, indent=2)

    print(f"Table info saved to: {info_path}")

    return output_dir


def main():
    """Main execution."""
    results = {
        'session': 84,
        'title': 'Void Catalog Download',
        'date': datetime.now().isoformat(),
        'status': 'STARTED'
    }

    # Download catalog
    catalogs = download_void_catalog()

    if catalogs is None:
        results['status'] = 'FAILED'
        results['error'] = 'Could not download catalog'
        print("\n" + "=" * 70)
        print("SESSION #84 VOID DOWNLOAD FAILED")
        print("=" * 70)
        return results

    # Analyze tables
    tables_info = analyze_void_tables(catalogs)
    results['n_tables'] = len(catalogs)

    # Find galaxy table
    galaxy_table, table_idx = find_galaxy_table(catalogs)
    results['galaxy_table_idx'] = table_idx
    results['n_galaxies'] = len(galaxy_table)
    results['columns'] = list(galaxy_table.colnames)

    # Save
    output_dir = save_void_catalog(galaxy_table, tables_info)
    results['output_dir'] = str(output_dir)
    results['status'] = 'SUCCESS'

    # Summary
    print("\n" + "=" * 70)
    print("VOID CATALOG DOWNLOAD COMPLETE")
    print("=" * 70)
    print(f"""
    Downloaded {len(catalogs)} tables from Douglass+ 2023

    Galaxy table: {len(galaxy_table)} galaxies
    Columns: {galaxy_table.colnames}

    Files saved to: {output_dir}

    NEXT: Cross-match with ALFALFA catalog
    """)

    # Save results
    results_path = Path(__file__).parent / 'results' / 'session84_void_download.json'
    results_path.parent.mkdir(exist_ok=True)
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    return results


if __name__ == '__main__':
    main()
