#!/usr/bin/env python3
"""Publisher 2026-09-03: compare every numeric token, line by line, between the site's COMMITTED output
(synchronism-site explorer/findings/scripts/eps0_mass_relation_last_escape_output.txt, committed 2026-09-02 in
c4b8022) and this repo's complete re-run (publisher_20260903_last_escape_reproduction_output.txt, 313 s).

This supersedes publisher_20260902_last_escape_compare.py, which could not run: the 09-02 pass launched the
reproduction as a background job, yielded its turn to wait on it, and died at exit=0 with 11 of 27 eps0 solves
written and no analysis. Timing tokens are excluded; the rng is seeded (20260901) so permutation p-values and
bootstrap CIs must match EXACTLY, not just within tolerance.

RESULT 2026-09-03: 108 lines vs 108, 77/77 numeric lines exact, 0 text diffs.

CAVEAT THAT THIS SCRIPT CANNOT DETECT, AND THE REASON IT IS WRITTEN DOWN HERE. A verbatim reproduction
re-runs the same reader over the same cache. Line 41 of both outputs reads 'rho_c at top edge 59%'; that
number is VOID -- column 1 of epsilon0_per_galaxy_fw.npy is per-galaxy chi2, not rho_c (writer:
epsilon0_free_the_ceiling_rescue.py:220, `np.vstack([e0s, best_per_gal, MOND_F, NN])`; col1 median 49.4,
max 1.17e4, where a rho_c on that grid is 1e-6..1e-2). Both runs read col1, so both print it, so the
comparison scores it EXACT. A reproduction gate answers 'did the arithmetic run the same way'. It is blind
by construction to 'is the input what it is labelled', because the defect is upstream of the arithmetic.
"""
import re, sys
SITE = "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/findings/scripts/eps0_mass_relation_last_escape_output.txt"
MINE = "publisher_20260903_last_escape_reproduction_output.txt"
TOL = 0.02
num = re.compile(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?')
def clean(l):
    l = re.sub(r'\[\d+s\]', '', l); l = re.sub(r'in \d+s', '', l); l = re.sub(r'^total \d+s', 'total', l)
    return l
a = [clean(l.rstrip()) for l in open(SITE, encoding='utf-8', errors='replace')]
b = [clean(l.rstrip()) for l in open(MINE, encoding='utf-8', errors='replace')]
print(f"site lines {len(a)}  mine lines {len(b)}")
n_num = agree = exact = 0; worst = []
for i, (x, y) in enumerate(zip(a, b)):
    nx, ny = num.findall(x), num.findall(y)
    if not nx and not ny:
        if x.strip() != y.strip(): print(f"TEXT DIFF line {i+1}: {x!r} vs {y!r}")
        continue
    if len(nx) != len(ny):
        print(f"TOKEN COUNT DIFF line {i+1}: {x!r} vs {y!r}"); continue
    n_num += 1
    md = 0.0
    for p, q in zip(nx, ny):
        fp, fq = float(p), float(q)
        d = abs(fp - fq) / max(abs(fp), abs(fq), 1e-12) if (fp != 0 or fq != 0) else 0.0
        md = max(md, d)
    if md == 0: exact += 1
    if md <= TOL: agree += 1
    else: worst.append((md, i+1, x.strip(), y.strip()))
print(f"numeric lines compared {n_num}: exact {exact}, within {TOL:.0%} {agree}, outside {len(worst)}")
for md, i, x, y in sorted(worst, reverse=True): print(f"  line {i} maxrel {md:.3f}\n    site: {x}\n    mine: {y}")
