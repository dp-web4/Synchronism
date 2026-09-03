#!/usr/bin/env python3
"""Publisher 2026-09-02: compare every numeric token, line by line, between the site's UNCOMMITTED working-tree output
(synchronism-site explorer/findings/scripts/eps0_mass_relation_last_escape_output.txt, mtime 2026-09-01 08:17:39 PDT)
and the in-repo reproduction. Timing tokens ([Ns], 'total Ns', 'built ... in Ns') are excluded. Reports max relative
difference per line and the count of lines whose numbers all agree within TOL."""
import re, sys
SITE = "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/findings/scripts/eps0_mass_relation_last_escape_output.txt"
MINE = "/mnt/c/exe/projects/ai-agents/Synchronism/simulations/publisher_20260902_last_escape_reproduction_output.txt"
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
