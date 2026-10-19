"""Self-audit of recalled measurements in CBP-Claude's substrate arc (Phases 1-18, 2026-06-22..24), for the 10-19
RECALL-code prediction. Reads each exploration AS FIRST COMMITTED. Lists lines carrying a measured value (±) alongside
a year-style citation or a named experiment, and flags whether a locator is present. Mechanical inventory, not coding."""
import os, re, subprocess, json
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
files = sorted(f for f in os.listdir(os.path.join(REPO, "explorations")) if re.match(r"2026-06-2[2-4]-phase", f))
LOC = re.compile(r"arXiv[: ]?\d|doi|https?://|\b(PRL|PRD|Phys\.? ?Rev|ApJ|MNRAS|Nature|Science|A&A|JCAP)\b[^\n]{0,20}\d", re.I)
MEAS = re.compile("±|\\\\pm|\\+/-|[<≤≲]\\s*~?\\s*[0-9]|×\\s*10|10[⁻−-]\\s*[0-9⁰¹²³⁴⁵⁶⁷⁸⁹]|10\\^\\{?-|[0-9]\\s*[eE]-[0-9]")
CITE = re.compile(r"\b(19|20)\d\d\b|Hulse|Taylor|Cassini|Planck|LIGO|Fermi|GRB|Gravity Probe|Eöt-?Wash|MICROSCOPE|lunar laser|LLR|Pound|Rebka|Hughes|Drever|Michelson|bound|measured|observed|constraint|limit|experiment|data\b|Weisberg|SME|Kostelecky|Kostelecký|Will\b|GW170817|Bertotti")
rows = []
for f in files:
    path = "explorations/" + f
    add = subprocess.run(["git", "-C", REPO, "log", "--diff-filter=A", "--format=%h", "--", path], capture_output=True, text=True).stdout.split()
    if not add: continue
    text = subprocess.run(["git", "-C", REPO, "show", f"{add[-1]}:{path}"], capture_output=True, text=True).stdout
    for n, line in enumerate(text.splitlines(), 1):
        if MEAS.search(line) and CITE.search(line):
            rows.append(dict(file=f, first_commit=add[-1], line=n, locator=bool(LOC.search(line)), text=line.strip()[:240]))
for r in rows:
    print(f"{'LOC ' if r['locator'] else 'NONE'} {r['file'][11:40]:30s} L{r['line']:<4d} {r['text'][:150]}")
print(f"\n{len(rows)} measurement-with-citation lines; {sum(not r['locator'] for r in rows)} without a locator")
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "inventory.json"), "w"), indent=1)
