# The site's chemistry table was never regenerated from this archive; the DE local-fluid bound is one-sided in kind

*Site maintainer, 2026-11-02. Sources: site explorer finding 2026-10-31
(`synchronism-site/explorer/findings/chemistry-table-rows-are-not-traceable-and-the-z-null-fails-on-the-real-rows.md`,
PREREG site `6a1c4b6`) and the 2026-11-02 simulated-visitor log (graduate and researcher personas).
No bucket moves. Count 6, range 3–8, Bucket 0 = 0.*

## 1. Chemistry: the record was wrong in both directions, and the verdict holds on firmer ground

The site's 23-row "chemistry correlation" table dates from the site's initial commit (2026-02-21) and was never
regenerated from this archive. The explorer traced it row by row:

- The rows are **23 properties**, each a correlation over a different material set from a different session
  (n = 19–33). There is no 23-material sample.
- In the 16 high-r rows, **no cited session is about the property**, and several r values belong to other
  properties (e.g. "thermal conductivity 0.961" is redox E°, #72; "electronegativity 0.979" is chemical softness, #118).
- **A documented failure was shown as a success for eight months:** electrical conductivity appeared as +0.955
  "Coherence helps"; session #81 records σ vs γ_phonon at **r = −0.414, "WRONG SIGN"**.
- The polynomial-in-Z null, run for the first time on the 10 recoverable property rows, **fails** (|r(P, Z)| ≤ 0.53,
  quadratic R ≤ 0.51, against framework r 0.60–0.98). The 2026-05-10 "a Z-polynomial matches within 0.07" is refuted
  on data, not only mis-described (the 10-31 block already withdrew it as a synthetic-curve run).
- The high r values are identities: γ_phonon = 2T/θ_D, so r(P, 1/γ_ph) ≡ r(P, θ_D); and γ_optical = 13.6 eV / IE
  (`electron_affinity_coherence.py:80`), so "χ vs 1/γ_optical" is Pauling χ vs ionization energy (0.94–0.96), a relation
  Mulliken's definition builds in.
- The elastic rows rest on diamond: jackknife-min r falls to 0.60 (E), 0.63 (G), 0.37 (B).

**Ledger consequence.** The Bucket 3 row "γ ~ 1 universal boundary → Debye θ_D" is right but incomplete: the optical
channel reparametrizes ionization energy, not θ_D. Suggested row text: *"Reparametrizes θ_D (phonon channel,
γ_ph = 2T/θ_D, an identity) and ionization energy (optical channel, γ_opt = IE_H/IE, an identity); ~86% of the '89%'
restate the Debye model."* And the chemistry failures listed anywhere in the archive should include the conductivity
wrong-sign result (#81), which the compilation layer dropped.

**Process point (record-facing).** This is the fourth instance of a compiled number standing for months against its
primary (07-10, 10-09 Desmond, 10-10 DESI QSO row, now). This one is larger: a whole table. The site has now added a
per-row *Trace* column (traced / corrected / untraced) and keeps untraced rows visible rather than deleting them.
Recommendation: any table the archive or site quotes should carry a pointer to the script that produced each number.

## 2. DE local-fluid horn: the γ < ½ side is gradient-unstable, not amplitude-bounded

The 09-15 block states the local-fluid exclusion as |ε| ≲ 10⁻⁵ (ε = 2γ − 1). The explorer's own run already shows the
two sides differ in kind (c_s² ≈ 0.223 ε): for γ < ½, c_s² < 0, and in linear theory every scale below
~aH/√|c_s²| grows exponentially, so on the adiabatic reading **no γ < ½ is stable at any distance from ½**; 10⁻⁵ only
marks where the blow-up reaches P(k) scales. For γ > ½ it is Jeans suppression and 10⁻⁵ is an amplitude bound.
The DESI mean-density fit γ = 0.487 sits on the unstable side. (A graduate persona also placed the SPARC 0.489 there;
that γ belongs to C_g, keyed on acceleration, so it is not on this axis — the 09-24 (ii) conflation again.)
The escape is unchanged and two-sided: non-adiabatic c_s² (= 0, or an entropy-carrying second field). Verdict
unchanged; the statement should be one-sided in kind.

## 3. For the project lead (all gate on you)

1. **The count rule as worded excludes one of the six.** The site's inclusion rule says "executed against a registered
   criterion on data the framework did not fit"; row #6 (Bell/CHSH) is a simulation, and it sits in the rule-invariant
   floor of 3. Today the site says so beside the rule and beside the landing-page headline, and adds in plain words:
   *of the three that survive every rule, one (TEST-10) is an external-data refutation of content the framework itself
   asserts.* This is the 10-07 "5 on data + 1 construction check" recommendation and the 10-29 floor sentence, now forced
   by a contradiction in the printed rule. Please rule on it.
2. **Two results an outside researcher called "worth writing up on their own".** The researcher persona singled out
   (a) the LLR Ġ/G bound on any G_eff = G·f(ρ_s), |d ln G_eff/d ln ρ_s| ≲ 4×10⁻⁹–4×10⁻⁸ (09-24 / 10-27), and (b) the
   γ·A unidentifiability of density-keyed fits at SPARC's x ~ 7×10⁻⁵ (10-23). Both are framework-independent. It asked
   for two things before either is citable: the ρ_s-gradient model shown explicitly, and why Earth's and the Sun's
   environments see the same f in the LLR mapping. This bears on the 10-27 ranking criterion (which negative leads a
   preprint).
3. **Frame.** Today's top visitor items were again bookkeeping and ordering; the record-integrity finding of the week
   (a whole table never regenerated) came from the explorer reading primaries, not from the visitor loop. That fits the
   10-24 frame note: the visitor loop now polishes the record, and the record's deepest defects sit where no persona looks.
   An explorer sweep of every initial-commit data table on the site is seeded.
