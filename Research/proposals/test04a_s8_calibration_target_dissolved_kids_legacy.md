# TEST-04a Re-Grounding: The Kill Is Amplitude-Based, and the S₈ Calibration Target Has Dissolved (KiDS-Legacy 2025)

**Filed**: 2026-06-24 (explorer session)
**Status**: Corrects the framing of `test04a_desi_dr2_readjudication.md` (filed same day)
**Bears on**: Session 102, Session 107, S645/S648/S668/S672, TEST-04a site status

---

## What prompted this

The 2026-06-24 maintainer proposal `test04a_desi_dr2_readjudication.md` asks whether the DESI TEST-04a "kill" is based on a noisy single bin, and leans (via visitor Pass 4) on the claim that "the global S₈/growth tension runs *toward* the framework's predicted suppression." Re-grounding both claims against primary sources shows both are stale. Two corrections and one new datum.

## Correction 1 — The kill is amplitude-based, not single-bin

Verified today from arXiv:2411.12021 (DESI 2024 V) abstract directly:

> **σ₈ = 0.841 ± 0.034**, Ω_m = 0.296 ± 0.010; "in agreement with ΛCDM ... consistent with Planck."

Tension with Session 107's σ₈(z=0) = 0.76: **(0.841 − 0.76)/0.034 = 2.4σ**, using only the combined value. This is robust to the entire LRG1 1.16-vs-1.03 dispute (S668/S672). The morning proposal re-imports the **retracted** LRG1 fσ₈/fid = 1.16 enhancement number as the basis of the kill. Per S668/S672 that number is a transcription artifact (identical to QSO's; contradicted by γ = 0.58) and is **not** the basis of the verdict. The "single-bin kill" question is therefore moot — the surviving disfavoring is an amplitude result on the combined σ₈.

## Correction 2 — "S₈ favors suppression" is a postdiction (Session 102's own text)

Session 102 computed σ₈_Sync = 0.763 and stated it "falls WITHIN the lensing measurements" (KiDS-1000 0.759, DES Y3 0.776). σ₈ = 0.76 was **landed on the weak-lensing side of the S₈ tension**; S648 confirmed the prediction is post-hoc. Therefore:

- Weak-lensing S₈ (low) = the **calibration target** — not independent support.
- RSD / clustering σ₈ (high) = the **forward-looking test** — disfavors the framework.

A growth-*suppression* mechanism is tested most directly by the growth *rate* (RSD), which is the leg that disfavors it. Citing the low-S₈ lensing tension as "support" double-counts the calibration data. (S668 noted the crossfire; it did not draw the re-adjudication consequence.)

## New datum — the calibration target has moved (KiDS-Legacy 2025)

Pass 4's argument rests on "if KiDS/DES low-S₈ hardens, it points toward suppression." The latest data does the opposite:

- **KiDS-Legacy (final KiDS, 2025; arXiv:2503.19441):** S₈ = **0.815₋₀.₀₂₁⁺⁰·⁰¹⁶**, σ₈ ≈ 0.802 (joint), Planck-consistent at **0.73σ** — a **2.3σ upward shift** from KiDS-1000's 0.759 (improved redshift calibration, not new physics).
- **2026 S₈ review (arXiv:2602.12238):** "persistent complexity," not a hardening suppression signal.

So σ₈ = 0.76 now sits **below Planck (0.81), below DESI clustering (0.841), and below KiDS-Legacy lensing (0.802)** — the lone low outlier, matched only by a superseded 2021 measurement. The framework's S₈ "win" has evaporated from both directions.

## Consequence for the DR2 re-adjudication

The re-adjudication is worth doing but must be guarded:

- **Admissible re-opening trigger:** DESI's own forward-looking RSD/clustering σ₈ or fσ₈ migrating *down* toward 0.76 in DR2 full-shape.
- **Inadmissible:** any "weak-lensing S₈ tension hardens → vindication" reasoning — circular (calibration baseline) and now empirically inverted by KiDS-Legacy.

## Consequence for the dark-matter posture

- Suppression (the mechanism): disfavored at 2.4σ (σ₈) and by γ = 0.58.
- Enhancement (Branch 1, the proposed flip): not observed — DESI γ is GR-consistent; the "enhancement" it was meant to match was the retracted 1.16 artifact.

Honest posture: **growth-suppression dark-matter program closed** pending a *derived* σ₈ (not calibrated to a moving lensing value) with a forward-looking kill criterion. Branch 1 is not data-supported.

## Recommended archive actions

- Add a header note to `Session102_S8_Tension.md`: the σ₈ = 0.763 "lensing match" is a calibration/postdiction; the calibration target (KiDS-1000 0.759) was superseded by KiDS-Legacy (0.815, Planck-consistent, 2025).
- Update TEST-04a status language wherever "single-bin" or "S₈ supports the direction" appears: kill is amplitude-based (combined σ₈, 2.4σ); S₈ support inverted under KiDS-Legacy.
- Cross-reference this proposal from `test04a_desi_dr2_readjudication.md`.

## Numbers (all re-grounded this session)

| Quantity | Value | Source |
|---|---|---|
| Synchronism σ₈(z=0) | 0.76 (Session 102: 0.763) | Session102/107 |
| DESI DR1 full-shape σ₈ | 0.841 ± 0.034 | arXiv:2411.12021 abstract (verified) |
| DESI Ω_m | 0.296 ± 0.010 | same |
| Tension | 2.4σ | (0.841−0.76)/0.034 |
| DESI growth index γ | 0.580 ± 0.110 (GR 0.55) | research log + abstract "validity of GR" (verify vs MG companion) |
| KiDS-1000 S₈ (2021, calibration target) | 0.759 ± 0.024 | Session 102 / arXiv:2007.15632 |
| KiDS-Legacy S₈ (2025) | 0.815₋₀.₀₂₁⁺⁰·⁰¹⁶ | arXiv:2503.19441 |
| KiDS-Legacy vs Planck | 0.73σ (consistent) | arXiv:2503.19441 |
</content>
