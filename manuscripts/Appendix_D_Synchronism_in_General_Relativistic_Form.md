## Appendix D: Synchronism in General Relativistic Form (Draft Framework)

This appendix outlines a preliminary framework for embedding Synchronism within a relativistic description of gravity. The goal is not to replace General Relativity (GR), but to understand GR as an **emergent, effective description** of deeper coherence dynamics in a single-observer Synchronism model.

We proceed in three layers:
1. Ontological stance (what is fundamental in Synchronism),
2. Mapping to an effective GR-like description,
3. Worked toy limits that connect to the Newtonian Synchronism equations used in the main text.

This appendix is intentionally minimal and will be expanded as the relativistic Synchronism program develops.

---

### D.1 Ontological Stance: Single Observer and Coherence

Synchronism is a CFD(Computational Fluid Dynamics)-like representational model of the universe that takes the following view of “spacetime” and “observers”:

1. **Single observer (the grid)**  
   There is one true observer: the universe, represented as a discrete grid of cells updated in discrete ticks. All physical systems (galaxies, clocks, detectors, etc.) are **patterns** evolving on this grid.

2. **Time as update rate**  
   What we call “time” is the progression of grid ticks. Local “time dilation” is not time itself changing but a change in the **pattern evolution rate**: how many grid updates are available for a pattern’s internal dynamics versus its motion relative to other patterns.

3. **Coherence as a physical resource**  
   Coherence is a measure of how tightly the degrees of freedom of a pattern are constrained to evolve together. It is encoded in a dimensionless function \(C(\rho) \in (0,1]\) of local baryonic density \(\rho\). High coherence means many constraints and slower effective evolution; low coherence means fewer constraints and a greater effective gravitational response.

4. **Geometry as emergent description**  
   GR’s curved spacetime is treated as an **emergent macroscopic description** inferred by patterns that attempt to fit their observations using continuous fields and coordinates. In Synchronism, coherence gradients are fundamental; curvature is a mathematically efficient way to describe how patterns constrained by those gradients behave.

Thus, Synchronism does not start from a fundamental metric \(g_{\mu\nu}\). Instead, it starts from coherence \(C(\rho)\), grid dynamics, and a Markov Relevancy Horizon (MRH), and then asks: *What effective metric would an internal pattern reconstruct if it modeled its experience with GR?*

---

### D.2 Effective Newtonian Limit: Modified Poisson Equation

The main text adopts a phenomenological Newtonian limit in which Synchronism modifies the Poisson equation via an effective gravitational coupling:

\[
\nabla^2 \Phi = 4\pi G\,\frac{\rho}{C(\rho)} ,
\]

with
\[
C(\rho) = \tanh\!\left[\gamma \,\ln\!\left(\frac{\rho}{\rho_{\mathrm{crit}}} + 1\right)\right], \quad \gamma = 2,
\]
and
\[
\rho_{\mathrm{crit}} = A\,V_{\mathrm{flat}}^{B},
\]
where \(A\) and \(B\) are derived from Jeans stability and galaxy size–velocity scaling.

In this limit:

- The **effective gravitational constant** is
  \[
  G_{\mathrm{eff}}(\rho) = \frac{G}{C(\rho)}.
  \]
- The observed acceleration and circular velocity are
  \[
  g_{\mathrm{obs}} = \frac{g_{\mathrm{bar}}}{C(\rho)}, \qquad
  V_{\mathrm{obs}}(r) = \frac{V_{\mathrm{bar}}(r)}{\sqrt{C(\rho(r))}}.
  \]

This is the non-relativistic behavior that any relativistic Synchronism–GR formulation must reproduce in the weak-field, slow-motion limit.

> **[CORRECTION 2026-08-09] The two equations above are not the same law, and only the second one carries any of this program's evidence.**
>
> §D.2 presents \(\nabla^2\Phi = 4\pi G\,\rho/C\) (call it **L1**) and \(g_{\mathrm{obs}} = g_{\mathrm{bar}}/C\) (**L3**) as one statement, the second following "in this limit" from the first. It does not follow. \(g = g_{\mathrm{bar}}/C\) is the spherical Gauss solution of a *different* field equation,
> \[
> \nabla\cdot\!\left[C(\rho)\,\nabla\Phi\right] = 4\pi G\rho \qquad (\textbf{L2}),
> \]
> since integrating L2 over a ball gives \(C\,g\,r^2 = G M_{\mathrm{bar}}(<r)\) directly. L1 instead makes \(\rho/C\) the *source*, so \(g(r) = (G/r^2)\!\int_0^r 4\pi r'^2 \rho/C\,dr'\). The two coincide **iff \(C\) is spatially constant** — and \(C\) is this framework's entire content, so no choice of \(\gamma\) or \(\rho_{\mathrm{crit}}\) reconciles them.
>
> **Measured** (`simulations/publisher_20260809_appendixD_coupling_fork.py`, in-repo SPARC mass models, Lelli et al. 2016, 113 galaxies at Q ≤ 2 and \(i>30°\), self-consistent spherical construction \(\rho_{\mathrm{sph}} = (4\pi r^2)^{-1}\,d[V_{\mathrm{bar}}^2 r/G]/dr\) so that no scale-height convention enters): at the outermost measured point \(\log_{10}(g_{L3}/g_{L1})\) has **median +0.821 dex**, IQR [+0.515, +1.416], range [+0.023, +3.203]; **105/113 galaxies exceed 0.3 dex and 47/113 exceed 1 dex**. The constant-\(C\) control returns \(2.6\times10^{-15}\) dex, confirming the fork *is* the spatial variation of \(C\).
>
> The gap is **\(\gamma\)-invariant to 0.0000 dex** between \(\gamma = 2\) and \(\gamma = 0.489\) — but the reason bounds how far that fact travels. Every SPARC point sits at \(\rho/\rho_{\mathrm{crit}} \le 0.045\) (median \(6.1\times10^{-6}\)), deep in \(C\)'s linear regime where \(C \to \gamma\rho/\rho_{\mathrm{crit}}\) to within 2.4%; there \(1/C \propto 1/\gamma\) in *both* laws and \(\gamma\) cancels identically. The invariance holds **at this \(\rho_{\mathrm{crit}}\)** and is a consequence of the calibration being far enough off that no galaxy reaches \(C\)'s knee — the same fact as the already-recorded 2–5 dex miss of the density law against rotation curves.
>
> **Provenance.** L1 has exactly one implementation in this repository — `simulations/session72_spherical_toy_model.py`, written 2025-12-01 to complete §D.6 — and **zero citations in the whitepaper**. L3 appears in 26 simulation scripts and carries every number the whitepaper quotes. L2 has **zero implementations**; it was written down as prose in whitepaper §5.15 on 2026-08-04, by an author who did not know §D.2 existed. So the framework's only *stated* field equation has never been evaluated against data, and the force law that carries all of its evidence had no field equation implemented for it.
>
> **The eliminated reading is the better-fitting one.** Against the observed accelerations at the same radii, **L1 is closer than L2 = L3 in 113/113 galaxies** at both \(\gamma = 2\) and \(\gamma = 0.489\) (median offsets +4.64 vs +5.65 dex, and +5.25 vs +6.26 dex). Both are catastrophically wrong; the elimination below is **structural, not a goodness-of-fit verdict**. *(These absolute offsets use the spherical \(\rho_{\mathrm{sph}}\) construction and are not comparable to the whitepaper's "2–5 orders of magnitude", which uses \(\rho = \Sigma/2h\).)*
>
> **L1 is eliminated a priori, and the reason generalizes past it.** As \(\rho \to 0\), \(C \to \gamma\rho/\rho_{\mathrm{crit}}\), so \(\rho_{\mathrm{eff}} = \rho/C \to \rho_{\mathrm{crit}}/\gamma\), a **constant** (verified to \(<10^{-4}\) relative for \(\gamma \in [0.25, 5]\) down to \(\rho/\rho_{\mathrm{crit}} = 10^{-10}\)). L1's source never vanishes in vacuum; \(M_{\mathrm{eff}}(<r)\) diverges as \(r^3\). The sharper objection, which the divergence obscures: \(\rho_{\mathrm{crit}} = A V_{\mathrm{flat}}^B\) is a **per-galaxy** constant, so L1 assigns the same point of empty space a different source density for every galaxy one asks about — the equation has no single-valued source. That is ill-posedness, not divergence, and **it survives into L2 = L3**, which the divergence argument does not. Contracting the same limit with \(u_\mu u_\nu\) for dust carries it into §D.3: \(\mathcal{T}_{\mu\nu} \to (\rho_{\mathrm{crit}}/\gamma)\,u_\mu u_\nu\) in vacuum, a nonzero stress–energy filling all space whose direction requires a preferred vacuum 4-velocity.
>
> **Recommended resolution, stated but not applied here** (this is a physics decision, not an editorial one): §D.2 and §D.6.1 should adopt L2, whose spherical solution is *already* §D.2's own second line, and drop \(\rho/C\) as a Poisson source. Nothing downstream in this program depends on L1. Registered at `Research/proposals/appendix_D_states_two_force_laws_and_only_one_carries_evidence_20260809.md`.

---

### D.3 Emergent Einstein-Like Equations

To connect with GR, we introduce an **effective field equation** that treats curvature as a coarse-grained description of coherence-constrained pattern evolution. The simplest Ansatz consistent with the modified Poisson equation is:

\[
G_{\mu\nu}[g] = 8\pi G \,\mathcal{T}_{\mu\nu},
\]

where \(G_{\mu\nu}\) is the Einstein tensor derived from an effective metric \(g_{\mu\nu}\), and \(\mathcal{T}_{\mu\nu}\) is an **effective stress–energy tensor** defined by

\[
\mathcal{T}_{\mu\nu} \equiv \frac{T_{\mu\nu}}{C(\rho)},
\]
with \(T_{\mu\nu}\) the usual baryonic stress–energy tensor and \(\rho\) the local baryonic density measured in an appropriate frame.

This can be interpreted in two equivalent ways:

1. **Density-dependent coupling**
   \[
   G_{\mu\nu} = 8\pi G_{\mathrm{eff}}(\rho)\,T_{\mu\nu}, \quad G_{\mathrm{eff}}(\rho) = \frac{G}{C(\rho)}.
   \]

2. **Coherence-weighted source**
   \[
   G_{\mu\nu} = 8\pi G\,\frac{T_{\mu\nu}}{C(\rho)}.
   \]

In the weak-field, static limit,
- the \(00\)-component reduces to
  \[
  \nabla^2 \Phi \approx 4\pi G\,\frac{\rho}{C(\rho)},
  \]
- matching the Synchronism-modified Poisson equation.

**Important:** This is an **effective** equation. The fundamental dynamics live on the grid and are encoded in coherence \(C(\rho)\) and MRH, not in a fundamental metric. The metric and Einstein tensor are emergent objects that patterns infer when they coarse-grain over many grid cells.

---

### D.4 Light, Gravitational Waves, and the Speed Limit

In Synchronism, the universal speed \(c\) is understood as the **maximum pattern propagation rate permitted by the grid**. Patterns with minimal internal coherence overhead can approach this throughput limit; patterns with complex internal structure and high coherence cannot.

- **Photons and gravitational waves** are modeled as **simple patterns** (minimal internal degrees of freedom) that propagate at the grid’s throughput limit \(c_{\rm grid}\), which internal witnesses infer as “the speed of light.” 
- **Massive composite objects** have deep internal coherence hierarchies and must spend a fraction of each grid update maintaining internal alignment, reducing their effective maximum speed below \(c_{\rm grid}\).

This picture explains why:

- Light and gravitational waves travel at essentially the same speed,
- Macroscopic objects cannot be accelerated arbitrarily close to \(c\),
- Relativistic time dilation can be understood as a shift in internal vs. external update allocation rather than “time changing.”

In an effective GR language, this is encoded in the usual null geodesics of the metric \(g_{\mu\nu}\):

- Both photons and gravitational waves follow (or approximate) null geodesics of the emergent metric.
- Synchronism demands that any relativistic extension respect observational constraints such as
  \[
  \left|\frac{c_{\rm GW} - c_{\rm EM}}{c}\right| \ll 10^{-15},
  \]
  as seen in events like GW170817.

---

### D.5 Gravity as Coherence Gradient and Emergent Geodesics

In the Synchronism ontology, gravity is not a fundamental force but a manifestation of **coherence gradients**:

- Regions of higher baryonic density tend to exhibit higher coherence.
- Coherence gradients translate into **constraints** on how patterns can evolve and move.
- Patterns “fall” along directions that reduce their coherence maintenance cost, which, in the macroscopic limit, look like trajectories toward massive objects.

We can sketch an effective connection between this picture and geodesics:

1. Define an effective “coherence potential” \(U(\mathbf{x})\), for example
   \[
   U(\mathbf{x}) \propto -\ln C(\rho(\mathbf{x})),
   \]
   so that regions of low coherence (large mass enhancement) correspond to lower potential.

2. Consider an effective action for a pattern’s center-of-mass worldline,
   \[
   S_{\mathrm{eff}} = \int \left[ -m\,\sqrt{-g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu} - \lambda\,U(\mathbf{x}) \right] d\tau,
   \]
   where the first term describes motion in the emergent metric and the second term encodes coherence gradients, with \(\lambda\) a coupling constant to be calibrated.

3. In the limit where the coherence perturbation can be absorbed into a redefinition of the metric components (e.g. through a conformal factor or potential-like \(g_{tt}\) perturbation), the equations of motion reduce to **geodesics** of an effective metric \(g_{\mu\nu}[C]\).

In this way, Synchronism’s “coherence gradient” picture and GR’s “objects follow geodesics” picture become two languages describing the same effective behavior of patterns under coherence constraints.

---

### D.6 Static, Spherically Symmetric Toy Model

As a concrete testbed, consider a static, spherically symmetric baryonic mass distribution \(\rho(r)\) with total baryonic mass \(M_{\mathrm{bar}}\). Synchronism prescribes a radial coherence profile \(C(\rho(r))\) and a modified effective density
\[
\rho_{\mathrm{eff}}(r) = \frac{\rho(r)}{C(\rho(r))}.
\]

#### D.6.1 Newtonian limit

In the Newtonian Synchronism limit:
\[
\frac{1}{r^2}\frac{d}{dr}\left[r^2 \frac{d\Phi}{dr}\right] = 4\pi G\,\rho_{\mathrm{eff}}(r)
= 4\pi G\,\frac{\rho(r)}{C(\rho(r))},
\]
with \(\Phi(r)\) the gravitational potential. The circular velocity is then
\[
V_{\mathrm{obs}}^2(r) = r\,\frac{d\Phi}{dr}.
\]

#### D.6.2 Effective metric Ansatz

To connect with GR, we can consider an effective, static, spherically symmetric metric of the form
\[
ds^2 = -e^{2\phi(r)} c^2 dt^2 + e^{2\lambda(r)} dr^2 + r^2 d\Omega^2,
\]
where \(\phi(r)\) and \(\lambda(r)\) are functions to be determined.

The standard GR equations relate \(\phi(r)\) and \(\lambda(r)\) to \(\rho(r)\) and pressure \(p(r)\). In Synchronism, we instead postulate that the **effective** source entering the Einstein tensor is \(\rho_{\mathrm{eff}}(r)\), leading to equations structurally similar to the Tolman–Oppenheimer–Volkoff system but with \(\rho \rightarrow \rho_{\mathrm{eff}}\).

In the weak-field limit (\(|\phi| \ll 1\), \(|\lambda| \ll 1\)), this reduces to
\[
\phi(r) \approx \frac{\Phi(r)}{c^2}, \qquad \nabla^2 \Phi = 4\pi G\,\rho_{\mathrm{eff}}(r),
\]
thus recovering the modified Poisson equation. Autonomous derivations can refine this Ansatz, derive explicit forms of \(\phi(r)\) and \(\lambda(r)\) for given \(C(\rho(r))\), and compare with the Schwarzschild solution in the appropriate limits.

---

### D.7 Summary and Future Directions

This appendix sketches how Synchronism can be embedded in an effective GR-like framework:

- **Fundamental level:** single observer (grid), coherence \(C(\rho)\), MRH, and pattern dynamics.
- **Newtonian effective level:** modified Poisson equation
  \[
  \nabla^2 \Phi = 4\pi G\,\frac{\rho}{C(\rho)}.
  \]
- **Relativistic effective level:** Einstein-like equations with coherence-weighted stress–energy
  \[
  G_{\mu\nu} = 8\pi G\,\frac{T_{\mu\nu}}{C(\rho)},
  \]
  recovering the Newtonian limit above and allowing photons and gravitational waves to propagate as null patterns in an emergent metric.

Open tasks for future work include:

1. Constructing a covariant definition of \(\rho\) and \(C(\rho)\), possibly via scalar invariants built from \(T_{\mu\nu}\).
2. Ensuring energy–momentum conservation, potentially via additional internal fields representing coherence dynamics.
3. Exploring cosmological solutions (Friedmann-like equations) with coherence-dependent coupling.
4. Deriving detailed predictions for gravitational lensing, gravitational waves, and strong-field tests (e.g., binary pulsars, black hole shadows) in Synchronism.

These steps will turn the sketch presented here into a fully specified Synchronism–GR framework capable of being confronted with the full range of relativistic gravitational data.