# Markov Chains as a Formal Lens on Coherence and Emergence

**Status:** exploration seed — translation first, physics claims later (if ever).

Markov-chain mathematics looks unusually well matched to several Synchronism concepts, especially **entity**, **coherence**, **emergence**, and the **Markov Relevancy Horizon (MRH)**.

The important caution is ontological: using a Markov chain does **not** require the Synchronism substrate itself to be fundamentally stochastic. A deterministic tick-by-tick system can yield an effective Markov description after coarse-graining, incomplete knowledge, or projection onto a smaller state space. The Markov machinery is therefore best treated first as a **witness-level mathematical description**.

## 1. Tick dynamics → transition dynamics

For a Markov process,

\[
P(X_{t+1}=x'\mid X_t=x,X_{t-1},\ldots)=P(X_{t+1}=x'\mid X_t=x).
\]

At a chosen MRH, a distribution over states evolves as

\[
p_{t+1}=p_t P.
\]

This makes the transition matrix/operator \(P\) a candidate mathematical object for describing how intent-pattern classes transform from tick to tick at that level of description.

## 2. Coherence as metastability — but with a stronger criterion

A first approximation to a coherent pattern is a region \(C\) of state space that tends to keep the dynamics inside itself:

\[
P(X_{t+1}\in C\mid X_t\in C)\approx 1.
\]

But persistence alone is too weak. A frozen state would score highly without exhibiting meaningful emergence.

A stronger definition is **time-scale separation**:

- microstates inside \(C\) mix/rearrange relatively quickly;
- escape from \(C\) is relatively slow.

Define

\[
\mathcal C(C)=\frac{\tau_{\rm exit}(C)}{\tau_{\rm mix}(C)}.
\]

Large \(\mathcal C\) means the detailed microscopic configuration can change many times while the macroscopic identity persists. That looks much closer to Synchronism's definition of an entity as a recurring pattern over sequences of ticks.

For a substochastic transition matrix restricted to \(C\), its dominant eigenvalue \(\lambda_C<1\) gives a natural persistence scale. Near metastability,

\[
\tau_{\rm exit}\sim\frac{1}{1-\lambda_C}.
\]

This suggests a route from the qualitative phrase **stable intent pattern** to something measurable.

## 3. Emergence as successful coarse-graining

Let many microstates be grouped into macrostate classes:

\[
\{x_1,x_2,\ldots\}\rightarrow\{A,B,C,\ldots\}.
\]

If the projected macrostate process is itself approximately Markov, the coarse description has gained a degree of **predictive closure**. In the exact case this is Markov-chain **lumpability**.

A strong-lumpability condition for a partition \(\{C_a\}\) is that any two microstates in the same block have the same total transition probability into every target block:

\[
\sum_{k\in C_b}P_{ik}=\sum_{k\in C_b}P_{jk}
\quad\forall i,j\in C_a,\;\forall b.
\]

Then the macro-dynamics can be written without remembering which microstate inside the macrostate was occupied.

That gives a potentially useful operational statement:

> **Emergence occurs when a collection of lower-level states admits a coarse-graining that becomes a stable predictive state-space in its own right.**

Approximate lumpability is probably more relevant than exact lumpability for physical systems.

## 4. Markov blanket, MRH, and conditional independence

A Markov blanket \(B\) around internal variables \(I\) separates them statistically from external variables \(E\):

\[
I\perp E\mid B,
\]

or equivalently

\[
I(I;E\mid B)=0
\]

in conditional-mutual-information language.

This is close to, but not identical with, Synchronism's MRH. The current MRH definition is the boundary beyond which correlations become negligible **for a pattern's dynamics**.

That suggests an explicitly dynamical, approximate version:

\[
I(I_{t+1};E_t\mid I_t,B_t)\le \epsilon.
\]

The MRH would then be the smallest practical boundary for which information outside the boundary contributes less than tolerance \(\epsilon\) to prediction of the pattern's next state.

This may be a useful formal distinction:

- **Markov blanket:** the statistical conditional-independence boundary;
- **MRH:** the scale/context at which residual dependence becomes irrelevant for the question being asked.

In other words, MRH may behave like an **approximate, task-relative Markov blanket** rather than simply being another name for one.

## 5. A possible hierarchy

Putting the pieces together:

\[
\text{micro transition dynamics}
\rightarrow
\text{metastable state clusters}
\rightarrow
\text{approximately lumpable macro dynamics}
\rightarrow
\text{new MRH / effective blanket}
\rightarrow
\text{repeat}.
\]

This supplies a mathematical candidate for Synchronism's fractal hierarchy: each level need not merely be visually or conceptually coarser; it must earn its existence by exhibiting persistence and predictive closure over a characteristic range of timescales.

There is established mathematics nearby: metastable Markov processes, lumpability and chain reduction, Markov Stability across timescales, conditional independence / Markov blankets, and causal-emergence measures. The useful question is therefore **not whether Synchronism can rename these ideas**, but whether putting them together under the Synchronism ontology yields a sharper operational definition, a useful cross-domain transfer, or eventually a discriminator.

A full exploration arc is proposed in:

[`explorations/2026-08-17-markov-coherence-emergence-arc-proposal.md`](../explorations/2026-08-17-markov-coherence-emergence-arc-proposal.md)

The first toy calculation is in:

[`explorations/2026-08-17-markov-phase1-metastability-and-lumpability.md`](../explorations/2026-08-17-markov-phase1-metastability-and-lumpability.md)

with a reproducible script at:

[`simulations/markov_coherence_toy.py`](../simulations/markov_coherence_toy.py)

## References / adjacent formalisms

- Delvenne, Yaliraki & Barahona (2010), *Stability of graph communities across time scales*, PNAS 107:12755–12760. DOI: 10.1073/pnas.0903215107
- Hoel, Albantakis & Tononi (2013), *Quantifying causal emergence shows that macro can beat micro*, PNAS 110:19790–19795. DOI: 10.1073/pnas.1314922110
- Geiger & Temmel (2014), *Lumpings of Markov chains, entropy rate preservation, and higher-order lumpability*, Journal of Applied Probability / arXiv:1212.4375
- Friston et al. (2020), *Markov blankets, information geometry and stochastic thermodynamics*, Proceedings of the Royal Society A 476:20190806.

These are comparison tools and mathematical neighbors, not evidence for Synchronism.
