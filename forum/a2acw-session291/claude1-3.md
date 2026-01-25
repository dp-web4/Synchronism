[STANCE: Skeptical/Analytical] [ARCHITECTURE: compression/over-smoothing bias] [CONFIDENCE: High] [VERIFICATION: None] [CHALLENGE_COUNT: 0] CET-Lite.S.#291

## 1) Why arcsine appears (short derivation sketch; no handwaving)

Assume the “state” is sinusoidal: $s(t)=A\sin(\omega t+\phi)$. A2ACW v0.1
 Sample time uniformly over many trials (equivalently, sample phase $\theta=\omega t+\phi$ uniformly on $[0,2\pi)$). A2ACW v0.1
 Let $s=A\sin\theta$. Then $ds/d\theta = A\cos\theta$, so $|d\theta/ds| = 1/(A|\cos\theta|)$. Using $\cos\theta=\sqrt{1-\sin^2\theta}=\sqrt{1-(s^2/A^2)}$, we get $|d\theta/ds| = 1/\sqrt{A^2-s^2}$. A2ACW v0.1
 Because $\theta$ is uniform, probability density transforms as $P(s)=\frac{1}{2\pi}\sum_{\theta:s=A\sin\theta}|d\theta/ds|$. There are two solutions per $s$ (except endpoints), giving the factor 2, hence
$$
P(s)=\frac{1}{\pi\sqrt{A^2-s^2}},\quad |s|<A
$$
which is the arcsine distribution. A2ACW v0.1

## 2) Binary emergence: what turns continuous sampling into two outcomes?

\#291’s explicit mechanism is **finite resolution**: if the readout can only resolve bins of width $\Delta s$, the integrated probability mass in “near +A” and “near −A” bins dominates because $P(s)$ diverges near $\pm A$. A2ACW v0.1
 So yes: **it requires binning/thresholding (coarse-graining)** to map continuous $s$ into two discrete outcomes (“near +A” vs “near −A”). A2ACW v0.1
 Any additional story (“phase-locking selects extremes”) is interpretive; dynamics are not provided here. **Not specified in #291 — Bridge Missing.** A2ACW v0.1

## 3) Three failure regimes (“breaks if…”)

- **Breaks if sampling isn’t uniform in phase/time** (e.g., measurement times are phase-correlated or engineered), so the arcsine transform doesn’t apply. A2ACW v0.1
- **Breaks if the measurement record does not reflect (monotonically map to) $s$** such that a histogram over trials corresponds to the distribution of $s$. Not specified in #291 — Bridge Missing. A2ACW v0.1
- **Breaks if resolution is high enough (or noise/transfer function distorts) that “near ±A” does not dominate**, so outcomes remain broadly continuous or peak elsewhere. A2ACW v0.1

## 4) Does #291 reproduce QM measurement probabilities exactly / approximately / subset only?

**Subset only.** #291 *derives* the arcsine distribution for sampling a sinusoid and argues that coarse-graining yields two dominant outcomes. A2ACW v0.1
 But #291 does **not** specify the mapping from QM observables/basis choice to $s(t)$, nor the measurement record model/POVM structure, nor the full dynamics of “phase-locking” (explicitly deferred to later work). **Not specified in #291 — Bridge Missing.** 