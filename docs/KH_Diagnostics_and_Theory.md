# Diagnostics & Theory Linkage for the Mini–KHI Project

This file spells out **what we measure, how we measure it, and which analytic formula each metric should satisfy**.  It doubles as the spec for the assertions inside `analysis/khi_analysis.py`.

---

## 1  Metric Catalogue

| ID | Symbol                  | Physical meaning                          | Extraction recipe                                                                                                                  | Units             |
| -- | ----------------------- | ----------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| M1 | $\sigma_{\mathrm{num}}$ | Linear e‑fold rate of the primary KH mode | • Stripe sample of $v_y$ at the interface (4‑cell band)<br>• FFT → amplitude $A(t)$<br>• Fit slope of $\ln A$ over first 2 e‑folds | 1/$t$             |
| M2 | $\Delta(t)$             | Shear‑layer thickness                     | Domain‑average density ρ(x,y,t).  Find y where ρ crosses mid‑value $(ρ_b+ρ_t)/2$; width between top & bottom crossings.            | fraction of $L_y$ |
| M3 | $k_{\max}(t)$           | Dominant billow wavenumber                | Span‑average vorticity ω(x) → 1‑D FFT; index of max magnitude (exclude k=0).                                                       | dimensionless     |
| M4 | $ΔE/E_0$                | Total‑energy drift                        | History file: kinetic + internal (+magnetic).  Relative change between t=0 and t=t\_end.                                           | dimensionless     |

---

## 2  Analytic Benchmarks

| Metric                                 | Theory                                                                               | Notes                                                                                   |
| -------------------------------------- | ------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------- |
| **Growth rate** $\sigma_{\mathrm{th}}$ | $\displaystyle \sigma = \frac{k\,ΔV}{2}\sqrt{1-A^2}$                                 | For hydro runs with Mach ≲ 0.5 (compressibility correction < 15 %).                     |
| **Magnetic reduction**                 | $\displaystyle \frac{\sigma}{\sigma_{0}} = \sqrt{1-\left(\tfrac{2v_A}{ΔV}\right)^2}$ | Use baseline hydro $\sigma_0$ and measured $v_A = B_x/\sqrt{ρ_b+ρ_t}$.                  |
| **Mixing‑layer slope** $s$             | $\Delta(t) = s\,ΔV\,t$ with $s\approx0.18$ for low Atwood                            | Expect $s\downarrow$ with $A↑$ and non‑zero $v_A$.  No hard assertion, but trend check. |
| **Pairing**                            | First billow‑merger time $t_p \approx 3/\sigma$ when $k=2→1$                         | High‑A or MHD run should *not* pair before t\_end.                                      |

---

## 3  Assertion Matrix (Automated Tests)

| Test ID | Expression                                         | Pass threshold                                      |    |                                 |
| ------- | -------------------------------------------------- | --------------------------------------------------- | -- | ------------------------------- |
| T1      | (                                                  | \sigma\_{\text{num}}/\sigma\_{\text{th}} - 1        | )  | < **0.15** for hydro A=0.33 run |
| T2      | Hydro high‑A run: $\sigma_{\text{num}}$ / baseline | < **0.6** (≥ 40 % suppression)                      |    |                                 |
| T3      | MHD run: $\sigma_{\text{num}}/\sigma_{0}$          | Matches magnetic reduction formula within **±0.15** |    |                                 |
| T4      | Energy drift                                       | (                                                   | ΔE | /E\_0 < 0.01)                   |
| T5      | Pairing event                                      | Detected only in hydro A=0.33 run; not in others    |    |                                 |

These tests are encoded in `analysis/tests/test_diagnostics.py` and are executed by `pytest` in CI.

---

## 4  Theory Refresher (one‑paragraph explanations)

### 4.1  Atwood‑controlled growth

Linear KH growth arises from pressure continuity across a perturbed vortex sheet.  Solving Laplace's equation for the velocity potential and matching kinematic & dynamic conditions yields
$\sigma = \tfrac12 kΔV \sqrt{1-A^2}\,,$  where $A$ reduces the pressure perturbation that can accelerate the heavier fluid.  As $A\to1$ the sheet is marginally stable.

### 4.2  Magnetic‑tension stabilisation

Insert a uniform $B_x$ and repeat the jump calculation; magnetic pressure adds a restoring term $k^2v_A^2$ which subtracts from $σ^2$.  Modes satisfy $ΔV>2v_A$ to remain unstable; otherwise they become travelling waves.

### 4.3  Mixing‑layer entrainment

After non‑linear roll‑up the layer grows self‑similarly: entrainment velocity $E = d\Delta/dt$ is proportional to $ΔV$.  Experiments give coefficient $s≈0.18$ for small $A$; larger density or magnetic tension starves small‑scale vortex shedding, lowering $s$.

### 4.4  Vortex pairing criterion

In 2‑D, billows merge when their self‑induced velocity offsets half the background shear; empirically $t_p\approx3/σ$ for initial mode number $k=2$.  Any mechanism that slows growth (high A, tension) delays this cascade.

---

## 5  How to expand (future work)

\* Add viscosity to measure Re‑dependent cut‑off.
\* Switch to 3‑D to capture rib vortices and sub‑harmonic pairing paths.
\* Compare compressible cases (M≈2) to low‑Mach benchmark.

--- 