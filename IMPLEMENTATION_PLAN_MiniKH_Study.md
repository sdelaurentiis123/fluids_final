## Mini–Kelvin–Helmholtz (KH) Study – Implementation Plan

*Document owner — <YOUR NAME>, Staff Research Scientist / Staff SWE*

---

### 0  Table of Contents  
1. Executive Summary  
2. Scope & Non-Goals  
3. Existing Code Audit  
4. Work Breakdown Structure (WBS)  
5. Detailed Implementation Tasks  
6. Verification & Validation Strategy  
7. Deliverables & Directory Layout  
8. Timeline & Resource Allocation  
9. Risk Register & Mitigations  
10. Appendix A – Physics Reference  

---

### 1  Executive Summary
We will conduct three 2-D Athena++ simulations to quantify the stabilising roles of
1) density contrast and 2) in-plane magnetic tension on the Kelvin–Helmholtz (KH) instability.
The project produces: 3 HDF5 dumps, a reproducible analysis pipeline, three publication-quality
figures, and a ≤3-page research note.  CI gates will fail when diagnostics deviate
from theory by >15 % or if energy drifts by >1 %.

---

### 2  Scope & Non-Goals
**In-scope**  
• Add a *minimal* problem generator `kh_bfield.cpp` (≈20 LOC diff to `kh.cpp`).  
• Prepare input deck `inputs/kh_local.in` parametrisable via CLI.  
• Bash runner `scripts/run_suite.sh` (3 jobs).  
• Python analysis `analysis/khi_analysis.py` with CLI + pytest coverage.  
• GitHub Action (optional) that executes the suite & asserts physics tests.  
• Documentation (`PRD_MiniKH_Project.md`, this plan, and research note template).  

**Out-of-scope**  
• GPU optimisation; 3-D simulations; AMR studies; full MHD energy budgets.  

---

### 3  Existing Code Audit
Grep confirms the stock Kelvin–Helmholtz problem generator already lives at
```src/pgen/kh.cpp``` (538 LOC).  Relevant excerpt:
```6:6:src/pgen/kh.cpp
//! \brief Problem generator for KH instability.
```
No magnetic-field toggle currently exists for iprob=2.  No dedicated analysis package is
present in the repository (`analysis/` is absent).

---

### 4  Work Breakdown Structure (WBS)
1. **Repo-level housekeeping**  
   1.1 Check that `bin/athena` builds with present options.  
   1.2 Copy PRD to root & update PI field.  
2. **C++: Problem Generator**  
   2.1 Duplicate `src/pgen/kh.cpp → kh_bfield.cpp`.  
   2.2 Inject uniform \(B_x\) initialisation controlled by new parameter `problem/bx0`.  
   2.3 Register file in `CMakeLists.txt` / build system.  
3. **Inputs**  
   3.1 Draft `inputs/kh_local.in` (baseline iprob=2 tanh shear).  
   3.2 Expose runtime overrides (`rho0`, `mag_field`).  
4. **Automation Scripts**  
   4.1 `scripts/run_suite.sh` (iterates over assoc-array matrix).  
5. **Python Analysis**  
   5.1 Scaffold `analysis/` package (pip-installable, flake8 clean).  
   5.2 Implement diagnostics: growth rate fit, layer-width slope, FFT peak, energy drift.  
   5.3 Produce `results.csv`, plots, exit-code assertions.  
   5.4 Unit tests using synthetic data fixtures.  
6. **Continuous Integration**  
   6.1 GH Action: build Athena++, run suite, run analysis, cache artifacts.  
7. **Documentation & Publishing**  
   7.1 Auto-generate README badges from CI.  
   7.2 Draft research note skeleton (`doc/mini_kh_note.md`).  

---

### 5  Detailed Implementation Tasks
| # | Deliverable | Owner | Est. hrs | Blocking | Notes |
|---|-------------|-------|----------|----------|-------|
|2.1|`kh_bfield.cpp` created|C++ dev|1|—|Copy kh.cpp; remove unused iprobs; add BX setter|
|2.2|Parameter plumbed|C++ dev|0.5|2.1|`pin->GetReal("problem","bx0")`|
|2.3|Build rules updated|C++ dev|0.5|2.1|Makefile / CMake|
|3.1|Input deck authoring|Sim engineer|1|2.1|iprob=2, vflow=0.5|
|4.1|Runner script|Sim eng|0.5|3.1|Associative array per PRD|
|5.1|Package scaffold|Py dev|1|—|`analysis/__init__.py`|
|5.2|Diagnostics implemented|Py dev|3|5.1|Use `h5py`, `numpy`, `scipy`, `matplotlib`|
|5.3|Plots + CSV|Py dev|1|5.2|Save under `results/`|
|5.4|Pytests|Py QA|1|5.2|Synthetic fixtures|
|6.1|CI pipeline|DevOps|1|4.1,5.2|Self-hosted or GH Runners|
|7.2|Research note|PI|2|5.3|Markdown + equations|

_Total ≈ 12 hrs concentrated effort._

---

### 6  Verification & Validation Strategy
1. **Code correctness**:  
   • `clang-tidy` (C++) & `flake8` + `pytest` (Python).  
2. **Physical correctness**:  
   • Numerical growth rate σ\_num must satisfy |σ/σ\_th – 1| ≤ 15 %.  
   • Energy drift ≤ 1 %.  
   • Vortex pairing only in hydro
a033 run.  
   CI fails if assertions throw.  
3. **Reproducibility**: all outputs seeded (`ISEED=1234`), runner script prints Athena++ commit hash.

---

### 7  Deliverables & Directory Layout
```
PRD_MiniKH_Project.md
IMPLEMENTATION_PLAN_MiniKH_Study.md   ← this file
doc/mini_kh_note.md
src/pgen/kh_bfield.cpp
inputs/kh_local.in
scripts/run_suite.sh
analysis/
  __init__.py
  khi_analysis.py
  tests/
results/  (CI artifact)
```

---

### 8  Timeline & Resource Allocation (3 × 1-day sprints)
| Day | Sprint focus | Milestones |
|-----|--------------|------------|
|1|C++ + input deck|Problem generator compiles · baseline run passes|
|2|Automation & Python|Suite runs · diagnostics validated on baseline|
|3|CI + docs|CI green · research note draft committed|

---

### 9  Risk Register & Mitigations
| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
|Incorrect B-field initialisation|M|H|Compare to analytical Alfvén speed; unit test grid cell Bx|
|FFT peak aliasing|L|M|Zero-pad & cross-check with analytic k|2D FFT|
|CI walltime|M|M|Use low resolution (256×128) & O2 build|
|Athena++ build quirks on GH runners|M|L|Pre-install deps; cache `obj/`|

---

### 10  Appendix A – Physics Reference
See §2 of PRD for growth-rate and mixing-layer formulae. 