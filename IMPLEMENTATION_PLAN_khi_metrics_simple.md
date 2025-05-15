## KHI Metrics Suite – Concise Implementation Plan

*Document Owner — <<YOUR NAME>>, Staff SWE / Computational Physicist*

---

### 1. Executive Summary
The goal is to provide an out-of-the-box Python pipeline (≈300 LOC) that ingests existing Athena++ `.athdf` snapshots and emits a tidy `metrics.csv` plus a handful of publication-quality plots summarising Kelvin–Helmholtz growth, mixing, and vortex dynamics.  The plan assumes simulations have already been run; **no C++ or simulation work is required.**  Where possible we reuse the proven functions already present in `analysis/khi_analysis.py`.

### 2. Scope & Non-Goals
In-scope:
• CLI tool `analysis/khi_analysis.py` callable via `python -m analysis.khi_analysis --run_dir <PATH>`
• Minimal config file (`config.json`) for batch processing multiple runs
• Helper shell script `run_all.sh` for batch execution
• Unit tests + GitHub CI workflow to keep pipeline green
• README usage examples and sample plot

Out-of-scope:
• Re-running simulations, generating new `.athdf`s
• Extending Athena++ source code
• 3-D or AMR support (2-D, uniform grid only)

### 3. Existing Code Audit
`analysis/khi_analysis.py` already implements:  
• M1 growth-rate fit (`calculate_metric_M1_growth_rate`)  
• M2 shear-layer thickness  
• M3 dominant wavenumber  
• M4 energy drift  

Each routine takes a directory of `.athdf` files and optional plotting directory.  The code paths are battle-tested for both 2-D and 3-D snapshots.

Gaps:
1. No single CLI wrapper that orchestrates all metrics for a given run
2. No batch runner
3. No `requirements.txt`
4. No unit tests / CI
5. No aggregation plot across multiple runs

### 4. Work Breakdown Structure (WBS)
1. **Environment**  
   1.1 Author `analysis/requirements.txt`  
   1.2 Add import sanity test in CI
2. **CLI Wrapper**  
   2.1 Add `if __name__ == "__main__":` block to `analysis/khi_analysis.py`  
   2.2 Parse `--run_dir`, `--plot_dir`, `--json_out`  
   2.3 Call existing metric functions and dump CSV / PNG
3. **Batch Runner**  
   3.1 Define `config.json` (see PRD §6)  
   3.2 Write `scripts/run_all.sh` that loops over runs, calling the CLI
4. **Plots**  
   4.1 Re-use per-metric plots produced by existing routines  
   4.2 Add `summary_plot.png` that stitches them together via `matplotlib.gridspec`
5. **Validation**  
   5.1 Fabricate mini `.athdf` fixtures with synthetic arrays (100 × 50 grid)  
   5.2 Pytest: L2-norm difference between metric outputs and known truth < 1e-8  
   5.3 Flake8 compliance
6. **CI**  
   6.1 GitHub Action: set up Python 3.11, install deps, run pytest, execute CLI on fixture data, upload plots as job artefacts

### 5. Deliverables & Directory Layout
```
analysis/
  khi_analysis.py      (extended with CLI)
  requirements.txt
  tests/
    test_metrics.py
scripts/
  run_all.sh           (wrapper around CLI + config)
config.json            (list of run directories & labels)
results/               (auto-created; holds CSV & PNG for each run)
IMPLEMENTATION_PLAN_khi_metrics_simple.md  ← this file
tasks_simple.md
```

### 6. Verification & Acceptance
• Unit tests green on GH CI  
• `python -m analysis.khi_analysis --run_dir ./output/kh_hydro_A033_256x128` exits 0 and writes `metrics.csv` + plots  
• For provided reference run (`kh_hydro_A033_256x128`) growth-rate fit σ_num ≈ 1.0 ± 0.15  
• All plots render without Matplotlib warnings

### 7. Timeline (1-Day Sprint)
Morning: Environment + CLI (Tasks 0-2)  
Afternoon: Batch runner, plots, unit tests, CI (Tasks 3-6)

---
Ready for execution.  Refer to `tasks_simple.md` for actionable check-boxes. 