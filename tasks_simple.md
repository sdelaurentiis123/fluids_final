# KHI Metrics Suite – Simplified Task List

This checklist tracks the minimal steps needed to extract, analyse, and visualise key Kelvin–Helmholtz instability diagnostics from existing Athena++ `.athdf` outputs using the `analysis/khi_analysis.py` toolkit.

---

## Phase 0 – Environment
- [ ] **0.1 Create Python venv** – `python -m venv .venv && source .venv/bin/activate`
- [ ] **0.2 Install requirements** – `pip install -r analysis/requirements.txt`
- [ ] **0.3 Sanity test imports** – `python - <<<'import h5py, numpy, scipy, matplotlib, pandas'` exits 0

## Phase 1 – Data Discovery & Preparation
- [ ] **1.1 Locate ATHDF snapshots** – verify expected path(s) in `output/*/*.athdf`
- [ ] **1.2 (Optional) Symlink run directory** to `data/` for cleaner CLI usage

## Phase 2 – Metric Extraction CLI
- [ ] **2.1 Add CLI entry‐point** – expose `python -m analysis.khi_analysis --run_dir <PATH>`
- [ ] **2.2 Hook existing metric functions** (`calculate_metric_M1_growth_rate`, `...M2`, `...M3`, `...M4`)
- [ ] **2.3 Write `metrics.csv`** with time-series and scalar outputs

## Phase 3 – Plots
- [ ] **3.1 Plot summary panel** (`summary_plot.png`) combining KE, growth fit, vorticity, mixing variance
- [ ] **3.2 Plot FFT peak** (`fft_vy_peak_lambda.png`)

## Phase 4 – Automation Helpers
- [ ] **4.1 `run_all.sh`** – iterate over multiple run folders listed in `config.json`
- [ ] **4.2 Aggregate comparison plot** – cross-run overlay of growth rates & mixing metrics

## Phase 5 – Validation
- [ ] **5.1 Unit test fixtures** in `analysis/tests/` (mock `.athdf` arrays)
- [ ] **5.2 CI job** – `pytest` + run one tiny `.athdf` sample to ensure pipelines stay green

---

Keep PRs focused – one phase at a time.  Each completed phase should deliver runnable code and pass CI. 