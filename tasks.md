# Mini–Kelvin–Helmholtz (KH) Study – Tasks & Tracking

**Project:** Density- & Magnetic-Stabilisation of 2-D Kelvin–Helmholtz Billows (Mini Study)  
**PRD:** [PRD_MiniKH_Project.md](./PRD_MiniKH_Project.md)  
**Implementation Plan:** [IMPLEMENTATION_PLAN_MiniKH_Study.md](./IMPLEMENTATION_PLAN_MiniKH_Study.md)  
**Complexity Level:** 3 (Comprehensive Planning and Implementation)

---

## Meta Tasks

- [ ] **M1**: Ensure all team members have read and understood the PRD and Implementation Plan.
- [ ] **M2**: Set up a shared communication channel for the project (e.g., Slack channel, regular brief meetings).
- [ ] **M3**: Ensure version control (Git) is set up, with a dedicated feature branch for this project (e.g., `feature/mini-kh-study`).

---

## Phase 1: Repository & Environment Setup (Day 1)

- [ ] **1.1: Athena++ Build Confirmation**
    - **Description:** Verify that the current `main` or `develop` branch of Athena++ compiles successfully with the standard configuration required for this project (e.g., MHD, HDF5 output).
    - **Owner:** Sim Engineer
    - **Dependencies:** None
    - **Verification:**
        - `make clean && ./configure.py <options> && make -jN` completes without errors.
        - `bin/athena` executable is created.
    - **Status:** To Do

- [ ] **1.2: PRD Integration**
    - **Description:** Copy the `PRD_MiniKH_Project.md` to the repository root. Update the PI field to the designated project lead.
    - **Owner:** PI / Project Lead
    - **Dependencies:** None
    - **Verification:**
        - `PRD_MiniKH_Project.md` exists at repo root.
        - PI name is correctly updated in the document.
    - **Status:** To Do

- [ ] **1.3: Initial Commit of Planning Documents**
    - **Description:** Commit `IMPLEMENTATION_PLAN_MiniKH_Study.md` and this `tasks.md` to the feature branch.
    - **Owner:** Project Lead
    - **Dependencies:** M3
    - **Verification:** Files are committed and pushed to the remote feature branch.
    - **Status:** To Do

---

## Phase 2: C++ Problem Generator (`kh_bfield.cpp`) (Day 1)

- [ ] **2.1: Duplicate and Basic Setup**
    - **Description:** Duplicate `src/pgen/kh.cpp` to `src/pgen/kh_bfield.cpp`.
    - **Owner:** C++ Dev
    - **Dependencies:** 1.1
    - **Verification:**
        - `src/pgen/kh_bfield.cpp` exists and is a copy of `kh.cpp`.
    - **Status:** To Do

- [ ] **2.2: Add Magnetic Field Initialization**
    - **Description:** Modify `kh_bfield.cpp` to include uniform \(B_x\) initialization. This should be controlled by a new input parameter `problem/bx0`. The relevant problem setup is `iprob=2` (Frank et al. 1996 tanh profile).
    - **Owner:** C++ Dev
    - **Dependencies:** 2.1
    - **Details:**
        - Add `Real bx0 = pin->GetReal("problem","bx0");` in `MeshBlock::ProblemGenerator`.
        - Ensure `MAGNETIC_FIELDS_ENABLED` is handled.
        - Set `pfield->b.x1f(k,j,i) = bx0;` for all cells.
        - Set `pfield->b.x2f` and `pfield->b.x3f` to 0.0.
        - If `NON_BAROTROPIC_EOS`, add magnetic energy `0.5*bx0*bx0` to `phydro->u(IEN,...)`.
        - Consider if other `iprob` options in the new file should be removed or clearly marked as not supporting the new `bx0` parameter for this specific study. For this study, focus modifications on `iprob=2`.
    - **Verification:**
        - **Numerical:** Code compiles. `bx0` can be set in the input file and is read by the problem generator.
        - **Physical (Unit Test Idea):** A minimal test input that initializes with `bx0` and dumps the initial state. Check if `cons(IB1)` or `prim(B1)` in HDF5 output matches `bx0`.
    - **Status:** To Do

- [ ] **2.3: Build System Integration**
    - **Description:** Register `kh_bfield.cpp` in the build system (e.g., `CMakeLists.txt` or Makefile system if applicable).
    - **Owner:** C++ Dev
    - **Dependencies:** 2.1
    - **Verification:** Athena++ compiles successfully with the new problem generator selectable.
    - **Status:** To Do

- [ ] **2.4: Code Style and Linting**
    - **Description:** Ensure `kh_bfield.cpp` adheres to Athena++ coding standards and passes `clang-tidy` or other linting checks.
    - **Owner:** C++ Dev
    - **Dependencies:** 2.2
    - **Verification:** Linting tools report no major issues.
    - **Status:** To Do

---

## Phase 3: Input Deck (`inputs/kh_local.in`) (Day 1)

- [ ] **3.1: Draft Input Deck**
    - **Description:** Create `inputs/kh_local.in` configured for the baseline `hydro_A033` case (iprob=2, \( \Delta V = 1 \), density ratio 2, \(B_x=0\)).
    - **Owner:** Sim Engineer
    - **Dependencies:** 2.3 (to know the problem generator name, e.g., `kh_bfield`)
    - **Details:**
        - `<job>`: `problem_id = kh_bfield` (or similar, after 2.3)
        - `<problem>`: `iprob = 2`, `vflow = 0.5` (so \( \Delta V = 1 \)), `drat = 2.0` (for density ratio, if `iprob=1` structure is borrowed, or ensure tanh profile uses this ratio), `amp = 0.01` (perturbation amplitude), `bx0 = 0.0`.
        - `<mesh>`: `nx1 = 256`, `nx2 = 128`, `nx3 = 1`. `x1min=0, x1max=1`, `x2min=-0.5, x2max=0.5`. Reflecting boundaries for x1, periodic for x2 (or vice-versa based on typical KH setup, check `kh.cpp` conventions). PRD uses Lx=Ly=1, implies x1max=1, x2max=1, typically centered. Assuming `x1` is shear direction, `x2` is gradient.
            *Correction*: PRD implies `x1max=1.0, x2max=1.0`. Usually, KH has periodic in X (flow direction) and outflow/reflecting in Y (gradient direction). Let's assume `x1` is flow (periodic), `x2` is gradient (outflow or reflecting). Let domain be `x1=[0,1], x2=[-0.5,0.5]` or `x1=[0,1], x2=[0,1]` if layers are at y=0.25, 0.75. PRD uses `Frank et al. 1996` (iprob=2 in `kh.cpp`) which has `tanh((pcoord->x2v(j))/a)` implying layer at `x2=0`. So `x2min=-0.5, x2max=0.5` is appropriate.
        - `<meshblock>`: `nx1 = 256`, `nx2 = 128`.
        - `<hydro>`: `gamma = 1.6666666666666667` (for ideal gas, unless specified), `cfl_number = 0.3`.
        - `<time>`: `tlim` calculated based on \( t_{KH} = 1/\sigma_{hydro} \). For `k=2` (implies fundamental mode wavelength \( \lambda = L_x/1 = 1 \), so \( k = 2\pi/ \lambda = 2\pi \)), \( A=0.33 \). \( \sigma = (k \Delta V / 2) \sqrt{1-A^2} \). Need to confirm \( k \). `k=2` in PRD usually means \( \lambda = L_x / 2 \), so \( k_{phys} = 2 \pi / (L_x/2) = 4\pi/L_x \). If \(L_x=1\), \(k_{phys}=4\pi\).
            The formula provided \( \sigma=\frac{k\,\Delta V}{2}\sqrt{1-A^2} \) where \(k=2\) is stated. This \(k\) might be a mode number, not wavenumber. If it's mode number 2 on domain length 1, \( \lambda = 1/2 \), so true wavenumber is \(2\pi/\lambda = 4\pi\).
            Let's use the PRD's "k=2" directly as \(k_{wavenumber}=2\). Then \(\sigma = (2 \times 1 / 2) \sqrt{1-0.33^2} = \sqrt{1-0.1089} = \sqrt{0.8911} \approx 0.944\). So \(t_{KH} = 1/0.944 \approx 1.059\). Stop time = \(2 t_{KH} \approx 2.118\).
        - `<output>`: HDF5 dumps at regular intervals.
    - **Verification:** Input deck is parsable by Athena++. Baseline simulation (`hydro_A033`) runs using `kh_bfield.cpp`.
    - **Status:** To Do

- [ ] **3.2: Parameterization for CLI Overrides**
    - **Description:** Ensure `kh_local.in` can take `rho_t` (top layer density, assuming bottom is 1 or vice-versa, or `drat` for density ratio) and `mag_field` (for `bx0`) as command-line parameters. The PRD uses `rho0` (likely density of one layer if the other is 1.0) and `mag_field`. Adapt to `problem/drat` and `problem/bx0`.
    - **Owner:** Sim Engineer
    - **Dependencies:** 3.1
    - **Verification:**
        - Running `./bin/athena -i inputs/kh_local.in problem/drat=10.0 problem/bx0=0.2` correctly overrides these parameters.
        - Check in the log output or initial HDF5 dump.
    - **Status:** To Do

---

## Phase 4: Automation Script (`scripts/run_suite.sh`) (Day 2)

- [ ] **4.1: Create Run Suite Script**
    - **Description:** Develop `scripts/run_suite.sh` as per PRD §4. It should iterate through the three simulation cases, launching Athena++ with appropriate parameters and logging output.
    - **Owner:** Sim Engineer
    - **Dependencies:** 3.2
    - **Details:**
        - Create `logs/` directory if it doesn't exist.
        - Use `declare -A params` as in PRD. Adapt parameter names to match `kh_local.in` (e.g. `problem/drat` instead of `rho0`, `problem/bx0` instead of `mag_field`).
    - **Verification:**
        - Script is executable.
        - Runs all three simulations: `hydro_A033`, `hydro_A082`, `mhd_A033`.
        - Output logs are created in `logs/` directory (e.g., `logs/hydro_A033.log`).
        - HDF5 files are generated for each run.
    - **Status:** To Do

- [ ] **4.2: Initial Test Runs**
    - **Description:** Perform short test runs of all 3 cases using the `run_suite.sh` script to ensure they start and produce initial outputs.
    - **Owner:** Sim Engineer
    - **Dependencies:** 4.1
    - **Verification:** Each run starts, produces a few HDF5 files, and completes without crashing prematurely. Check logs for obvious errors.
    - **Status:** To Do

---

## Phase 5: Python Analysis (`analysis/khi_analysis.py`) (Day 2)

- [ ] **5.1: Scaffold Analysis Package**
    - **Description:** Create the `analysis/` directory. Set up a basic Python package structure (e.g., `__init__.py`). Install necessary libraries (`h5py`, `numpy`, `scipy`, `matplotlib`, `pytest`, `flake8`). Create `requirements.txt` or use a `pyproject.toml`.
    - **Owner:** Python Dev
    - **Dependencies:** None
    - **Verification:**
        - `analysis/` directory and `__init__.py` exist.
        - `pip install -r requirements.txt` (or equivalent) works.
        - `flake8 analysis/` runs without errors on empty files.
    - **Status:** To Do

- [ ] **5.2: HDF5 Data Loading**
    - **Description:** Implement functions in `khi_analysis.py` to load relevant data (density `rho`, velocity `vx`, `vy`, magnetic field `Bx`, `By`, pressure/energy if needed) from Athena++ HDF5 files.
    - **Owner:** Python Dev
    - **Dependencies:** 5.1, 4.2 (needs sample HDF5 files)
    - **Verification:**
        - Can load data arrays from HDF5 for a given time snapshot.
        - Dimensions and basic statistics of loaded arrays are reasonable.
    - **Status:** To Do

- [ ] **5.3: Implement Diagnostic: Growth Rate (\(\sigma_{num}\))**
    - **Description:**
        - Calculate perturbation amplitude A(t). This often involves tracking the transverse kinetic energy or the amplitude of the \(v_y\) velocity component at the interface, or interface displacement. For `iprob=2`, the perturbation is on \(v_y\).
        - Fit \( \ln A(t) \) slope over the initial exponential growth phase (e.g., first 2 e-folds, or first \(t_{KH}\)).
        - Calculate theoretical growth rate \(\sigma_{th}\) based on PRD formulae for each run.
    - **Owner:** Python Dev / Physicist
    - **Dependencies:** 5.2
    - **Verification:**
        - **Numerical:** Function returns a numerical value for \(\sigma_{num}\) and \(\sigma_{th}\).
        - **Physical (hydro_A033):** For `hydro_A033`, \(\sigma_{num}/\sigma_{th} \in [0.85, 1.15]\).
        - Plot \( \ln A(t) \) vs \(t\) to visually inspect the linear growth phase.
    - **Status:** To Do

- [ ] **5.4: Implement Diagnostic: Mixing Layer Width (\(\Delta(t)\))**
    - **Description:**
        - Define and calculate mid-density contour width \(\Delta(t)\) (e.g., region where density is between \(\rho_t + 0.1(\rho_b - \rho_t)\) and \(\rho_b - 0.1(\rho_b - \rho_t)\), or based on vorticity layer).
        - Calculate the slope of \(\Delta(t)\) vs \(t\) in the self-similar mixing phase (later times).
    - **Owner:** Python Dev / Physicist
    - **Dependencies:** 5.2
    - **Verification:**
        - **Numerical:** Function returns slope of \(\Delta(t)\).
        - **Physical:** Slope for `hydro_A033` should be around \(s \approx 0.18 \times \Delta V\). Slopes for `hydro_A082` and `mhd_A033` should be less than `hydro_A033`.
        - Plot \(\Delta(t)\) vs \(t\).
    - **Status:** To Do

- [ ] **5.5: Implement Diagnostic: Wavenumber Peak (\(k(t)\))**
    - **Description:**
        - For a given time, take a 1D cut of vorticity \(\omega_z(x)\) or \(v_y(x)\) at \(y=0\) (interface).
        - Perform FFT to get power spectrum \(P(k_x)\).
        - Identify the peak wavenumber \(k_{peak}(t)\).
        - Track \(k_{peak}(t)\) over time to detect vortex pairing (e.g., \(k_{peak}\) halves).
    - **Owner:** Python Dev / Physicist
    - **Dependencies:** 5.2
    - **Verification:**
        - **Numerical:** Function returns \(k_{peak}(t)\).
        - **Physical:** `hydro_A033` should show evidence of \(k_{peak}\) decreasing (pairing). `mhd_A033` should show suppressed pairing compared to `hydro_A033`. `hydro_A082` might also show different pairing behavior.
        - Plot \(k_{peak}(t)\) vs \(t\), or plot power spectra at different times.
    - **Status:** To Do

- [ ] **5.6: Implement Diagnostic: Energy Conservation (\(\Delta E/E_0\))**
    - **Description:** Read total energy from Athena++ history (.hst) files for each run. Calculate percentage drift \( (E(t) - E_0)/E_0 \).
    - **Owner:** Python Dev
    - **Dependencies:** 5.1, 4.2 (needs .hst files)
    - **Verification:**
        - **Numerical:** Function returns max energy drift.
        - **Physical:** For all runs, \( |\Delta E / E_0| \le 1\% \).
    - **Status:** To Do

- [ ] **5.7: Output Generation: `results.csv` and Plots**
    - **Description:**
        - Script should output `results.csv` with columns: `run_id, sigma_num, sigma_theory, sigma_ratio, slope_Delta, paired?` (paired? can be boolean or notes on \(k_{peak}\) evolution), `max_energy_drift_percent`.
        - Generate `fig_sigma.png` (e.g., \( \ln A(t) \) plots for all runs, with fits).
        - Generate `fig_delta.png` (\(\Delta(t)\) plots for all runs).
        - Generate `fig_k_peak.png` (\(k_{peak}(t)\) plots or representative spectra).
        - (Optional) Add a composite plot showing 2D density/vorticity snapshots at key times for each run to illustrate KH billows, stabilization, and pairing.
    - **Owner:** Python Dev
    - **Dependencies:** 5.3, 5.4, 5.5, 5.6
    - **Verification:**
        - `results.csv` is created with correct data.
        - All specified PNG/PDF plots are generated and look publication-quality.
    - **Status:** To Do

- [ ] **5.8: Script CLI and Assertions**
    - **Description:** `khi_analysis.py` should be runnable as a script (e.g., `python analysis/khi_analysis.py --run_dir path/to/simulation_outputs`). It should perform assertions based on PRD §5 and exit with non-zero code if assertions fail.
    - **Owner:** Python Dev
    - **Dependencies:** 5.7
    - **Verification:**
        - Script runs via CLI.
        - Correctly passes for runs meeting criteria.
        - Exits with error code for runs failing criteria (test with dummy data if needed).
    - **Status:** To Do

- [ ] **5.9: Pytest Unit Tests**
    - **Description:** Write unit tests for key analysis functions using synthetic/mock data or a very small, known HDF5 output. Test data loading, core calculations (e.g., \(\sigma_{th}\) calculation, a simple FFT peak find).
    - **Owner:** Python Dev / Py QA
    - **Dependencies:** 5.2 - 5.6
    - **Verification:** `pytest analysis/` passes. Code coverage is reasonable.
    - **Status:** To Do

- [ ] **5.10: Flake8 Compliance**
    - **Description:** Ensure all Python code in `analysis/` is `flake8` clean.
    - **Owner:** Python Dev
    - **Dependencies:** 5.1 - 5.9
    - **Verification:** `flake8 analysis/` reports no errors.
    - **Status:** To Do

---

## Phase 6: Full Simulation Runs & Analysis (Day 2-3)

- [ ] **6.1: Execute Full Simulations**
    - **Description:** Run all three simulations to completion using `scripts/run_suite.sh`.
    - **Owner:** Sim Engineer
    - **Dependencies:** 4.1, (ensure 3.1 `tlim` is correct)
    - **Verification:** All simulations complete. Final HDF5 files and .hst files are present for all runs.
    - **Status:** To Do

- [ ] **6.2: Perform Full Analysis**
    - **Description:** Run `khi_analysis.py` on the completed simulation outputs.
    - **Owner:** Python Dev / Sim Engineer
    - **Dependencies:** 5.8, 6.1
    - **Verification:**
        - `results.csv` and all figures are generated.
        - Assertions pass for all runs as per PRD success gates:
            - \(\sigma_{num}/\sigma_{th} \in [0.85,1.15]\) (for hydro_A033, potentially for others where applicable).
            - \(\Delta E/E_0 \le 1\%\).
            - Pairing suppressed in `mhd_A033` run relative to `hydro_A033`.
            - Mixing layer width \(\Delta(t)\) slope drops from `hydro_A033` to `hydro_A082`/`mhd_A033`.
    - **Status:** To Do

- [ ] **6.3: Review and Iterate (if necessary)**
    - **Description:** Review analysis results. If assertions fail or results are unexpected, debug simulations or analysis scripts. This may involve adjusting simulation parameters, fixing bugs in `kh_bfield.cpp` or `khi_analysis.py`.
    - **Owner:** Entire Team
    - **Dependencies:** 6.2
    - **Verification:** Final results are consistent with physical expectations and PRD success criteria.
    - **Status:** To Do

---

## Phase 7: Continuous Integration (Optional but Recommended) (Day 3)

- [ ] **7.1: Setup GitHub Action Workflow**
    - **Description:** Create a GitHub Action workflow (`.github/workflows/ci.yml`).
    - **Workflow Steps:**
        1. Checkout code.
        2. Setup Python. Install dependencies.
        3. Build Athena++ (consider caching dependencies or `obj/` directory if build is too long, or use a pre-built Docker container with Athena++).
        4. Run `scripts/run_suite.sh` (possibly with very short `tlim` for CI speed, or use mock data).
        5. Run `analysis/khi_analysis.py`.
        6. (Optional) Upload `results/` directory as an artifact.
    - **Owner:** DevOps / CI Lead
    - **Dependencies:** 4.1, 5.8
    - **Verification:** CI pipeline runs successfully on push/PR to the feature branch. Job fails if `run_suite.sh` or `khi_analysis.py` (assertions) fail.
    - **Status:** To Do

---

## Phase 8: Documentation & Publishing (Day 3)

- [ ] **8.1: Draft Research Note**
    - **Description:** Draft the research note (`doc/mini_kh_note.md` or LaTeX) based on the `results.csv` and figures. Incorporate physics interpretations.
    - **Owner:** PI / Physicist
    - **Dependencies:** 6.2
    - **Format:** Markdown or LaTeX, ≤3 pages.
    - **Content:** Abstract, Introduction (KH, stabilization), Methods (Athena++, setup), Results (σ, Δ, k, E figs + discussion), Conclusion.
    - **Verification:** Draft is complete and includes all key results and figures.
    - **Status:** To Do

- [ ] **8.2: README Updates (Optional)**
    - **Description:** If CI is implemented, add badges to the main `README.md` for build status. Add a brief section to `README.md` describing the MiniKH study and linking to the PRD and results.
    - **Owner:** Project Lead
    - **Dependencies:** 7.1 (for badges)
    - **Verification:** README is updated.
    - **Status:** To Do

- [ ] **8.3: Final Code Review & Merge**
    - **Description:** Perform a final review of all new/modified code (`kh_bfield.cpp`, `kh_local.in`, `run_suite.sh`, `khi_analysis.py`, CI workflow). Merge feature branch into main/develop.
    - **Owner:** Entire Team
    - **Dependencies:** All preceding tasks.
    - **Verification:** Code review approved. Branch merged.
    - **Status:** To Do

- [ ] **8.4: Prepare Deliverables Package**
    - **Description:** Collate all deliverables: `results.csv`, three PNG/PDF plots, `analysis/` python code, and the research note.
    - **Owner:** Project Lead
    - **Dependencies:** 6.2, 8.1
    - **Verification:** All listed deliverables are present and organized.
    - **Status:** To Do

---

**End of Project** 