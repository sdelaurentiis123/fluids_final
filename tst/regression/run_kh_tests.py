import subprocess
import os
import shutil
import json # Added for saving params

# Define the Athena++ executable and input file
ATHENA_EXECUTABLE = "./bin/athena"
BASE_INPUT_FILE = "inputs/kh_local.in"
OUTPUT_BASE_DIR = "output"

# Define test cases
test_cases = [
    {
        'id': 'hydro_A033',
        'drat': 2.0,
        'bx0': 0.0,
        'by0': 0.0,
        'bz0': 0.0,
        'resolution': {'nx1': 256, 'nx2': 128, 'nx3': 1}
    },
    {
        'id': 'mhd_B0.1_res128x64',
        'drat': 2.0,
        'bx0': 0.1,
        'by0': 0.05,
        'bz0': 0.0,
        'resolution': {'nx1': 128, 'nx2': 64, 'nx3': 1}
    },
    {
        'id': 'mhd_B0.1_res64x32_3D',
        'drat': 1.0,
        'bx0': 0.1,
        'by0': 0.0,
        'bz0': 0.01,
        'resolution': {'nx1': 64, 'nx2': 32, 'nx3': 32}
    },
    {
        'id': 'hydro_A033_512x256',
        'drat': 2.0,
        'bx0': 0.0,
        'by0': 0.0,
        'bz0': 0.0,
        'resolution': {'nx1': 512, 'nx2': 256, 'nx3': 1}
    },
    {
        'id': 'mhd_B0.1_hires_2D_256x128',
        'drat': 1.0,
        'bx0': 0.1,
        'by0': 0.0,
        'bz0': 0.0,
        'resolution': {'nx1': 256, 'nx2': 128, 'nx3': 1}
    },
    {
        'id': 'mhd_B0.1_3D_hires_128x64x64',
        'drat': 1.0,
        'bx0': 0.1,
        'by0': 0.0,
        'bz0': 0.01,
        'resolution': {'nx1': 128, 'nx2': 64, 'nx3': 64}
    },
    # {
    #     'id': 'hydro_A033_1024x1024',
    #     'drat': 2.0,
    #     'bx0': 0.0,
    #     'by0': 0.0,
    #     'bz0': 0.0,
    #     'resolution': {'nx1': 1024, 'nx2': 1024, 'nx3': 1}
    # },
    # {
    #     'id': 'mhd_B0.1_1024x1024_2D',
    #     'drat': 2.0,
    #     'bx0': 0.1,
    #     'by0': 0.05,
    #     'bz0': 0.0,
    #     'resolution': {'nx1': 1024, 'nx2': 1024, 'nx3': 1}
    # },
    # {
    #     'id': 'mhd_B0.1_512x512x32_3D',
    #     'drat': 1.0,
    #     'bx0': 0.1,
    #     'by0': 0.0,
    #     'bz0': 0.01,
    #     'resolution': {'nx1': 512, 'nx2': 512, 'nx3': 32}
    # },
]

def run_test(test_params):
    '''Runs a single Athena++ test case.'''
    test_id = test_params['id']
    drat_val = test_params.get('drat')
    bx0_val = test_params.get('bx0')
    by0_val = test_params.get('by0', 0.0)
    bz0_val = test_params.get('bz0', 0.0)
    resolution = test_params.get('resolution')

    res_str_parts = []
    if resolution and 'nx1' in resolution:
        res_str_parts.append(str(resolution['nx1']))
    if resolution and 'nx2' in resolution:
        res_str_parts.append(str(resolution['nx2']))
    if resolution and 'nx3' in resolution and resolution.get('nx3', 1) > 1:
        res_str_parts.append(str(resolution['nx3']))
    res_str = "x".join(res_str_parts) if res_str_parts else "default_res"
    
    print(f"--- Running Test: {test_id} ({res_str}) ---")

    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    output_dir_name = f"kh_{test_id}_{res_str}"
    absolute_output_dir_path = os.path.join(project_root, OUTPUT_BASE_DIR, output_dir_name)

    if os.path.exists(absolute_output_dir_path):
        print(f"Cleaning existing output directory: {absolute_output_dir_path}")
        shutil.rmtree(absolute_output_dir_path)
    os.makedirs(absolute_output_dir_path, exist_ok=True)
    print(f"Output will be in: {absolute_output_dir_path}")

    command = [
        ATHENA_EXECUTABLE,
        "-i", BASE_INPUT_FILE,
        "-d", os.path.join(OUTPUT_BASE_DIR, output_dir_name)
    ]

    if drat_val is not None: command.append(f"problem/drat={drat_val}")
    if bx0_val is not None: command.append(f"problem/bx0={bx0_val}")
    if by0_val is not None: command.append(f"problem/by0={by0_val}")
    if bz0_val is not None: command.append(f"problem/bz0={bz0_val}")

    if resolution:
        nx1 = resolution.get('nx1')
        nx2 = resolution.get('nx2')
        nx3 = resolution.get('nx3', 1) # Default to 1 if not specified for 2D cases
        if nx1: command.append(f"mesh/nx1={nx1}")
        if nx2: command.append(f"mesh/nx2={nx2}")
        if nx3: command.append(f"mesh/nx3={nx3}")
        
        # Set MeshBlock size to be the same as total mesh for simplicity in these tests
        # This ensures divisibility. For larger runs, one might use smaller blocks.
        if nx1: command.append(f"meshblock/nx1={nx1}")
        if nx2: command.append(f"meshblock/nx2={nx2}")
        if nx3: command.append(f"meshblock/nx3={nx3}")

    # --- Store parameters for analysis script ---
    # These are the parameters actually used for the simulation run
    run_params_for_analysis = {
        'test_id': test_id,
        'problem_generator': 'kh_bfield', # Hardcoded for now, could be from input file config
        'drat': drat_val if drat_val is not None else 1.0, # Default if not in test_params
        'bx0': bx0_val if bx0_val is not None else 0.0,
        'by0': by0_val if by0_val is not None else 0.0,
        'bz0': bz0_val if bz0_val is not None else 0.0,
        'resolution': resolution if resolution else {},
        'vflow': 0.5, # From kh_local.in, passed for analytic formulae
        'amp': 0.01,   # From kh_local.in, perturbation amplitude
        'domain_Lx': 1.0, # From kh_local.in <mesh> x1max-x1min
        'domain_Ly': 1.0, # From kh_local.in <mesh> x2max-x2min
        # Add other relevant params from BASE_INPUT_FILE if needed for theory checks
    }
    # Ensure defaults from kh_local.in are reflected if not overridden by test_params
    # Example: if test_params doesn't specify 'drat', kh_local.in has drat=1.0.
    # The GetOrAddReal in C++ defaults drat to 1.0. bx0 to 0.0.
    # The problem block in kh_local.in also defines these, which are used if not overridden.
    # The run_params_for_analysis tries to capture the *effective* parameters.

    # If specific resolution values were not in test_params, but are in kh_local.in, reflect those.
    # This is slightly redundant if resolution dict is always complete, but good for robustness.
    if not run_params_for_analysis['resolution'].get('nx1'):
        # These would be the defaults from kh_local.in if not overridden by CLI
        run_params_for_analysis['resolution']['nx1'] = 256 
        run_params_for_analysis['resolution']['nx2'] = 128
        run_params_for_analysis['resolution']['nx3'] = 1

    params_filepath = os.path.join(absolute_output_dir_path, "params.json")
    try:
        with open(params_filepath, 'w') as f_params:
            json.dump(run_params_for_analysis, f_params, indent=4)
        print(f"Saved simulation parameters to {params_filepath}")
    except IOError as e:
        print(f"Error saving simulation parameters to {params_filepath}: {e}")
        # Decide if this should be a fatal error for the test run

    print(f"Executing command: {' '.join(command)} from CWD: {project_root}")

    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False, cwd=project_root) # check=False initially
        
        stdout_content = process.stdout if process.stdout else ""
        stderr_content = process.stderr if process.stderr else ""

        if process.returncode != 0 or "FATAL ERROR" in stdout_content or "FATAL ERROR" in stderr_content:
            print(f"Error running test {test_id} (Return Code: {process.returncode}):")
            print("STDOUT:")
            print(stdout_content)
            print("STDERR:")
            print(stderr_content)
            return False
        else:
            print("STDOUT tail:")
            stdout_lines = stdout_content.splitlines()
            for line in stdout_lines[-20:]:
                print(line)
            print(f"Test {test_id} completed successfully.")

            # --- Call plotting script for the last HDF5 file ---
            # Find the last HDF5 file in the output directory
            # Athena HDF5 files are typically named like problem_id.outN.xxxxx.athdf
            # We need to list files, filter for .athdf, and sort to find the last one.
            try:
                output_files = os.listdir(absolute_output_dir_path)
                hdf5_files = sorted([f for f in output_files if f.endswith('.athdf')])
                hst_files = sorted([f for f in output_files if f.endswith('.hst')]) # Find .hst files
                
                # Plot last HDF5 snapshot
                if hdf5_files:
                    last_hdf5_file = os.path.join(absolute_output_dir_path, hdf5_files[-1])
                    snapshot_plotting_script_path = os.path.join(project_root, "vis/python/plot_kh_snapshot.py")
                    
                    print(f"Attempting to plot snapshot {last_hdf5_file} using {snapshot_plotting_script_path}")
                    snap_plot_command = ["python", snapshot_plotting_script_path, last_hdf5_file]
                    snap_plot_process = subprocess.run(snap_plot_command, capture_output=True, text=True, check=False, cwd=project_root)
                    
                    snap_plot_stdout = snap_plot_process.stdout if snap_plot_process.stdout else ""
                    snap_plot_stderr = snap_plot_process.stderr if snap_plot_process.stderr else ""

                    if snap_plot_process.returncode != 0 or "Error:" in snap_plot_stdout or "Error:" in snap_plot_stderr:
                        print(f"Snapshot plotting for {test_id} reported an error (Code: {snap_plot_process.returncode}):")
                        print("Snapshot Plot STDOUT:"); print(snap_plot_stdout)
                        print("Snapshot Plot STDERR:"); print(snap_plot_stderr)
                    else:
                        print(f"Snapshot plotting for {test_id} ({last_hdf5_file}) successful.")
                        print("Snapshot Plot STDOUT:"); print(snap_plot_stdout)
                else:
                    print(f"No HDF5 files found in {absolute_output_dir_path} to plot snapshot for test {test_id}.")

                # Plot HST metrics (e.g., energy drift)
                if hst_files:
                    hst_file_to_plot = os.path.join(absolute_output_dir_path, hst_files[0]) # Typically only one .hst file
                    hst_plotting_script_path = os.path.join(project_root, "vis/python/plot_hst_metrics.py")
                    
                    print(f"Attempting to plot HST metrics for {hst_file_to_plot} using {hst_plotting_script_path}")
                    hst_plot_command = ["python", hst_plotting_script_path, hst_file_to_plot]
                    hst_plot_process = subprocess.run(hst_plot_command, capture_output=True, text=True, check=False, cwd=project_root)
                    
                    hst_plot_stdout = hst_plot_process.stdout if hst_plot_process.stdout else ""
                    hst_plot_stderr = hst_plot_process.stderr if hst_plot_process.stderr else ""

                    if hst_plot_process.returncode != 0 or "Error:" in hst_plot_stdout or "Error:" in hst_plot_stderr:
                        print(f"HST metrics plotting for {test_id} reported an error (Code: {hst_plot_process.returncode}):")
                        print("HST Plot STDOUT:"); print(hst_plot_stdout)
                        print("HST Plot STDERR:"); print(hst_plot_stderr)
                    else:
                        print(f"HST metrics plotting for {test_id} ({hst_file_to_plot}) successful.")
                        print("HST Plot STDOUT:"); print(hst_plot_stdout)
                        print("HST Plot STDERR:"); print(hst_plot_stderr)
                else:
                    print(f"No HST files found in {absolute_output_dir_path} for HST metrics plotting for test {test_id}.")

            except Exception as e_plot:
                print(f"An error occurred during post-simulation plotting for {test_id}: {e_plot}")
            
            # --- Call analysis script ---
            try:
                analysis_script_path = os.path.join(project_root, "analysis/khi_analysis.py")
                print(f"Attempting to run analysis script {analysis_script_path} for test {test_id}")
                
                analysis_command = [
                    "python", 
                    analysis_script_path,
                    "--test_id", test_id,
                    "--output_dir", absolute_output_dir_path,
                    # Potentially add --problem_id if it varies and isn't always kh_bfield
                    # For now, khi_analysis.py defaults problem_id to "kh_bfield"
                ]
                
                # KHI analysis output should go to a dedicated file
                analysis_output_file = os.path.join(absolute_output_dir_path, f"{test_id}_analysis_results.txt")
                print(f"Full analysis results will be saved to: {analysis_output_file}")
                
                analysis_process = subprocess.run(analysis_command, capture_output=True, text=True, check=False, cwd=project_root)
                
                analysis_stdout = analysis_process.stdout if analysis_process.stdout else ""
                analysis_stderr = analysis_process.stderr if analysis_process.stderr else ""

                # Always save stdout to the analysis results file
                try:
                    with open(analysis_output_file, 'w') as f:
                        f.write(analysis_stdout)
                        if analysis_stderr:
                            f.write("\n\n--- STDERR ---\n\n")
                            f.write(analysis_stderr)
                except Exception as e:
                    print(f"Error saving analysis results to file: {e}")

                if analysis_process.returncode != 0 or "Error:" in analysis_stdout or "Error:" in analysis_stderr or "Traceback" in analysis_stderr:
                    print(f"KHI analysis for {test_id} reported an error (Code: {analysis_process.returncode}):")
                    print("Analysis STDERR:"); print(analysis_stderr)
                    # Don't print full stdout on error, just a summary
                    analysis_summary = "\n".join([
                        line for line in analysis_stdout.splitlines() 
                        if any(marker in line for marker in ["Result:", "M1", "M2", "M3", "M4", "T1", "T2", "T3", "T4", "T5"])
                    ])
                    print("Analysis Summary:"); print(analysis_summary)
                else:
                    print(f"KHI analysis for {test_id} completed successfully.")
                    # Print just the benchmark results as a summary
                    benchmark_results = "\n".join([
                        line for line in analysis_stdout.splitlines() 
                        if "Result:" in line
                    ])
                    print("Analysis Benchmark Results:"); print(benchmark_results)

            except Exception as e_analyze:
                print(f"An error occurred during KHI analysis script execution for {test_id}: {e_analyze}")

            return True # Test, snapshot plot, HST plot, and analysis all attempted
    except FileNotFoundError:
        print(f"Error: Athena executable not found at {os.path.join(project_root, ATHENA_EXECUTABLE)}")
        print("Make sure you have compiled Athena++.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred while running test {test_id}: {e}")
        return False

if __name__ == "__main__":
    all_tests_passed = True
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root_check = os.path.dirname(os.path.dirname(script_dir))
    
    print(f"Script run from: {os.getcwd()}")
    print(f"Calculated project root for CWD: {project_root_check}")
    print(f"Note: This script should be run from the project root directory.")

    base_output_abs_path = os.path.join(project_root_check, OUTPUT_BASE_DIR)
    if not os.path.exists(base_output_abs_path):
        print(f"Creating base output directory: {base_output_abs_path}")
        os.makedirs(base_output_abs_path)

    for test_case in test_cases:
        if not run_test(test_case):
            all_tests_passed = False
            print(f"!!! Test Failed: {test_case['id']} !!!")
    
    if all_tests_passed:
        print("\nAll tests completed successfully.")
    else:
        print("\nSome tests FAILED.") 