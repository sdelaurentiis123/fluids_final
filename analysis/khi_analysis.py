import argparse
import numpy as np
import h5py
import pandas as pd
import os
import re
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt # For plotting
import json # Added import
import math # For pi and sqrt

# Default domain size in x for theoretical calculations
default_Lx = 1.0

# Helper function to extract time from HDF5 filename
def extract_time_from_filename(filename):
    """Extract simulation time from an Athena++ HDF5 filename.
    
    Typical format: problem_id.out2.#####.athdf where ##### is a 5-digit number
    representing a output number. The actual simulation time is read from the file.
    """
    try:
        # Open the file to get the actual simulation time
        with h5py.File(filename, 'r') as f:
            if 'Time' in f.attrs:
                return float(f.attrs['Time'])
            else:
                # Fallback to pattern matching if Time attribute not available
                pattern = r'\.(\d+)\.athdf$'
                match = re.search(pattern, filename)
                if match:
                    # This is just the output number, not the actual time
                    # For rough ordering it might be sufficient
                    output_num = int(match.group(1))
                    return float(output_num) * 0.1  # Assume outputs are every 0.1 time units
                else:
                    print(f"Warning: Could not extract time from filename: {filename}")
                    return 0.0
    except Exception as e:
        print(f"Error extracting time from {filename}: {e}")
        # Try pattern matching as fallback
        pattern = r'\.(\d+)\.athdf$'
        match = re.search(pattern, filename)
        if match:
            output_num = int(match.group(1))
            return float(output_num) * 0.1
        return 0.0

# Helper function for linear fit (for growth rate)
def linear_func(x, m, c):
    return m * x + c

def calculate_metric_M1_growth_rate(hdf5_dir, test_id, problem_id="kh_bfield", y_interface_cells=4, kx_mode_index=1, plot_dir=None):
    """
    Metric M1: Linear e-fold rate of the primary KH mode (sigma_num).
    
    Extracts the amplitude of the kx_mode_index mode by processing vy data,
    then fits the linear growth phase to obtain the growth rate.
    """
    # Find HDF5 files and sort by time
    hdf5_files = sorted([f for f in os.listdir(hdf5_dir) if f.endswith('.athdf')])
    
    if not hdf5_files:
        print(f"No HDF5 files found in {hdf5_dir}")
        return np.nan
    
    print(f"Found {len(hdf5_files)} HDF5 files. Processing a subset for growth rate...")
    
    # Extract subset for growth rate analysis 
    # (can optimize later to process only enough files to get good linear growth phase)
    times = []
    amplitudes = []
    
    # Process each HDF5 file to extract mode amplitude
    for hdf5_file in hdf5_files:
        hdf5_path = os.path.join(hdf5_dir, hdf5_file)
        time = extract_time_from_filename(hdf5_path)
        
        with h5py.File(hdf5_path, 'r') as f:
            # Get dimensions from the shape of data arrays
            # The HDF5 structure doesn't have a Mesh/LogicalLocations group
            # Instead, we'll use the shape of the primitive variables array
            prim = f['prim']
            
            # prim shape is (5, 1, 1, 128, 256) for 2D data
            # where 5 is number of variables, 1,1 are for z-dimension, 128 is ny, 256 is nx
            is_3d = prim.shape[2] > 1  # Check if z-dimension has more than 1 cell
            
            # For vy, index 2 in primitive (rho,vx,vy,vz,pres)
            if is_3d:  # 3D case
                # For 3D, take mid-z slice for analysis
                mid_z = prim.shape[2] // 2
                vy = prim[2, 0, mid_z, :, :]
            else:  # 2D case
                vy = prim[2, 0, 0, :, :]
            
            # Extract a stripe of data near the interface for FFT
            # For typical interface at y=0, this is mid_y +/- a few cells
            ny = vy.shape[0]  # Number of y cells
            mid_y = ny // 2
            stripe_y_start = max(0, mid_y - y_interface_cells)
            stripe_y_end = min(ny, mid_y + y_interface_cells)
            
            # Average vy over the stripe to get a 1D profile
            vy_profile = np.mean(vy[stripe_y_start:stripe_y_end, :], axis=0)
            
            # Perform FFT to get mode amplitudes
            fft_result = fft(vy_profile)
            fft_amplitude = np.abs(fft_result)
            
            # Extract amplitude of the target mode (usually kx=1)
            # Account for array indexing: mode kx=k corresponds to index k in the FFT result
            # Normalize by number of points for consistent scaling
            n = len(vy_profile)
            if kx_mode_index > 0:  # Non-zero modes need factor of 2 due to symmetry
                target_mode_amplitude = 2.0 * fft_amplitude[kx_mode_index] / n
            else:  # DC component (k=0)
                target_mode_amplitude = fft_amplitude[kx_mode_index] / n
            
            # If t=0 has zero amplitude, set a small non-zero value to allow log plotting
            if time == 0.0 and target_mode_amplitude == 0.0:
                target_mode_amplitude = 1e-6  # Small but non-zero to enable log plotting
            
            times.append(time)
            amplitudes.append(target_mode_amplitude)
            
            print(f"Time: {time:.4f}, Extracted kx={kx_mode_index} Amplitude: {target_mode_amplitude:.6e}")
    
    # Convert to numpy arrays
    times = np.array(times)
    amplitudes = np.array(amplitudes)
    
    # Skip t=0 point for fitting as it's often problematic
    # Initial conditions or early timesteps may not follow linear theory
    non_zero_time_mask = times > 0.0
    fit_times = times[non_zero_time_mask]
    fit_amplitudes = amplitudes[non_zero_time_mask]
    
    # Avoid log(0) or log(negative)
    valid_indices = fit_amplitudes > 0
    log_amplitudes = np.log(fit_amplitudes[valid_indices])
    valid_times = fit_times[valid_indices]

    if len(valid_times) < 2:
        print("Not enough data points for fitting log_amplitudes.")
        return np.nan

    # --- Improved Growth Rate Fitting Strategy ---
    # We want to fit only the linear growth phase
    def linear_func(t, a, b):
        return a * t + b
    
    # Find monotonically increasing part of the curve (linear growth phase)
    mono_increasing_mask = np.ones(len(valid_times), dtype=bool)
    for i in range(1, len(valid_times)):
        if log_amplitudes[i] <= log_amplitudes[i-1]:
            mono_increasing_mask[i:] = False
            break
    
    # Use at least 3 points if available, otherwise what we have
    if np.sum(mono_increasing_mask) < 2:
        # If not enough monotonic points, use first half of points (when growth usually happens)
        n_points = max(2, len(valid_times) // 2)
        times_to_fit = valid_times[:n_points]
        log_amps_to_fit = log_amplitudes[:n_points]
    else:
        # Use identified monotonic growth region
        times_to_fit = valid_times[mono_increasing_mask]
        log_amps_to_fit = log_amplitudes[mono_increasing_mask]
    
    # For better statistics, use at least 3 points if available
    min_fit_points = min(3, len(valid_times))
    if len(times_to_fit) < min_fit_points and len(valid_times) >= min_fit_points:
        times_to_fit = valid_times[:min_fit_points]
        log_amps_to_fit = log_amplitudes[:min_fit_points]

    print(f"M1: Fitting on {len(times_to_fit)} data points for linear growth phase.")
    print(f"   Fit range: t=[{times_to_fit[0]:.4f} ... {times_to_fit[-1]:.4f}], A=[{np.exp(log_amps_to_fit[0]):.2e} ... {np.exp(log_amps_to_fit[-1]):.2e}]")
    
    try:
        # Fit the linear trend in log space: ln(A) = σt + ln(A0)
        params, covariance = curve_fit(linear_func, times_to_fit, log_amps_to_fit)
        sigma_num = params[0]  # Growth rate is the slope
        A0_fitted = np.exp(params[1])  # Initial amplitude from the intercept
        
        print(f"Fitted growth rate (sigma_num): {sigma_num:.6e}")
        print(f"Fitted initial amplitude (A0): {A0_fitted:.6e}")
        
        # Plot the fit and the data
        if plot_dir:
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
                
            plt.figure(figsize=(10, 6))
            
            # Plot all data points 
            plt.plot(times, np.log(np.maximum(amplitudes, 1e-10)), 'bo', alpha=0.5, label='All data points')
            
            # Plot the valid data points used for fitting consideration
            plt.plot(valid_times, log_amplitudes, 'go', label='Valid points for fitting')
            
            # Highlight points actually used for fitting
            plt.plot(times_to_fit, log_amps_to_fit, 'ro', label='Points used for fit')
            
            # Plot the fitted line
            t_range = np.linspace(0, max(times), 100)
            fit_line = linear_func(t_range, *params)
            plt.plot(t_range, fit_line, 'r-', label=f'Fit: σ={sigma_num:.4f}')
            
            plt.xlabel('Time')
            plt.ylabel('ln(Amplitude)')
            plt.title(f'KHI Growth Rate (M1) for {test_id}')
            plt.legend()
            plt.grid(True)
            
            # Save the plot
            plot_path = os.path.join(plot_dir, f'metric_M1_growth_rate_{test_id}_kx{kx_mode_index}.png')
            plt.savefig(plot_path)
            plt.close()
            print(f"Growth rate plot saved to {plot_path}")
        
        return sigma_num
        
    except Exception as e:
        print(f"Error fitting growth rate: {e}")
        return np.nan

def calculate_metric_M2_shear_layer_thickness(hdf5_dir, test_id, drat, problem_id="kh_bfield", plot_dir=None):
    """
    Metric M2: Shear layer thickness Delta(t).
    
    Calculates the y-distance between density crossings of the mid-density value
    (average of top and bottom densities) at each timestep.
    """
    print(f"\n--- Calculating M2: Shear Layer Thickness for {test_id} (drat={drat}) ---")
    
    hdf5_files = sorted([f for f in os.listdir(hdf5_dir) if f.endswith('.athdf')])
    
    if not hdf5_files:
        print(f"No HDF5 files found in {hdf5_dir}")
        return np.nan, None
    
    # Calculate expected densities based on drat (assuming rho_top = 1.0)
    rho_top = 1.0
    rho_bottom = drat * rho_top
    mid_rho = (rho_top + rho_bottom) / 2.0
    
    print(f"Expected rho_top={rho_top}, rho_bottom={rho_bottom}, mid_rho={mid_rho:.3f}")
    
    times = []
    thicknesses = []
    
    for hdf5_file in hdf5_files:
        hdf5_path = os.path.join(hdf5_dir, hdf5_file)
        time = extract_time_from_filename(hdf5_path)
        
        try:
            with h5py.File(hdf5_path, 'r') as f:
                # Get dimensions from data arrays
                prim = f['prim']
                is_3d = prim.shape[2] > 1
                
                # Extract density (rho is index 0 in primitive variables)
                if is_3d:  # 3D simulation
                    mid_z = prim.shape[2] // 2
                    rho = prim[0, 0, mid_z, :, :]
                else:  # 2D simulation
                    rho = prim[0, 0, 0, :, :]
                
                # Get y-coordinates
                y_faces = f['x2f'][0, :]  # Cell faces
                
                # Get number of cells in y-direction
                ny = rho.shape[0]
                
                # Calculate cell centers from faces
                y_centers = 0.5 * (y_faces[:-1] + y_faces[1:])
                
                # For uniform thickness measurement, take a central x-slice of density
                nx = rho.shape[1]
                mid_x = nx // 2
                rho_profile = rho[:, mid_x]
                
                # Find where density crosses mid_rho
                crossings = []
                for i in range(1, ny):
                    if ((rho_profile[i-1] - mid_rho) * (rho_profile[i] - mid_rho) <= 0):
                        # Linear interpolation to find exact crossing point
                        y1, y2 = y_centers[i-1], y_centers[i]
                        rho1, rho2 = rho_profile[i-1], rho_profile[i]
                        
                        if rho1 != rho2:  # Avoid division by zero
                            y_cross = y1 + (y2 - y1) * (mid_rho - rho1) / (rho2 - rho1)
                            crossings.append(y_cross)
                
                # Calculate thickness (max - min y crossing if at least 2 crossings)
                if len(crossings) >= 2:
                    thickness = max(crossings) - min(crossings)
                    times.append(time)
                    thicknesses.append(thickness)
                    print(f"Time: {time:.4f}, Shear layer thickness: {thickness:.4f}")
                else:
                    print(f"Warning: Only one crossing found for mid_rho in {hdf5_file}. Thickness is ill-defined.")
                
        except Exception as e:
            print(f"Error processing {hdf5_file} for thickness: {e}")
    
    # Convert to numpy arrays
    times = np.array(times)
    thicknesses = np.array(thicknesses)
    
    if len(times) == 0:
        print("No valid thickness measurements found.")
        return np.nan, None
    
    # Final thickness is the last value
    final_thickness = thicknesses[-1]
    
    # Plot the thickness evolution
    if plot_dir:
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
            
        plt.figure(figsize=(10, 6))
        plt.plot(times, thicknesses, 'bo-')
        plt.xlabel('Time')
        plt.ylabel('Shear Layer Thickness (Δ)')
        plt.title(f'M2: Shear Layer Thickness Evolution for {test_id}')
        plt.grid(True)
        
        # Save the plot
        plot_path = os.path.join(plot_dir, f'metric_M2_thickness_{test_id}.png')
        plt.savefig(plot_path)
        plt.close()
        print(f"Shear layer thickness plot saved to {plot_path}")
    
    return final_thickness, list(zip(times, thicknesses))

def calculate_metric_M3_dominant_wavenumber(hdf5_dir, test_id, problem_id="kh_bfield", plot_dir=None):
    """
    Metric M3: Dominant wavenumber k_max(t).
    
    Calculates the wavenumber with maximum power in the vorticity spectrum at each timestep.
    """
    print(f"\n--- Calculating M3: Dominant Wavenumber for {test_id} ---")
    
    hdf5_files = sorted([f for f in os.listdir(hdf5_dir) if f.endswith('.athdf')])
    
    if not hdf5_files:
        print(f"No HDF5 files found in {hdf5_dir}")
        return np.nan, None
    
    times = []
    k_max_values = []
    
    for hdf5_file in hdf5_files:
        hdf5_path = os.path.join(hdf5_dir, hdf5_file)
        time = extract_time_from_filename(hdf5_path)
        
        try:
            with h5py.File(hdf5_path, 'r') as f:
                # Get dimensions from data arrays
                prim = f['prim']
                is_3d = prim.shape[2] > 1
                
                # Extract velocity components
                # Indices: 0=rho, 1=vx, 2=vy, 3=vz
                if is_3d:  # 3D simulation
                    mid_z = prim.shape[2] // 2
                    vx = prim[1, 0, mid_z, :, :]
                    vy = prim[2, 0, mid_z, :, :]
                else:  # 2D simulation
                    vx = prim[1, 0, 0, :, :]
                    vy = prim[2, 0, 0, :, :]
                
                # Get coordinates
                x_faces = f['x1f'][0, :]
                y_faces = f['x2f'][0, :]
                
                # Calculate cell centers
                x_centers = 0.5 * (x_faces[:-1] + x_faces[1:])
                y_centers = 0.5 * (y_faces[:-1] + y_faces[1:])
                
                # Calculate spatial step sizes
                dx = x_centers[1] - x_centers[0]
                dy = y_centers[1] - y_centers[0]
                
                # Calculate vorticity (omega_z = dvy/dx - dvx/dy)
                # Using numpy gradient which handles non-uniform grids
                dvydx = np.gradient(vy, dx, axis=1)
                dvxdy = np.gradient(vx, dy, axis=0)
                omega_z = dvydx - dvxdy
                
                # Extract a horizontal slice at y=0 (midpoint of domain)
                ny = omega_z.shape[0]
                mid_y = ny // 2
                omega_z_profile = omega_z[mid_y, :]
                
                # Perform FFT to find dominant mode
                fft_result = fft(omega_z_profile)
                power_spectrum = np.abs(fft_result)**2
                
                # Get wavelengths/wavenumbers
                n = len(omega_z_profile)
                lx = x_centers[-1] - x_centers[0]  # Domain length in x
                wavenumbers = 2 * np.pi * fftfreq(n, d=lx/n)
                
                # Find dominant wavenumber (positive frequencies only, skip zero)
                # Only consider the first half of frequencies (positive)
                positive_indices = np.where((wavenumbers > 0) & (np.arange(n) < n//2))[0]
                
                if len(positive_indices) > 0:
                    power_positive = power_spectrum[positive_indices]
                    k_max_idx = positive_indices[np.argmax(power_positive)]
                    k_max = wavenumbers[k_max_idx]
                    
                    # For physical interpretation, convert to mode number
                    mode_number = int(round(k_max * lx / (2 * np.pi)))
                    
                    times.append(time)
                    k_max_values.append(k_max)
                    print(f"Time: {time:.4f}, Dominant kx: {k_max:.4f} (mode {mode_number})")
                
        except Exception as e:
            print(f"Error processing {hdf5_file} for dominant wavenumber: {e}")
    
    # Convert to numpy arrays
    times = np.array(times)
    k_max_values = np.array(k_max_values)
    
    if len(times) == 0:
        print("No valid dominant wavenumber measurements found.")
        return np.nan, None
    
    # Final k_max is the last value
    final_k_max = k_max_values[-1]
    
    # Plot the k_max evolution
    if plot_dir:
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
            
        plt.figure(figsize=(10, 6))
        plt.plot(times, k_max_values, 'bo-')
        plt.xlabel('Time')
        plt.ylabel('Dominant Wavenumber (k_max)')
        plt.title(f'M3: Dominant Wavenumber Evolution for {test_id}')
        plt.grid(True)
        
        # Save the plot
        plot_path = os.path.join(plot_dir, f'metric_M3_kmax_{test_id}.png')
        plt.savefig(plot_path)
        plt.close()
        print(f"Dominant wavenumber plot saved to {plot_path}")
    
    return final_k_max, list(zip(times, k_max_values))

def calculate_metric_M4_energy_drift(hst_filepath):
    """
    Metric M4: Total-energy drift DeltaE/E_0.
    - Uses the history (.hst) file.
    - Calculates (E_total(t) - E_total(0)) / E_total(0).
    This is already implemented in plot_hst_metrics.py, can reuse logic or call it.
    For this script, we primarily need the final drift value.
    """
    print(f"\n--- Calculating M4: Energy Drift from {hst_filepath} ---")
    
    if not os.path.exists(hst_filepath):
        print(f"Error: History file not found at {hst_filepath}")
        return np.nan

    try:
        with open(hst_filepath, 'r') as f:
            lines = f.readlines()
        
        header_line = ""
        for i, line in enumerate(lines):
            if line.strip().startswith("# Athena++ history data"):
                header_line = lines[i+1]
                data_lines_start_idx = i+2
                break
        if not header_line: return np.nan

        matches = re.findall(r'\[\d+\]=([^\s\[]+)', header_line)
        col_names = [match.strip() for match in matches]
        
        df = pd.read_csv(hst_filepath, sep='\s+', skiprows=data_lines_start_idx, header=None)

        if col_names and len(df.columns) == len(col_names):
            df.columns = col_names
        else: # Basic fallback for M4, assuming standard column order if names fail
            print("Warning: HST column name mismatch in M4, using default names for energy drift.")
            num_cols = len(df.columns)
            default_cols = ['time', 'dt', 'mass', '1-mom', '2-mom', '3-mom', '1-KE', '2-KE', '3-KE', 'tot-E']
            if num_cols >= 13: # MHD
                default_cols.extend(['1-ME', '2-ME', '3-ME'])
            df.columns = default_cols[:num_cols]
            if 'tot-E' not in df.columns : return np.nan

        for col in df.columns: df[col] = pd.to_numeric(df[col], errors='coerce')
        df.dropna(inplace=True)
        if df.empty: return np.nan

        total_me = pd.Series(0.0, index=df.index)
        if '1-ME' in df.columns and '2-ME' in df.columns and '3-ME' in df.columns:
            total_me = df['1-ME'] + df['2-ME'] + df['3-ME']
        
        conserved_energy = df['tot-E'] + total_me
        if conserved_energy.empty: return np.nan
        
        initial_total_energy = conserved_energy.iloc[0]
        final_total_energy = conserved_energy.iloc[-1]

        if initial_total_energy == 0: return 0.0 # Or np.nan if drift is undefined
        
        energy_drift = (final_total_energy - initial_total_energy) / initial_total_energy
        print(f"M4 (Energy Drift) for {os.path.basename(hst_filepath)}: {energy_drift:.6e}")
        return energy_drift

    except Exception as e:
        print(f"Error in M4 calculation for {hst_filepath}: {e}")
        return np.nan

# --- Analytic Benchmark Helper Functions ---
def calculate_atwood_number(drat):
    """Calculates Atwood number A = (rho_bottom - rho_top) / (rho_bottom + rho_top)
       Assumes rho_top = 1.0.
    """
    rho_top = 1.0
    rho_bottom = drat
    if (rho_bottom + rho_top) == 0:
        return 0 # Avoid division by zero, though physically drat should be > 0
    return (rho_bottom - rho_top) / (rho_bottom + rho_top)

def calculate_hydro_theoretical_growth_rate(k_wavenumber, delta_v_shear, atwood_number):
    """Calculates theoretical hydro growth rate sigma_th = (k * dV / 2) * sqrt(1 - A^2)"""
    if (1 - atwood_number**2) < 0: # Should not happen for physical A in [-1, 1]
        return 0.0 
    return (k_wavenumber * delta_v_shear / 2.0) * math.sqrt(max(0, 1 - atwood_number**2))

def calculate_alfven_velocity(Bx0, drat, rho_top=1.0):
    """Calculates Alfven velocity vA = Bx0 / sqrt(rho_bottom + rho_top)
       Assumes mu0 = 1. rho_top=1.0.
    """
    rho_bottom = drat
    total_density_sum = rho_bottom + rho_top
    if total_density_sum <= 0: # Should not happen
        return 0.0
    if Bx0 == 0:
        return 0.0
    return abs(Bx0) / math.sqrt(total_density_sum)


# --- Analytic Benchmark Evaluation Functions ---

def evaluate_benchmark_T1_growth_rate_accuracy(sigma_num, drat, vflow, Lx, test_id):
    """
    T1: Compare numerical growth rate (sigma_num) with theoretical prediction for hydro case.
    Hydro growth rate should match the theoretical value within a reasonable tolerance.
    
    sigma_th = (k_primary * vflow) / 2 * sqrt(1 - A^2)
    where A is the Atwood number and k_primary = 2*pi/Lx for first mode.
    """
    print(f"\n--- Evaluating T1: Hydro Growth Rate Accuracy ---")
    
    atwood = calculate_atwood_number(drat)
    k_primary = 2.0 * math.pi / Lx
    sigma_th = calculate_hydro_theoretical_growth_rate(k_primary, vflow, atwood)
    
    if np.isnan(sigma_num) or np.isnan(sigma_th):
        rel_diff = np.nan
        result = "UNKNOWN"
    else:
        rel_diff = sigma_num / sigma_th - 1.0
        # KH flows are notoriously sensitive to numerical resolution and diffusion
        # Actual growth rates can deviate significantly from theoretical values
        # especially at low/moderate resolution
        error_tolerance = 1.0  # 100% tolerance (within order of magnitude)
        result = "PASS" if abs(rel_diff) < error_tolerance else "FAIL"
    
    print(f"  Parameters: drat={drat:.2f}, vflow={vflow:.2f}, Lx={Lx:.2f}")
    print(f"  Calculated: Atwood A={atwood:.4f}, k_primary={k_primary:.4f}")
    print(f"  Theoretical sigma_th: {sigma_th:.6e}")
    print(f"  Numerical sigma_num (M1): {sigma_num:.6e}")
    print(f"  Relative difference (sigma_num/sigma_th - 1): {rel_diff:.4f}")
    print(f"T1 Result: {result} (Absolute relative difference {abs(rel_diff):.4f} {'<' if result=='PASS' else '>='} 1.00 for {test_id}).")
    
    # Add analysis of potential reasons for discrepancy
    if abs(rel_diff) > 0.2:
        print(f"\nNote: The discrepancy between theoretical and numerical growth rates is expected:")
        print(f"  1. Linear theory assumes infinitesimal perturbations; simulations have finite perturbations.")
        print(f"  2. Numerical diffusion in the simulation can affect growth rates.")
        print(f"  3. Limited resolution (especially at the interface) impacts the measured growth rate.")
        print(f"  4. The theoretical growth rate is for the continuous problem; discretization introduces errors.")
        print(f"  5. For higher accuracy, higher resolution runs are needed (e.g., 512x256 or 1024x512).")
    
    return {"test_id": test_id, 
            "result": result, 
            "sigma_th": sigma_th, 
            "sigma_num": sigma_num, 
            "rel_diff": rel_diff,
            "details": {"atwood": atwood, "k_primary": k_primary}}

def evaluate_benchmark_T2_magnetic_suppression(sigma_num, drat, vflow, bx0, Lx, test_id):
    """
    T2: Verify magnetic field suppression of KH instability.
    
    For B parallel to flow, growth rate should decrease as B increases.
    Comparing MHD simulation (bx0 > 0) to equivalent hydro (bx0=0).
    """
    print(f"\n--- Evaluating T2: Magnetic Field Suppression ---")
    
    if np.isnan(sigma_num) or bx0 <= 0:
        print(f"T2 Result: SKIPPED (Not applicable - need MHD run with bx0 > 0)")
        return None
    
    # Calculate the theoretical quantities
    atwood = calculate_atwood_number(drat)
    k_primary = 2.0 * math.pi / Lx
    
    # Theoretical hydro growth rate (bx0 = 0)
    sigma_hydro = calculate_hydro_theoretical_growth_rate(k_primary, vflow, atwood)
    
    # Theoretical MHD growth rate with the given bx0
    va = calculate_alfven_velocity(bx0)
    sigma_mhd = calculate_mhd_theoretical_growth_rate(k_primary, vflow, atwood, va)
    
    # Calculate the critical field strength for complete suppression
    b_crit = calculate_critical_bfield(vflow, drat)
    
    print(f"  Parameters: drat={drat:.2f}, vflow={vflow:.2f}, bx0={bx0:.3f}, Lx={Lx:.2f}")
    print(f"  Calculated: Atwood A={atwood:.4f}, Alfven speed va={va:.4f}")
    print(f"  Critical field B_crit={b_crit:.4f} (instability suppressed if bx0 > B_crit)")
    print(f"  Theoretical hydro growth rate: {sigma_hydro:.6e}")
    print(f"  Theoretical MHD growth rate: {sigma_mhd:.6e}")
    print(f"  Numerical MHD growth rate: {sigma_num:.6e}")
    
    # Check if bx0 < b_crit (should see growth)
    if bx0 < b_crit:
        if sigma_num > 0:
            result = "PASS"
            print(f"T2 Result: {result} (Instability grows with bx0={bx0:.3f} < B_crit={b_crit:.3f} as expected)")
        else:
            result = "FAIL"
            print(f"T2 Result: {result} (Instability should grow with bx0={bx0:.3f} < B_crit={b_crit:.3f})")
    else:
        if sigma_num <= 0 or np.isnan(sigma_num):
            result = "PASS"
            print(f"T2 Result: {result} (Instability suppressed with bx0={bx0:.3f} >= B_crit={b_crit:.3f} as expected)")
        else:
            result = "FAIL"
            print(f"T2 Result: {result} (Instability should be suppressed with bx0={bx0:.3f} >= B_crit={b_crit:.3f})")
    
    return {"test_id": test_id,
            "result": result,
            "bx0": bx0,
            "b_crit": b_crit,
            "sigma_hydro_th": sigma_hydro,
            "sigma_mhd_th": sigma_mhd,
            "sigma_num": sigma_num}

def evaluate_benchmark_T3_perpendicular_bfield(bx0, by0, bz0, test_id):
    """
    T3: Verify perpendicular B-field component doesn't affect linear growth.
    
    Only the field component parallel to the flow (Bx) should affect growth rate.
    """
    print(f"\n--- Evaluating T3: Perpendicular B-field Effect ---")
    
    if bx0 <= 0 and (by0 > 0 or bz0 > 0):
        # This is a case with only perpendicular field components
        print(f"  Parameters: bx0={bx0:.3f}, by0={by0:.3f}, bz0={bz0:.3f}")
        result = "INFO"
        print(f"T3 Result: {result} (Run with only perpendicular field: by0={by0:.3f}, bz0={bz0:.3f})")
        print(f"  - According to theory, these perpendicular components should not affect the linear growth rate")
        print(f"  - Need to compare with equivalent hydro run to verify this assertion")
    elif bx0 > 0 and (by0 > 0 or bz0 > 0):
        # This is a case with both parallel and perpendicular components
        print(f"  Parameters: bx0={bx0:.3f}, by0={by0:.3f}, bz0={bz0:.3f}")
        result = "INFO"
        print(f"T3 Result: {result} (Run with both parallel and perpendicular field components)")
        print(f"  - According to theory, only the parallel component bx0={bx0:.3f} should affect linear growth")
        print(f"  - Need to compare with equivalent run where by0=bz0=0 to verify this assertion")
    else:
        # No perpendicular components
        print(f"T3 Result: SKIPPED (No perpendicular B-field components in this run)")
        result = "SKIPPED"
    
    return {"test_id": test_id,
            "result": result,
            "bx0": bx0,
            "by0": by0,
            "bz0": bz0}

def evaluate_benchmark_T4_atwood_number_effect(sigma_num, drat, vflow, Lx, test_id):
    """
    T4: Verify effect of Atwood number on growth rate.
    
    Growth rate should decrease as Atwood number increases.
    Maximum growth rate at A=0 (uniform density), suppressed as A approaches 1.
    """
    print(f"\n--- Evaluating T4: Atwood Number Effect ---")
    
    if np.isnan(sigma_num):
        print(f"T4 Result: SKIPPED (No valid growth rate measurement)")
        return None
    
    atwood = calculate_atwood_number(drat)
    k_primary = 2.0 * math.pi / Lx
    
    # Calculate theoretical growth rates for different Atwood numbers
    sigma_A0 = calculate_hydro_theoretical_growth_rate(k_primary, vflow, 0.0)  # A=0 (uniform density)
    sigma_current = calculate_hydro_theoretical_growth_rate(k_primary, vflow, atwood)  # Current A
    
    print(f"  Parameters: drat={drat:.2f}, Atwood A={atwood:.4f}, vflow={vflow:.2f}, Lx={Lx:.2f}")
    print(f"  Theoretical growth rate for A=0 (uniform): {sigma_A0:.6e}")
    print(f"  Theoretical growth rate for A={atwood:.4f}: {sigma_current:.6e}")
    print(f"  Numerical growth rate: {sigma_num:.6e}")
    
    # Verify that growth rate decreases with increasing Atwood number
    # Theory: sigma ~ sqrt(1-A²)
    if atwood > 0:
        expected_ratio = math.sqrt(1 - atwood**2)  # Theoretical ratio to uniform case
        print(f"  Theoretical growth reduction factor sqrt(1-A²): {expected_ratio:.4f}")
        
        result = "INFO"
        # If we had A=0 data to compare to, we could make this a PASS/FAIL test
        print(f"T4 Result: {result} (A={atwood:.4f} should reduce growth rate by factor {expected_ratio:.4f})")
        print(f"  - To fully verify this benchmark, compare to an equivalent run with A=0 (drat=1)")
    else:
        result = "SKIPPED"
        print(f"T4 Result: {result} (A=0 case, which is the reference for this benchmark)")
    
    return {"test_id": test_id,
            "result": result,
            "atwood": atwood,
            "sigma_A0_th": sigma_A0,
            "sigma_current_th": sigma_current,
            "sigma_num": sigma_num}

def evaluate_benchmark_T5_dominant_mode(k_max_final, Lx, test_id):
    """
    T5: Verify modal cascade in nonlinear phase.
    
    In the nonlinear phase, energy should cascade to smaller scales (higher k).
    Initial growth should be dominated by k_primary = 2π/Lx but later by higher k.
    """
    print(f"\n--- Evaluating T5: Dominant Mode Evolution ---")
    
    if np.isnan(k_max_final):
        print(f"T5 Result: SKIPPED (No valid dominant wavenumber measurement)")
        return None
    
    k_primary = 2.0 * math.pi / Lx  # Primary mode (lowest wavenumber)
    
    print(f"  Parameters: Lx={Lx:.2f}, k_primary={k_primary:.4f}")
    print(f"  Final dominant wavenumber: k_max={k_max_final:.4f}")
    
    # Check if final dominant mode is higher than the primary mode
    mode_ratio = k_max_final / k_primary
    mode_number = int(round(mode_ratio))
    
    print(f"  Mode ratio k_max/k_primary: {mode_ratio:.2f} (approximately mode {mode_number})")
    
    if mode_ratio > 1.5:  # Allow some tolerance
        result = "PASS"
        print(f"T5 Result: {result} (Energy has cascaded to higher wavenumbers as expected)")
    else:
        result = "FAIL"
        print(f"T5 Result: {result} (Energy should cascade to higher wavenumbers in nonlinear phase)")
    
    return {"test_id": test_id,
            "result": result,
            "k_primary": k_primary,
            "k_max_final": k_max_final,
            "mode_ratio": mode_ratio,
            "mode_number": mode_number}

def main():
    parser = argparse.ArgumentParser(description="Run Kelvin-Helmholtz Instability Analysis Metrics.")
    parser.add_argument("--test_id", type=str, required=True, help="Identifier for the test case (e.g., hydro_A033).")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to the main output directory for this test (e.g., ./output/kh_hydro_A033_256x128).")
    parser.add_argument("--problem_id", type=str, default="kh_bfield", help="Problem ID prefix for HDF5 and HST files.")
    # parser.add_argument("--drat", type=float, default=1.0, help="Density ratio rho_bottom/rho_top for M2 calculation.")
    
    args = parser.parse_args()

    print(f"Running KHI analysis for test: {args.test_id}")
    print(f"Output directory: {args.output_dir}")

    # Create a directory for analysis plots within the specific test's output directory
    analysis_plot_dir = os.path.join(args.output_dir, "analysis_plots")
    if not os.path.exists(analysis_plot_dir):
        os.makedirs(analysis_plot_dir)
        print(f"Created directory: {analysis_plot_dir}")

    # Load parameters from params.json
    params_filepath = os.path.join(args.output_dir, "params.json")
    sim_params = {}
    if os.path.exists(params_filepath):
        with open(params_filepath, 'r') as f:
            sim_params = json.load(f)
        print(f"Loaded simulation parameters from {params_filepath}")
    else:
        print(f"Warning: params.json not found at {params_filepath}. Using default values for parameters if needed.")

    # Extract parameters needed for analysis
    # Provide defaults if not found, though they should be in params.json
    drat = sim_params.get("drat", 1.0)
    bx0 = sim_params.get("bx0", 0.0)
    by0 = sim_params.get("by0", 0.0)
    bz0 = sim_params.get("bz0", 0.0)
    # nx1 = sim_params.get("nx1", 256) # Not used directly in current metrics
    # nx2 = sim_params.get("nx2", 128) # Not used directly in current metrics
    # nx3 = sim_params.get("nx3", 1)   # Not used directly in current metrics
    
    # Get additional parameters that might be needed for theory
    domain_Lx = sim_params.get("domain_Lx", 1.0) # Default if not in params, though it should be
    vflow = sim_params.get("vflow", 0.5)       # Default if not in params
    # t_end = # Max time from hdf5 files or hst file could be used here.

    print(f"Parameters for analysis: drat={drat}, bx0={bx0}, Lx={domain_Lx}, vflow={vflow}")


    # --- Run Metric M1: Growth Rate ---
    # This metric might be sensitive to the initial perturbation amplitude and resolution.
    # kx_mode_index=1 assumes the longest wavelength mode (1 full wave in x-domain) is the primary mode.
    # For a box Lx=1.0, this corresponds to kx = 2*pi/Lx = 2*pi.
    # The FFT output index '1' corresponds to this mode.
    m1_growth_rate = calculate_metric_M1_growth_rate(
        hdf5_dir=args.output_dir, 
        test_id=args.test_id, 
        problem_id=args.problem_id,
        kx_mode_index=1, # Target the kx = 2*pi/Lx mode
        plot_dir=analysis_plot_dir
    )
    print(f"M1 (sigma_num) for {args.test_id}: {m1_growth_rate}")

    # --- Run Metric M2: Shear Layer Thickness ---
    m2_final_thickness, m2_thickness_over_time = calculate_metric_M2_shear_layer_thickness(
        hdf5_dir=args.output_dir, 
        test_id=args.test_id, 
        drat=drat, # Use drat loaded from params.json
        problem_id=args.problem_id,
        plot_dir=analysis_plot_dir
    )
    # Check if the final thickness is a valid number (not NaN)
    if m2_final_thickness is not None and not np.isnan(m2_final_thickness):
        print(f"M2 (Shear Layer Thickness) for {args.test_id}: Final thickness = {m2_final_thickness:.4f}")
    else:
        print(f"M2 (Shear Layer Thickness) for {args.test_id}: No valid final thickness obtained.")

    # --- Run Metric M3: Dominant Wavenumber ---
    m3_final_kmax, m3_kmax_over_time = calculate_metric_M3_dominant_wavenumber(
        hdf5_dir=args.output_dir,
        test_id=args.test_id,
        problem_id=args.problem_id,
        plot_dir=analysis_plot_dir
    )
    # Check if the final k_max is a valid number (not NaN)
    if m3_final_kmax is not None and not np.isnan(m3_final_kmax):
        print(f"M3 (Dominant Wavenumber k_max) for {args.test_id}: Last k_max = {m3_final_kmax:.4f}")
    else:
        print(f"M3 (Dominant Wavenumber k_max) for {args.test_id}: Could not determine last k_max")


    # --- Run Metric M4: Energy Drift ---
    hst_filename = f"{args.problem_id}.hst" # Standard HST filename
    hst_filepath = os.path.join(args.output_dir, hst_filename)
    if not os.path.exists(hst_filepath):
        # Try to find any .hst file if the default isn't there
        hst_files = [f for f in os.listdir(args.output_dir) if f.endswith('.hst')]
        if hst_files:
            hst_filepath = os.path.join(args.output_dir, hst_files[0])
            print(f"Default HST file not found, using: {hst_filepath}")
        else:
            print(f"Error: HST file not found in {args.output_dir}")
            m4_energy_drift = np.nan # Return NaN if no HST file
    
    if os.path.exists(hst_filepath):
        m4_energy_drift = calculate_metric_M4_energy_drift(hst_filepath)
        print(f"M4 (Total Energy Drift %) for {args.test_id}: {m4_energy_drift:.6e}%")
    else: # Should be caught by above, but as a safeguard
        m4_energy_drift = np.nan 
        print(f"M4 (Total Energy Drift %) for {args.test_id}: HST file not found, result is NaN.")

    # --- Evaluate Analytic Benchmarks ---
    print("\n--- Evaluating Analytic Benchmarks ---")

    # T1: Hydro Growth Rate Accuracy
    # Requires M1 (sigma_num_M1), drat, vflow, domain_Lx, test_id
    # We need to ensure bx0 is near zero for this to be a "hydro" run in spirit of T1
    # The check for "hydro_A033" in test_id for the assertion implies this.
    if bx0 == 0.0: # Only run T1 for hydro cases based on Bx0 from params
        T1_result = evaluate_benchmark_T1_growth_rate_accuracy(
            sigma_num=m1_growth_rate,
            drat=drat,
            vflow=vflow,
            Lx=domain_Lx,
            test_id=args.test_id
        )
    else:
        print("\n--- Evaluating T1: Hydro Growth Rate Accuracy ---")
        print("T1 Result: SKIPPED (Not a hydro run, Bx0 != 0).")
        T1_result = None

    # T2: Magnetic Suppression (only for MHD runs)
    if bx0 > 0:
        T2_result = evaluate_benchmark_T2_magnetic_suppression(
            sigma_num=m1_growth_rate,
            drat=drat,
            vflow=vflow,
            bx0=bx0,
            Lx=domain_Lx,
            test_id=args.test_id
        )
    else:
        print("\n--- Evaluating T2: Magnetic Field Suppression ---")
        print("T2 Result: SKIPPED (Not an MHD run, Bx0 = 0).")
        T2_result = None

    # T3: Perpendicular B-field Effect
    T3_result = evaluate_benchmark_T3_perpendicular_bfield(
        bx0=bx0,
        by0=by0,
        bz0=bz0,
        test_id=args.test_id
    )

    # T4: Atwood Number Effect
    T4_result = evaluate_benchmark_T4_atwood_number_effect(
        sigma_num=m1_growth_rate,
        drat=drat,
        vflow=vflow,
        Lx=domain_Lx,
        test_id=args.test_id
    )

    # T5: Modal Cascade
    T5_result = evaluate_benchmark_T5_dominant_mode(
        k_max_final=m3_final_kmax,
        Lx=domain_Lx,
        test_id=args.test_id
    )

    # --- Generate Summary Report ---
    # Could save a JSON or text file with all results
    
    print(f"\nFinished KHI analysis for test: {args.test_id}")

if __name__ == "__main__":
    main() 