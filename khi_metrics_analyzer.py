import argparse
import json
import os
import re
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import math
import subprocess # Added for ffmpeg

# Helper function to plot 2D data
def plot_2d_field(field_data, x_coords, y_coords, title, filepath, xlabel="X1", ylabel="X2"):
    if field_data is None or np.all(np.isnan(field_data)):
        print(f"Skipping plot for {title} as data is None or all NaN.")
        return
    plt.figure(figsize=(8, 6))
    plt.imshow(field_data.T, origin='lower', aspect='auto',
               extent=[x_coords[0], x_coords[-1], y_coords[0], y_coords[-1]],
               cmap='viridis') # Or another cmap like 'RdBu_r' for divergent data
    plt.colorbar(label='Value')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(filepath)
    plt.close()
    print(f"Saved 2D field plot: {filepath}")

# New helper function for multi-panel snapshots
def plot_multipanel_snapshot(fields_data, x_coords, y_coords, panel_titles, main_title, filepath, xlabel="X1", ylabel="X2"):
    if len(fields_data) != 4 or len(panel_titles) != 4:
        print("Error: plot_multipanel_snapshot expects 4 fields and 4 panel_titles.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(17, 14), sharex=True, sharey=True) # Increased figsize
    fig.suptitle(main_title, fontsize=18, y=0.98) # Adjusted y for suptitle

    # Define specific colormaps for each panel for better visualization
    # vel1: shear, could use a diverging map if centered around 0, or sequential
    # vel2: KHI modes, diverging map is good
    # vorticity_z: diverging map
    # scalar: sequential map
    cmaps = ['viridis', 'RdBu_r', 'PiYG_r', 'coolwarm'] 
    
    # Determine common color limits for vel1 and vel2 if desired, or individual
    # For now, individual limits based on each field's data will be used by imshow by default

    for i, ax_row in enumerate(axes):
        for j, ax in enumerate(ax_row):
            idx = i * 2 + j
            field = fields_data[idx]
            title = panel_titles[idx]
            cmap = cmaps[idx]

            if field is None or np.all(np.isnan(field)):
                ax.text(0.5, 0.5, f"{title}\nData is None or all NaN",
                        ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.set_title(title, fontsize=14)
            else:
                # Transpose field for imshow
                im = ax.imshow(field.T, origin='lower', aspect='auto',
                               extent=[x_coords[0], x_coords[-1], y_coords[0], y_coords[-1]],
                               cmap=cmap)
                fig.colorbar(im, ax=ax, label='Value', fraction=0.046, pad=0.04) # Adjust colorbar size
                ax.set_title(title, fontsize=14)

            if i == 1:  # Bottom row
                ax.set_xlabel(xlabel, fontsize=12)
            if j == 0:  # Left column
                ax.set_ylabel(ylabel, fontsize=12)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust rect to prevent suptitle overlap
    plt.savefig(filepath)
    plt.close(fig)
    print(f"Saved multi-panel snapshot: {filepath}")

# Helper function to extract time from HDF5 filename or attributes
def extract_time_from_filename(filepath):
    try:
        with h5py.File(filepath, 'r') as f:
            if 'Time' in f.attrs:
                return float(f.attrs['Time'])
    except Exception:
        pass # Fallback to filename parsing
    
    # Match .out1.xxxxx.athdf or .out2.xxxxx.athdf
    match = re.search(r'\.out[12]\.(\d{5})\.athdf$', os.path.basename(filepath))
    if match:
        # Attempt to use the number as a sequence, but prefer Time attribute if available
        # For KHI, dt is 0.1, so file 00001 is t=0.1, 00010 is t=1.0.
        # This might not be robust if dt varies.
        # Let's assume the file number sequence directly maps to time steps if 'Time' attr fails
        # This part might need adjustment if file numbers are not directly proportional to time
        return int(match.group(1)) * 0.1 # Assuming 00000 is t=0, 00001 is t=0.1 for nominal time
    return 0

def calculate_kinetic_energies(rho, vel1, vel2, dV):
    ke_x = 0.5 * np.sum(rho * vel1**2 * dV)
    ke_y = 0.5 * np.sum(rho * vel2**2 * dV)
    return ke_x, ke_y

def calculate_vorticity_z(vel1, vel2, x1v, x2v):
    dv2_dx1 = np.gradient(vel2, x1v, axis=1)
    dv1_dx2 = np.gradient(vel1, x2v, axis=0)
    return dv2_dx1 - dv1_dx2

def calculate_vorticity_max(vort_z_field):
    return np.max(np.abs(vort_z_field))

def calculate_vortex_spacing(vel2_slice, x1f):
    nx = len(vel2_slice)
    if nx == 0:
        return np.nan
    Lx = x1f[-1] - x1f[0]
    freq = fftfreq(nx, d=Lx/nx) 
    wavenumbers = 2 * np.pi * freq
    yf = fft(vel2_slice)
    power_spectrum = np.abs(yf[1:nx//2])**2 
    positive_wavenumbers = wavenumbers[1:nx//2]
    if len(power_spectrum) == 0:
        return np.nan
    peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum)*0.1)
    if len(peaks) == 0:
        if np.max(power_spectrum) > 1e-9:
             dominant_k_index = np.argmax(power_spectrum) 
             dominant_k = positive_wavenumbers[dominant_k_index]
        else:
            return np.nan
    else:
        dominant_k = positive_wavenumbers[peaks[0]]
    if dominant_k == 0: return np.inf
    return 2 * np.pi / dominant_k

def calculate_mixing_variance(scalar_field):
    if scalar_field is None or np.all(np.isnan(scalar_field)):
        return np.nan
    return max(0.0, np.var(scalar_field))

def linear_fit_func(t, slope, intercept):
    return slope * t + intercept

def calculate_gamma_growth_rate(times, ke_y_series, fit_t_max=None):
    if len(times) < 2 or len(ke_y_series) < 2 or np.all(np.isnan(ke_y_series)):
        return np.nan, (np.nan, np.nan), None, None
    valid_mask = (ke_y_series > 1e-12) & (~np.isnan(ke_y_series))
    if np.sum(valid_mask) < 2:
        return np.nan, (np.nan, np.nan), None, None
    log_ke_y = np.log(ke_y_series[valid_mask])
    fit_times = times[valid_mask]
    if len(fit_times) < 2:
        return np.nan, (np.nan, np.nan), None, None
    if fit_t_max is None:
        if len(fit_times) > 3:
            diff_log_ke_y = np.diff(log_ke_y)
            if len(diff_log_ke_y) > 1:
                peak_growth_idx = np.argmax(diff_log_ke_y)
                max_ke_y_overall = np.max(ke_y_series[~np.isnan(ke_y_series)]) if np.any(~np.isnan(ke_y_series)) else 0
                current_ke_y_at_peak_growth = ke_y_series[valid_mask][peak_growth_idx] if max_ke_y_overall > 0 else 0
                if current_ke_y_at_peak_growth < 0.1 * max_ke_y_overall:
                    fit_end_idx = min(len(fit_times) -1, peak_growth_idx + int(len(fit_times)*0.3))
                else:
                    fit_end_idx = peak_growth_idx +1 
                fit_end_idx = max(1, fit_end_idx)
                fit_indices = slice(0, fit_end_idx + 1)
            else:
                 fit_indices = slice(0, len(fit_times))
        else:
            fit_indices = slice(0, len(fit_times))
    else:
        fit_indices = fit_times <= fit_t_max
        if np.sum(fit_indices) < 2:
             fit_indices = slice(0, len(fit_times))
    times_to_fit = fit_times[fit_indices]
    log_ke_y_to_fit = log_ke_y[fit_indices]
    if len(times_to_fit) < 2:
        return np.nan, (np.nan, np.nan), None, None
    try:
        params, covariance = curve_fit(linear_fit_func, times_to_fit, log_ke_y_to_fit)
        gamma = params[0] / 2.0
        return gamma, (params[0], params[1]), times_to_fit, log_ke_y_to_fit
    except RuntimeError:
        return np.nan, (np.nan, np.nan), None, None

def process_single_run(run_config, base_output_dir):
    run_id = run_config['id']
    athdf_dir = run_config['athdf_dir']
    label = run_config['label']
    print(f"\nProcessing run: {label} (ID: {run_id}) from {athdf_dir}")

    # Create output directories for this run
    run_main_output_dir = os.path.join(base_output_dir, run_id)
    initial_conditions_dir = os.path.join(run_main_output_dir, "initial_conditions")
    panel_snapshots_dir = os.path.join(run_main_output_dir, "panel_snapshots")
    os.makedirs(initial_conditions_dir, exist_ok=True)
    os.makedirs(panel_snapshots_dir, exist_ok=True)

    try:
        hdf5_files_raw = [os.path.join(athdf_dir, f) for f in os.listdir(athdf_dir) if f.endswith('.athdf') and f.startswith(run_id)]
    except FileNotFoundError:
        print(f"Error: Directory not found {athdf_dir} for run {run_id}. Skipping.")
        return None
    if not hdf5_files_raw:
        print(f"No HDF5 files found in {athdf_dir} for run {run_id}. Skipping.")
        return None

    file_time_pairs = sorted([(extract_time_from_filename(f_path), f_path) for f_path in hdf5_files_raw])
    sorted_hdf5_files = [f[1] for f in file_time_pairs]
    sim_times_nominal = np.array([f[0] for f in file_time_pairs])

    metrics_data = {
        'time': [], 'KE_x': [], 'KE_y': [], 'KE_total': [], 'vorticity_max': [], 
        'vortex_spacing': [], 'mixing_variance': []
    }

    num_files = len(sorted_hdf5_files)
    num_snapshots_target = 25 # Desired number of panel snapshots for the movie
    if num_files <= num_snapshots_target:
        snapshot_indices = list(range(num_files))
    else:
        snapshot_indices = np.linspace(0, num_files - 1, num_snapshots_target, dtype=int).tolist()
    
    # Ensure first and last frames are included if not already by linspace (should be)
    if 0 not in snapshot_indices and num_files > 0: snapshot_indices.insert(0,0)
    if (num_files -1) not in snapshot_indices and num_files > 0 : snapshot_indices.append(num_files-1)
    snapshot_indices = sorted(list(set(snapshot_indices)))

    frame_counter = 0 # For naming panel snapshot files sequentially

    for i, hdf5_file in enumerate(sorted_hdf5_files):
        try:
            with h5py.File(hdf5_file, 'r') as f:
                prim = f['prim'][()]
                # Debug: Print variable names and number of variables from attributes
                try:
                    var_names_attr = f.attrs.get('VariableNames')
                    if var_names_attr is None and 'prim' in f and hasattr(f['prim'], 'attrs'):
                        var_names_attr = f['prim'].attrs.get('VariableNames')
                    num_vars_attr = f.attrs.get('NumVariables')
                    print(f"Debug {hdf5_file}: HDF5 Attrs - NumVariables: {num_vars_attr}, VariableNames: {var_names_attr}")
                except Exception as e_attrs:
                    print(f"Debug {hdf5_file}: Could not read HDF5 attributes for var names: {e_attrs}")

                if prim.ndim != 5 or prim.shape[0] < 1:
                    print(f"Warning: 'prim' dataset in {hdf5_file} has unexpected shape: {prim.shape}. Expected 5D. Skipping file.")
                    # Fill with NaNs for this step
                    metrics_data['time'].append(sim_times_nominal[i] if i < len(sim_times_nominal) else np.nan)
                    for key in ['KE_x', 'KE_y', 'KE_total', 'vorticity_max', 'vortex_spacing', 'mixing_variance']:
                        metrics_data[key].append(np.nan)
                    continue
                
                rho = prim[0, 0, 0, :, :]
                # VariableNames were [b'rho' b'press' b'vel1' b'vel2' b'vel3' b'r0']
                # So, press is at index 1
                vel1 = prim[2, 0, 0, :, :] # vel1 at index 2
                vel2 = prim[3, 0, 0, :, :] # vel2 at index 3
                # vel3 would be at index 4, scalar at index 5

                num_prim_vars = prim.shape[0]
                nscalars_attr = f.attrs.get('NumScalars', 0)
                scal1 = np.full_like(rho, np.nan)
                print(f"Debug {hdf5_file}: prim.shape={prim.shape}, nscalars_attr={nscalars_attr}.")
                
                # Attempt to get variable names if they exist
                var_names_attr = f.attrs.get('VariableNames')
                if var_names_attr is None and 'prim' in f and hasattr(f['prim'], 'attrs'):
                    var_names_attr = f['prim'].attrs.get('VariableNames')

                # Convert to list of strings if not None
                variable_names = []
                if var_names_attr is not None:
                    try:
                        # Handle both list of bytes and list of strings
                        variable_names = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in var_names_attr]
                    except Exception as e:
                        print(f"Debug {hdf5_file}: Could not decode VariableNames: {e}")

                print(f"Debug {hdf5_file}: prim.shape={prim.shape}, nscalars_attr={nscalars_attr}, num_prim_vars={num_prim_vars}, VariableNames list: {variable_names}")

                scalar_idx = -1

                if nscalars_attr > 0:
                    idx_first_scalar = num_prim_vars - nscalars_attr
                    if idx_first_scalar >= 0 and idx_first_scalar < num_prim_vars:
                        scalar_idx = idx_first_scalar
                        print(f"Debug {hdf5_file}: Scalar index determined from nscalars_attr: {scalar_idx}")
                    else:
                        print(f"Warning {hdf5_file}: nscalars_attr={nscalars_attr} but derived idx_first_scalar={idx_first_scalar} is out of bounds for num_prim_vars={num_prim_vars}")
                
                if scalar_idx == -1 and variable_names: # If scalar not found via nscalars_attr, try var_names
                    # Common hydro variables (order can vary, this is a guess for typical Athena 2D)
                    # Common names for the first scalar are 's0', 'r0', 'scal0'
                    possible_scalar_names = ['s0', 'r0', 'scal0', 'scalar0'] 
                    # Check for pressure to estimate number of hydro variables more reliably
                    # Typical order: rho, press, vel1, vel2, vel3, (then scalars)
                    # Or: rho, vel1, vel2, vel3, press, (then scalars)
                    # Or: rho, vel1, vel2, vel3, E, (then scalars) - E is not press
                    
                    # Check for 'r0' specifically as seen in logs
                    if 'r0' in variable_names:
                        try:
                            scalar_idx = variable_names.index('r0')
                            print(f"Debug {hdf5_file}: Found 'r0' at index {scalar_idx} in VariableNames.")
                        except ValueError:
                            pass # Should not happen if 'r0' in variable_names
                    
                    if scalar_idx == -1: # If 'r0' not found or not primary check
                        for s_name in possible_scalar_names:
                            if s_name in variable_names:
                                try:
                                    potential_idx = variable_names.index(s_name)
                                    # Sanity check: is this index beyond typical hydro vars?
                                    # Assuming at least 4 hydro vars (rho, p, v1, v2) for 2D
                                    if potential_idx >= 4 and potential_idx < num_prim_vars : # ensure it's a valid index
                                        scalar_idx = potential_idx
                                        print(f"Debug {hdf5_file}: Found potential scalar '{s_name}' at index {scalar_idx} in VariableNames.")
                                        break 
                                except ValueError:
                                    continue
                
                if scalar_idx == -1: # Fallback if no explicit scalar info from attributes or names
                    # This case implies nscalars_attr was 0 or unreliable, and var_names didn't help.
                    # If prim has 5 or 6 vars for 2D (rho, vel1, vel2, vel3, press/E, [scal1]), guess scalar is last.
                    # For 2D, NHYDRO is typically 4 (rho, vx, vy, P) or 5 (rho, vx, vy, vz, P) if vel3 is there.
                    # The observed prim.shape[0] is 6. Hydro vars are rho, press, vel1, vel2, vel3. So scalar is index 5.
                    if num_prim_vars == (4 + 1 + 1) : # Assuming rho, p, v1, v2, (v3=0), s0
                         scalar_idx = 5 
                         print(f"Debug {hdf5_file}: Fallback: Assuming scalar is at index 5 for prim.shape[0]={num_prim_vars} (rho,p,v1,v2,v3,s0)")
                    elif num_prim_vars == (3 + 1 + 1): # Assuming rho, v1, v2, (v3=0), p, s0  (if pressure is 4th hydro var)
                         scalar_idx = 4 # This needs careful check of var order
                         print(f"Debug {hdf5_file}: Fallback: Assuming scalar is at index 4 for prim.shape[0]={num_prim_vars} (rho,v1,v2,p,s0)")
                    else:
                        print(f"Warning {hdf5_file}: Could not reliably determine scalar index. num_prim_vars={num_prim_vars}, nscalars_attr={nscalars_attr}, VariableNames='{variable_names}'. Scalar analysis will be skipped.")

                if scalar_idx != -1:
                    try:
                        scal1 = prim[scalar_idx, 0, 0, :, :]
                        print(f"Debug {hdf5_file}: Extracted scal1 from prim[{scalar_idx}]. Min={np.nanmin(scal1) if not np.all(np.isnan(scal1)) else 'all_nan'}, Max={np.nanmax(scal1) if not np.all(np.isnan(scal1)) else 'all_nan'}")
                        if np.all(scal1 == scal1[0,0]) and not np.isnan(scal1[0,0]): # Check if it's a constant field
                             print(f"Warning {hdf5_file}: Scalar field at index {scalar_idx} appears to be constant (value: {scal1[0,0]}). This might indicate an issue if mixing is expected.")
                        elif np.all(np.isnan(scal1)):
                             print(f"Warning {hdf5_file}: Scalar field at index {scalar_idx} is all NaN after extraction.")

                    except IndexError:
                        print(f"Error {hdf5_file}: Index {scalar_idx} out of bounds for prim array with shape {prim.shape} when trying to extract scalar.")
                        scal1 = np.full_like(rho, np.nan) 
                else:
                    # This means scal1 remains all NaN
                    print(f"Debug {hdf5_file}: No scalar index found. scal1 will be all NaN.")
                
                # Final check for scal1 before use
                if np.all(np.isnan(scal1)):
                    print(f"Warning {hdf5_file}: scal1 is all NaN before calculating mix_var.")
                
                sim_time_attr = f.attrs.get('Time', sim_times_nominal[i])
                x1f_raw, x2f_raw = f['x1f'][()], f['x2f'][()]
                x1v_raw, x2v_raw = f['x1v'][()], f['x2v'][()]
                x1f, x2f = x1f_raw.squeeze(), x2f_raw.squeeze()
                x1v, x2v = x1v_raw.squeeze(), x2v_raw.squeeze()

                if not (x1f.ndim == 1 and x1f.shape[0] > 0 and x2f.ndim == 1 and x2f.shape[0] > 0 and 
                        x1v.ndim == 1 and x1v.shape[0] > 0 and x2v.ndim == 1 and x2v.shape[0] > 0):
                    print(f"Warning: Coordinate data in {hdf5_file} has unexpected shape/size after squeeze. Skipping file.")
                    metrics_data['time'].append(sim_time_attr)
                    for key in ['KE_x', 'KE_y', 'KE_total', 'vorticity_max', 'vortex_spacing', 'mixing_variance']:
                        metrics_data[key].append(np.nan)
                    continue
                
                dx1, dx2 = x1f[1:] - x1f[:-1], x2f[1:] - x2f[:-1]
                if len(dx1) == 0 or len(dx2) == 0:
                    print(f"Warning: dx1 or dx2 is empty in {hdf5_file}. Skipping file.")
                    metrics_data['time'].append(sim_time_attr)
                    for key in ['KE_x', 'KE_y', 'KE_total', 'vorticity_max', 'vortex_spacing', 'mixing_variance']:
                        metrics_data[key].append(np.nan)
                    continue
                
                Lx, Ly = x1f[-1] - x1f[0], x2f[-1] - x2f[0]
                nx1_cells, nx2_cells = len(x1v), len(x2v)
                dV_cell = (Lx / nx1_cells) * (Ly / nx2_cells)

                ke_x, ke_y = calculate_kinetic_energies(rho, vel1, vel2, dV_cell)
                vort_z_field = calculate_vorticity_z(vel1, vel2, x1v, x2v)
                vort_max = calculate_vorticity_max(vort_z_field)
                mid_y_idx = vel2.shape[0] // 2
                vel2_slice_at_interface = vel2[mid_y_idx, :]
                vortex_space = calculate_vortex_spacing(vel2_slice_at_interface, x1f)
                mix_var = calculate_mixing_variance(scal1)

                metrics_data['time'].append(sim_time_attr)
                metrics_data['KE_x'].append(ke_x)
                metrics_data['KE_y'].append(ke_y)
                metrics_data['KE_total'].append(ke_x + ke_y)
                metrics_data['vorticity_max'].append(vort_max)
                metrics_data['vortex_spacing'].append(vortex_space)
                metrics_data['mixing_variance'].append(mix_var)

                # Initial conditions specific processing (first file)
                if i == 0:
                    print(f"--- Initial Conditions for {run_id} (t={sim_time_attr:.4f}) ---")
                    initial_ke_y = 0.5 * np.sum(rho * vel2**2 * dV_cell)
                    print(f"Initial KE_y: {initial_ke_y:.6e}")
                    if np.any(np.abs(vel2) > 1e-9):
                        print(f"Initial vel2 contains non-zero values. Max(abs(vel2)) = {np.max(np.abs(vel2)):.3e}")
                    else:
                        print("Initial vel2 appears to be zero everywhere.")
                    plot_2d_field(vel1, x1v, x2v, f"Initial vel1 - {label} (t={sim_time_attr:.2f})", 
                                  os.path.join(initial_conditions_dir, f"initial_vel1_{run_id}.png"))
                    plot_2d_field(vel2, x1v, x2v, f"Initial vel2 - {label} (t={sim_time_attr:.2f})", 
                                  os.path.join(initial_conditions_dir, f"initial_vel2_{run_id}.png"))
                    plot_2d_field(scal1, x1v, x2v, f"Initial scal1 - {label} (t={sim_time_attr:.2f})", 
                                  os.path.join(initial_conditions_dir, f"initial_scal1_{run_id}.png"))
                    print("--------------------------------------------")
                
                # Generate multi-panel snapshot if current index is in snapshot_indices
                if i in snapshot_indices:
                    time_str_for_title = f"t={sim_time_attr:.2f}"
                    # Use frame_counter for sequential filenames for ffmpeg
                    frame_filename_part = f"frame_{frame_counter:04d}"
                    
                    panel_filepath = os.path.join(panel_snapshots_dir, f"panel_snapshot_{run_id}_{frame_filename_part}.png")
                    
                    fields_to_plot = [vel1, vel2, vort_z_field, scal1]
                    panel_titles = [
                        f"vel1 (Shear Vel)\n{time_str_for_title}", 
                        f"vel2 (Transverse Vel)\n{time_str_for_title}", 
                        f"Vorticity (z-comp)\n{time_str_for_title}", 
                        f"Scalar Concentration\n{time_str_for_title}"
                    ]
                    main_plot_title = f"KHI Snapshot - {label} - Time: {sim_time_attr:.3f} (Frame: {frame_counter})"
                    
                    plot_multipanel_snapshot(fields_to_plot, x1v, x2v, panel_titles, main_plot_title, panel_filepath)
                    frame_counter += 1

        except Exception as e:
            print(f"Error processing file {hdf5_file}: {e}")
            metrics_data['time'].append(sim_times_nominal[i] if i < len(sim_times_nominal) else np.nan)
            for key in ['KE_x', 'KE_y', 'KE_total', 'vorticity_max', 'vortex_spacing', 'mixing_variance']:
                 metrics_data[key].append(np.nan)
            continue
    
    metrics_df = pd.DataFrame(metrics_data)
    metrics_df = metrics_df.sort_values(by='time').reset_index(drop=True)
    actual_sim_times_sorted = metrics_df['time'].values

    gamma, fit_params, fitted_t, fitted_log_ke_y = calculate_gamma_growth_rate(actual_sim_times_sorted, metrics_df['KE_y'].values)
    metrics_df['gamma_growth_rate'] = gamma 

    csv_path = os.path.join(run_main_output_dir, 'metrics.csv')
    metrics_df.to_csv(csv_path, index=False)
    print(f"Metrics saved to {csv_path}")

    growth_fit_path = os.path.join(run_main_output_dir, 'growth_fit.txt')
    with open(growth_fit_path, 'w') as f_growth:
        f_growth.write(f"Run ID: {run_id} ({label})\n")
        f_growth.write(f"Calculated gamma_growth_rate: {gamma}\n")
        if fit_params and not np.any(np.isnan(fit_params)):
             f_growth.write(f"Fit parameters (slope, intercept) for log(KE_y) vs t: {fit_params[0]}, {fit_params[1]}\n")
        if fitted_t is not None and len(fitted_t) > 0:
            f_growth.write(f"Fitted over t_start={fitted_t[0]:.4f}, t_end={fitted_t[-1]:.4f}\n")
    print(f"Growth fit info saved to {growth_fit_path}")

    plt.style.use('seaborn-v0_8-whitegrid')
    if fitted_t is not None and fitted_log_ke_y is not None and not np.isnan(gamma) and fit_params and not np.any(np.isnan(fit_params)):
        plt.figure(figsize=(10, 6))
        plt.plot(actual_sim_times_sorted, np.log(np.maximum(metrics_df['KE_y'].values, 1e-16)), 'o', label='log(KE_y) data', alpha=0.7)
        plt.plot(fitted_t, fitted_log_ke_y, '-', color='red', linewidth=2, label=f'Fit: log(KE_y) ~ {fit_params[0]:.2f}*t + {fit_params[1]:.2f}\nGamma = {gamma:.3f}')
        plt.xlabel('Time'); plt.ylabel('log(KE_y)'); plt.title(f'Transverse Kinetic Energy Growth - {label}'); plt.legend();
        plt.savefig(os.path.join(run_main_output_dir, f'growth_fit_KEy_{run_id}.png'))
        plt.close()
        print(f"Growth fit plot saved to {os.path.join(run_main_output_dir, f'growth_fit_KEy_{run_id}.png')}")

    fig_summary, axes_summary = plt.subplots(5, 1, figsize=(12, 18), sharex=True)
    fig_summary.suptitle(f'KHI Metrics Summary - {label}', fontsize=16)
    axes_summary[0].plot(metrics_df['time'], metrics_df['KE_x'], label='KE_x')
    axes_summary[0].set_ylabel('Longitudinal KE'); axes_summary[0].legend()
    axes_summary[1].plot(metrics_df['time'], metrics_df['KE_y'], label='KE_y')
    axes_summary[1].set_ylabel('Transverse KE'); axes_summary[1].legend()
    axes_summary[2].plot(metrics_df['time'], metrics_df['KE_total'], label='KE_total')
    axes_summary[2].set_ylabel('Total KE'); axes_summary[2].legend()
    axes_summary[3].plot(metrics_df['time'], metrics_df['vorticity_max'], label='Max Vorticity (abs)')
    axes_summary[3].set_ylabel('Max |ω_z|'); axes_summary[3].legend()
    axes_summary[4].plot(metrics_df['time'], metrics_df['mixing_variance'], label='Mixing Var (scal1)')
    axes_summary[4].set_ylabel('Var(scal1)'); axes_summary[4].set_xlabel('Time'); axes_summary[4].legend()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    summary_plot_path = os.path.join(run_main_output_dir, f'summary_plot_{run_id}.png')
    plt.savefig(summary_plot_path); plt.close(fig_summary)
    print(f"Summary plot saved to {summary_plot_path}")

    if not metrics_df['vortex_spacing'].isnull().all():
        plt.figure(figsize=(10,6))
        plt.plot(metrics_df['time'], metrics_df['vortex_spacing'], marker='.', linestyle='-')
        plt.xlabel('Time'); plt.ylabel('Dominant Vortex Wavelength (λ_x)'); plt.title(f'Vortex Spacing Evolution - {label}')
        plt.ylim(bottom=0, top= (metrics_df['vortex_spacing'].max() * 1.2 if metrics_df['vortex_spacing'].notnull().any() else 1.0) ) # Ensure y starts at 0
        plt.savefig(os.path.join(run_main_output_dir, f'fft_vy_peak_lambda_{run_id}.png')); plt.close()
        print(f"Vortex spacing plot saved to {os.path.join(run_main_output_dir, f'fft_vy_peak_lambda_{run_id}.png')}")
    else:
        print(f"Skipping vortex spacing plot for {label} as data is all NaN.")

    return metrics_df

def main():
    parser = argparse.ArgumentParser(description="Analyze KHI simulation data based on a sweep configuration.")
    parser.add_argument("--config", type=str, default="sweep_config.json", help="Path to the JSON config file defining the simulation runs.")
    parser.add_argument("--output_dir", type=str, default="./khi_metrics_analysis_results", help="Base directory to save analysis results (CSVs, plots).")
    args = parser.parse_args()
    try:
        with open(args.config, 'r') as f: sweep_config = json.load(f)
    except FileNotFoundError: print(f"Error: Configuration file {args.config} not found."); return
    except json.JSONDecodeError: print(f"Error: Could not decode JSON from {args.config}."); return
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Analysis results will be saved in: {args.output_dir}")
    all_runs_metrics = {}
    for run_cfg in sweep_config.get('runs', []):
        run_metrics_df = process_single_run(run_cfg, args.output_dir)
        if run_metrics_df is not None: all_runs_metrics[run_cfg['id']] = run_metrics_df
    
    print("\n--- Generating Movies (requires ffmpeg) ---")
    ffmpeg_executable = "ffmpeg" # Assume ffmpeg is in PATH
    try:
        # Check if ffmpeg is available
        subprocess.run([ffmpeg_executable, '-version'], capture_output=True, check=True, text=True)
        print("ffmpeg found.")
        ffmpeg_available = True
    except (FileNotFoundError, subprocess.CalledProcessError) as e_ffmpeg_check:
        print(f"Warning: ffmpeg not found or not working. Movies will not be generated. Error: {e_ffmpeg_check}")
        ffmpeg_available = False

    if ffmpeg_available:
        for run_cfg in sweep_config.get('runs', []):
            run_id = run_cfg['id']
            label = run_cfg['label']
            current_run_output_dir = os.path.join(args.output_dir, run_id)
            panel_snapshots_dir_for_movie = os.path.join(current_run_output_dir, "panel_snapshots")
            movie_output_path = os.path.join(current_run_output_dir, f"movie_KHI_{run_id}.mp4")
            
            if not os.path.isdir(panel_snapshots_dir_for_movie):
                print(f"Panel snapshots directory not found for {label}: {panel_snapshots_dir_for_movie}. Skipping movie.")
                continue
            
            # Correctly list files using the full path pattern for glob or filtering
            # Files are named like: panel_snapshot_A0.1_frame_0000.png
            input_pattern_glob = f"panel_snapshot_{run_id}_frame_*.png"
            # Check if files exist before attempting ffmpeg
            existing_frames = sorted([f for f in os.listdir(panel_snapshots_dir_for_movie) if f.startswith(f"panel_snapshot_{run_id}_frame_") and f.endswith(".png")])

            if not existing_frames:
                print(f"No panel snapshot PNG files matching pattern found in {panel_snapshots_dir_for_movie} for {label}. Skipping movie.")
                continue
            
            # The input pattern for ffmpeg needs to be relative to its CWD, or use full paths.
            # If ffmpeg is run with cwd=panel_snapshots_dir_for_movie, then pattern is simpler.
            # Otherwise, provide full path to pattern.
            # Let's provide the full path pattern to ffmpeg's -i argument.
            input_ffmpeg_pattern = os.path.join(panel_snapshots_dir_for_movie, f"panel_snapshot_{run_id}_frame_%04d.png")

            ffmpeg_cmd = [
                ffmpeg_executable, '-y', 
                '-framerate', '10', 
                '-i', input_ffmpeg_pattern,
                '-c:v', 'libx264', 
                '-r', '20',         
                '-pix_fmt', 'yuv420p', 
                '-crf', '23', 
                movie_output_path
            ]
            
            try:
                print(f"Attempting to generate movie for {label}...")
                print(f"Command: {' '.join(ffmpeg_cmd)}") # For debugging the command
                
                # Running ffmpeg from the base directory, so paths in command should be correct.
                result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True, check=False) 
                
                if result.returncode == 0:
                    print(f"Movie successfully generated for {label}: {movie_output_path}")
                else:
                    print(f"Error generating movie for {label}. Return code: {result.returncode}")
                    print(f"FFmpeg stdout:\n{result.stdout}")
                    print(f"FFmpeg stderr:\n{result.stderr}")
            except Exception as e:
                print(f"An error occurred during movie generation for {label}: {e}")
    
    print("\nKHI Metrics Analysis Finished.")

if __name__ == "__main__":
    main() 