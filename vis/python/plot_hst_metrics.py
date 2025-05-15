import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import re # Import regex module

def plot_energy_drift(hst_filepath):
    """
    Reads an Athena++ history (.hst) file, calculates and plots the 
    relative total energy drift over time.

    Args:
        hst_filepath (str): Path to the .hst history file.
    """
    if not os.path.exists(hst_filepath):
        print(f"Error: History file not found at {hst_filepath}")
        return

    try:
        # Read the HST file. Athena++ .hst files are space-delimited.
        # The first line is a header starting with #, second line defines columns.
        with open(hst_filepath, 'r') as f:
            lines = f.readlines()
        
        header_line = ""
        column_names_line = ""
        data_lines_start_idx = 0
        for i, line in enumerate(lines):
            if line.strip().startswith("# Athena++ history data"):
                header_line = lines[i+1] # Next line should be column names
                data_lines_start_idx = i+2
                break
        
        if not header_line:
            print(f"Error: Could not find column definition header in {hst_filepath}")
            return

        # Extract column names using regex to be more robust
        # Pattern: find occurrences of e.g., "[1]=time", "[12]=whatever-name"
        matches = re.findall(r'\[\d+\]=([^\s\[]+)', header_line)
        col_names = [match.strip() for match in matches]

        if not col_names:
            print(f"Error: No column names extracted from header: {header_line}")
            # Fallback to trying to read without explicit column names, hoping pandas infers correctly
            # Or provide a default set if the number of columns is known/standard
            # For now, let's proceed and see if the fallback logic below handles it, or error out if too few cols.
            # This path should ideally not be taken if the header is standard.

        # Debug print for extracted column names
        print(f"Extracted column names: {col_names}")

        # Load data using pandas, skipping header rows.
        # We will assign column names after loading if extraction was successful.
        df = pd.read_csv(hst_filepath, sep='\s+', skiprows=data_lines_start_idx, header=None)

        if col_names and len(df.columns) == len(col_names):
            df.columns = col_names
            print(f"Assigned column names: {list(df.columns)}")
        elif col_names and len(df.columns) != len(col_names):
            print(f"Warning: Number of columns in data ({len(df.columns)}) does not match extracted names ({len(col_names)}).")
            print(f"  Data columns: {list(df.columns)}")
            print(f"  Extracted names: {col_names}")
            # Attempt to assign to a subset if possible, or use fallback
            if len(df.columns) < len(col_names):
                df.columns = col_names[:len(df.columns)]
                print(f"Assigned subset of names: {list(df.columns)}")
            else: # More data columns than names, use fallback
                print("Proceeding with fallback column naming.")
                col_names = [] # Force fallback
        
        if not col_names: # If col_names is empty (extraction failed or mismatch led to reset)
            # Fallback if column name extraction is tricky or failed
            print(f"Using fallback logic for column names. Number of data columns: {len(df.columns)}")
            base_cols = ['time', 'dt', 'mass', '1-mom', '2-mom', '3-mom', '1-KE', '2-KE', '3-KE', 'tot-E']
            if len(df.columns) >= 10:
                temp_cols = base_cols[:min(len(base_cols), len(df.columns))]
                if len(df.columns) >= 13: # Likely has ME
                    temp_cols.extend(['1-ME','2-ME','3-ME'])
                # Add dummy names for any extra cols
                temp_cols.extend([f'col{i+len(temp_cols)}' for i in range(len(df.columns) - len(temp_cols))])
                df.columns = temp_cols[:len(df.columns)]
                col_names = list(df.columns)
                print(f"Using inferred column names: {col_names}")
            else:
                print("Error: Cannot reliably assign column names.")
                return

        # Ensure numeric types
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df.dropna(inplace=True) # Drop rows that couldn't be cast to numeric

        if 'time' not in df.columns or 'tot-E' not in df.columns:
            print("Error: Essential columns 'time' or 'tot-E' not found.")
            return

        time_col = df['time']
        
        # Calculate total kinetic energy if components are present
        ke_cols = ['1-KE', '2-KE', '3-KE']
        if all(col in df.columns for col in ke_cols):
            total_ke = df[ke_cols].sum(axis=1)
        else: # Fallback: if individual KEs not found, assume tot-E might be just thermal if KE isn't included by default or if it's a very simple setup
            print("Warning: KE components not found, assuming 'tot-E' includes kinetic energy or KE is negligible for this context.")
            # This case is tricky without knowing exactly what tot-E contains. Usually it is E_int + E_kin.

        # Calculate total magnetic energy if components are present
        me_cols = ['1-ME', '2-ME', '3-ME']
        total_me = pd.Series(0.0, index=df.index) # Initialize with zeros
        if all(col in df.columns for col in me_cols):
            print("Found Magnetic Energy columns.")
            total_me = df[me_cols].sum(axis=1)
        else:
            print("Magnetic Energy columns (1-ME, 2-ME, 3-ME) not found. Assuming hydro run or ME is zero.")

        # tot-E from Athena usually means thermal + kinetic
        # Total conserved energy E_total = tot-E (thermal+kinetic) + magnetic_energy
        conserved_energy = df['tot-E'] + total_me

        if conserved_energy.empty:
            print("Error: Conserved energy array is empty after calculations.")
            return

        initial_total_energy = conserved_energy.iloc[0]
        if initial_total_energy == 0:
            print("Warning: Initial total energy is zero. Cannot compute relative drift. Plotting absolute drift.")
            energy_drift = conserved_energy - initial_total_energy
            drift_label = 'Absolute Energy Drift'
        else:
            energy_drift = (conserved_energy - initial_total_energy) / initial_total_energy
            drift_label = 'Relative Energy Drift (ΔE/E₀)'

        # --- Generate Plot ---
        plt.figure(figsize=(10, 6))
        plt.plot(time_col, energy_drift)
        
        plt.xlabel('Time (t)')
        plt.ylabel(drift_label)
        plt.title('Total Energy Conservation Over Time')
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()

        # --- Save Plot ---
        # Save in a 'plots' subdirectory relative to the hst file's location.
        # The hst file is already in the specific test's output directory.
        hst_dir = os.path.dirname(hst_filepath)
        filename_no_ext = os.path.splitext(os.path.basename(hst_filepath))[0]
        plot_dir = os.path.join(hst_dir, "plots") # Consistent with plot_kh_snapshot.py
        os.makedirs(plot_dir, exist_ok=True)
        
        output_plot_path = os.path.join(plot_dir, f"{filename_no_ext}_energy_drift.png")
        plt.savefig(output_plot_path)
        print(f"Energy drift plot saved to {output_plot_path}")
        plt.close()

    except FileNotFoundError:
        print(f"Error: HST file not found at {hst_filepath}")
    except Exception as e:
        print(f"An error occurred while processing {hst_filepath}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot energy drift from an Athena++ history file.')
    parser.add_argument('hst_file', type=str, help='Path to the .hst history file.')
    parser.add_argument('--delim_whitespace', action='store_true', help='Use delim_whitespace for pd.read_csv (deprecated)') # For testing old behavior
    
    args = parser.parse_args()
    plot_energy_drift(args.hst_file) 