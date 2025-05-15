import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_snapshot(hdf5_filepath):
    """
    Reads density data from an Athena++ HDF5 file and generates a 2D pseudocolor plot.

    Args:
        hdf5_filepath (str): Path to the HDF5 output file.
    """
    if not os.path.exists(hdf5_filepath):
        print(f"Error: File not found at {hdf5_filepath}")
        return

    try:
        with h5py.File(hdf5_filepath, 'r') as f:
            # --- Read data ---
            # Assuming 'prim' (primitive variables) are saved and density is the first one (index 0)
            # The dataset shape is typically (num_variables, nx3, nx2, nx1) or (num_variables, nx2, nx1) for 2D
            
            prim_data = f['prim']
            dataset_shape = prim_data.shape
            print(f"Initial prim data shape: {dataset_shape}") # Debug print

            if len(dataset_shape) == 5: # Potentially (num_vars, num_blocks_z, num_blocks_y, nz_block, nx_block) or (num_vars, nz_glob, ny_glob, nx_block_DEPRECATED?, nx_glob) -> needs careful check
                # This handles cases where Athena++ might output with extra dimensions for block structure, even in serial.
                # Example: (nvar, 1, 1, nx2, nx1) for a 2D problem run as a single block.
                # Example: (nvar, 1, n_block_y, n_block_x, nz_cell_block, ny_cell_block, nx_cell_block)
                # For our current test cases (single block):
                # 2D case might be (nvar, 1, 1, nx2, nx1)
                # 3D case might be (nvar, 1, 1, nx3, nx2, nx1) - no, this would be 6D.
                # Let's assume for single block output it might be: (nvar, 1, 1, nz, ny, nx) if nz=1 for 2D.
                # Or (nvar, nz_total_actual, ny_total_actual, nx_total_actual) if those dimensions are squeezed out by h5py.
                # The error was (5, 1, 1, 128, 256) for a 256x128 2D run. Here 5=nvar.
                # So, prim_data[0, 0, 0, :, :] would get density.
                # For 3D (64x32x32) test: (5, 1, 32, 32, 64) -> prim_data[nvar_idx, 0, slice_idx_z, :, :]
                
                nvar_dim = dataset_shape[0]
                # Check if it's the 2D case that looks like (nvar, 1, 1, ny, nx)
                if dataset_shape[1] == 1 and dataset_shape[2] == 1 and len(dataset_shape) == 5: # Likely 2D data with extra block dims
                    density = prim_data[0, 0, 0, :, :] 
                    print(f"Extracted 2D density from 5D array (0,0,0,:,:). Density shape: {density.shape}")
                # Check if it's the 3D case that looks like (nvar, 1, nz, ny, nx) - this was previous error for 3D.
                # Error for 3D (64x32x32) was (5, 1, 32, 32, 64). This means (nvar, block_z=1, nz_cells, ny_cells, nx_cells)
                elif dataset_shape[1] == 1 and len(dataset_shape) == 5: # Likely 3D data with an extra block_z dim
                    slice_idx_z = dataset_shape[2] // 2 # Slice through the actual z-cells
                    density = prim_data[0, 0, slice_idx_z, :, :]
                    print(f"Extracted 2D slice from 5D array (3D data with block_z=1) (0,0,{slice_idx_z},:,:). Density shape: {density.shape}")
                else:
                    print(f"Error: Unhandled 5D data shape {dataset_shape}.")
                    return

            elif len(dataset_shape) == 4: # 3D data (num_vars, nz, ny, nx)
                # Take a slice at the mid-plane in z (k=0 if nz=1, or mid-plane if nz > 1)
                slice_idx_z = dataset_shape[1] // 2
                density = prim_data[0, slice_idx_z, :, :]
                print(f"Extracted 2D slice from 3D data (z-index: {slice_idx_z}). Density shape: {density.shape}")
            elif len(dataset_shape) == 3: # 2D data (num_vars, ny, nx)
                density = prim_data[0, :, :]
                print(f"Using 2D data. Density shape: {density.shape}")
            else:
                print(f"Error: Unexpected data shape {dataset_shape}. Expected 3 or 4 dimensions.")
                return

            # --- Read coordinates ---
            # x1f, x2f, x3f are face-centered coordinates
            # x1v, x2v, x3v are cell-centered coordinates (usually not directly in HDF5, derive from faces)
            x1f_data = f['x1f']
            x2f_data = f['x2f']
            print(f"Initial x1f shape: {x1f_data.shape}, x2f shape: {x2f_data.shape}")

            # Adjust for potential extra dimensions in coordinate arrays if they exist
            # Based on Athena++ HDF5 format for single-block output, x1f, x2f, x3f should be 1D.
            # However, if they are not, this attempts a simple squeeze or specific slice.
            if x1f_data.ndim > 1:
                # Try to squeeze out dimensions of size 1, or take a specific slice if a known pattern.
                # If shape is (1, 1, N+1), slice [0, 0, :]
                if x1f_data.shape[0] == 1 and x1f_data.shape[1] == 1 and x1f_data.ndim == 3:
                    x1f = x1f_data[0,0,:]
                    print(f"Sliced x1f from {x1f_data.shape} to {x1f.shape}")
                else: # General squeeze, might be too aggressive if other dims are not 1
                    x1f = np.squeeze(x1f_data)
                    print(f"Squeezed x1f from {x1f_data.shape} to {x1f.shape}")
            else:
                x1f = x1f_data[:]

            if x2f_data.ndim > 1:
                if x2f_data.shape[0] == 1 and x2f_data.shape[1] == 1 and x2f_data.ndim == 3:
                    x2f = x2f_data[0,0,:]
                    print(f"Sliced x2f from {x2f_data.shape} to {x2f.shape}")
                else:
                    x2f = np.squeeze(x2f_data)
                    print(f"Squeezed x2f from {x2f_data.shape} to {x2f.shape}")
            else:
                x2f = x2f_data[:]
            
            # Calculate cell-centered coordinates for plotting
            # These should now be 1D arrays if slicing/squeezing was correct
            if x1f.ndim == 1 and x1f.size > 1:
                x1v = 0.5 * (x1f[:-1] + x1f[1:])
            else:
                print(f"Error: x1f is not a 1D array with more than 1 element after processing. Shape: {x1f.shape}")
                x1v = np.array([]) # Empty array to ensure len(x1v) is 0 and triggers error below
            
            if x2f.ndim == 1 and x2f.size > 1:
                x2v = 0.5 * (x2f[:-1] + x2f[1:])
            else:
                print(f"Error: x2f is not a 1D array with more than 1 element after processing. Shape: {x2f.shape}")
                x2v = np.array([])

            print(f"Processed x1f shape: {x1f.shape}, x1v shape: {x1v.shape}")
            print(f"Processed x2f shape: {x2f.shape}, x2v shape: {x2v.shape}")

            if density.shape != (len(x2v), len(x1v)):
                 print(f"Warning: Density shape {density.shape} does not match derived cell-centered coordinate lengths ({len(x2v)}, {len(x1v)}). Adjusting...")
                 # This might happen if data is (nx1, nx2) vs (nx2, nx1) due to read/transpose.
                 # Or if slicing logic for 3D needs adjustment.
                 # A common convention is (x2_dim, x1_dim) for images with imshow
                 if density.shape == (len(x1v), len(x2v)):
                     density = density.T # Transpose if necessary
                     print(f"Transposed density to {density.shape} to match coordinates.")
                 else:
                     print("Error: Could not reconcile density shape with coordinate dimensions after transpose attempt.")
                     # Attempt to use file attributes for dimensions if available
                     nx1_attr = f.attrs.get('RootGridSize')[0]
                     nx2_attr = f.attrs.get('RootGridSize')[1]
                     print(f"Attributes: nx1_attr={nx1_attr}, nx2_attr={nx2_attr}") # Debug print for attributes

                     if density.shape == (nx2_attr, nx1_attr) and len(x2v) == nx2_attr and len(x1v) == nx1_attr:
                        print("Shapes match attributes RootGridSize. Proceeding.")
                     elif density.shape == (nx1_attr, nx2_attr) and len(x1v) == nx1_attr and len(x2v) == nx2_attr:
                        density = density.T
                        print(f"Transposed density based on attributes to {density.shape}.")
                     else:
                        print(f"Still mismatched. Density: {density.shape}, Expected based on x1v,x2v: ({len(x2v)},{len(x1v)}) or based on attrs: ({nx2_attr}, {nx1_attr})")
                        return


            # --- Read simulation time ---
            sim_time = f.attrs.get('Time', 0.0)

            # --- Generate Plot ---
            plt.figure(figsize=(10, 8))
            
            # Create a meshgrid for pcolormesh
            # pcolormesh expects X, Y to define the corners of the quadrilaterals.
            # So we should use the face-centered coordinates.
            X1f, X2f = np.meshgrid(x1f, x2f)

            # Ensure density array is (N_x2, N_x1) for pcolormesh if X1f, X2f are (N_x2+1, N_x1+1)
            # The density array should have dimensions one less than the coordinate arrays.
            # h5py reads as (nz, ny, nx) or (ny, nx) for density.
            # matplotlib pcolormesh with X, Y from meshgrid(x1f, x2f) expects C (density) to be (len(x2f)-1, len(x1f)-1)
            
            # Check if density dimensions match (len(x2f)-1, len(x1f)-1)
            expected_density_shape_for_plot = (len(x2f)-1, len(x1f)-1)
            if density.shape != expected_density_shape_for_plot:
                print(f"Adjusting density shape {density.shape} to expected {expected_density_shape_for_plot} for pcolormesh.")
                # This could be a transpose, or if we read cell-centered data but plot with face-coords
                # For now, assume it's correctly (N_x2_cells, N_x1_cells)
                if density.shape == (expected_density_shape_for_plot[1],expected_density_shape_for_plot[0]): # (N_x1_cells, N_x2_cells)
                    density = density.T
                    print(f"Transposed density for plotting to {density.shape}")
                else:
                    print("Error: Density shape cannot be easily reconciled for pcolormesh with face coordinates.")
                    print(f"Density shape: {density.shape}, x1f shape: {x1f.shape}, x2f shape: {x2f.shape}")
                    # Fallback to imshow with cell-centered coordinates if pcolormesh is problematic
                    print("Attempting fallback to imshow with cell-centered coordinates.")
                    X1v, X2v = np.meshgrid(x1v, x2v)
                    plt.imshow(density, aspect='auto', origin='lower',
                               extent=[x1v.min(), x1v.max(), x2v.min(), x2v.max()],
                               cmap='viridis')

            else:
                 # If shapes are fine, proceed with pcolormesh
                plt.pcolormesh(X1f, X2f, density, cmap='viridis', shading='auto')


            plt.colorbar(label='Density')
            plt.xlabel(f'x1 (Length Unit)')
            plt.ylabel(f'x2 (Length Unit)')
            plt.title(f'Density Snapshot at t = {sim_time:.2f}')
            plt.gca().set_aspect('equal', adjustable='box') # Make aspect ratio equal
            plt.tight_layout()

            # --- Save Plot ---
            base_dir = os.path.dirname(hdf5_filepath)
            filename_no_ext = os.path.splitext(os.path.basename(hdf5_filepath))[0]
            # Put plots in a 'plots' subdirectory relative to the HDF5 file's location
            plot_dir = os.path.join(base_dir, "plots")
            os.makedirs(plot_dir, exist_ok=True)
            
            output_plot_path = os.path.join(plot_dir, f"{filename_no_ext}_density.png")
            plt.savefig(output_plot_path)
            print(f"Plot saved to {output_plot_path}")
            plt.close()

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot density from an Athena++ HDF5 snapshot.')
    parser.add_argument('hdf5_file', type=str, help='Path to the HDF5 output file (e.g., kh_hydro_A033_256x128.out2.00000.athdf).')
    
    args = parser.parse_args()
    
    plot_snapshot(args.hdf5_file) 