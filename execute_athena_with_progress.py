import subprocess
import sys
import time # Added for a brief pause
import re   # Added for extracting log file path
from tqdm import tqdm

def extract_log_path(command_string):
    # Regex to find the log file path, typically after '>' and before '2>&1'
    # Example: ... > ./A0.1/A0.1.log 2>&1
    match = re.search(r'>\s*([^\s]+)\s*2>&1', command_string)
    if match:
        return match.group(1)
    return None

def main():
    commands = [line.strip() for line in sys.stdin if line.strip()]

    if not commands:
        print("No commands received to execute.")
        return

    print(f"Found {len(commands)} simulation commands to execute.")
    
    try:
        with tqdm(total=len(commands), unit="sim", desc="Overall Sweep Progress") as pbar_overall:
            for i, cmd_str in enumerate(commands):
                run_label = f"Run {i+1}/{len(commands)}"
                pbar_overall.set_description(f"{run_label}")
                print(f"\n--- Starting {run_label}: {cmd_str} ---")

                log_file_path = extract_log_path(cmd_str)
                if not log_file_path:
                    print(f"Warning: Could not extract log file path from command: {cmd_str}")
                    print("Skipping log tailing for this command.")
                
                tail_process = None
                try:
                    # Start Athena simulation
                    # stdout/stderr are captured by redirection in cmd_str, 
                    # but Popen's pipes are not used here to avoid blocking or complex async handling.
                    athena_process = subprocess.Popen(cmd_str, shell=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    if log_file_path:
                        # Wait a very short moment for the log file to be created
                        time.sleep(0.5) 
                        print(f"Tailing log file: {log_file_path} (Ctrl+C for this log, then again for sweep if needed)")
                        # Start tailing the log file
                        # Using sys.executable for python interpreter to run tail command
                        tail_cmd = [sys.executable, "-c", 
                                    f"import subprocess, sys; p = subprocess.Popen(['tail', '-f', '{log_file_path}'], stdout=sys.stdout, stderr=sys.stderr); p.wait()"]
                        # Tail process will output directly to the console.
                        tail_process = subprocess.Popen(tail_cmd) 

                    # Wait for Athena simulation to complete
                    athena_process.wait()

                    if athena_process.returncode != 0:
                        print(f"\nError during {run_label}. Athena exited with code: {athena_process.returncode}")
                        print(f"Command was: {cmd_str}")
                        print(f"Check log file {log_file_path if log_file_path else '(unknown)'} for details.")
                    else:
                        print(f"\n{run_label} completed successfully.")

                except Exception as e:
                    print(f"\nFailed to execute or monitor {run_label}: {cmd_str}\nError: {e}")
                finally:
                    if tail_process and tail_process.poll() is None: # Check if tail process is still running
                        print(f"Stopping log tail for {log_file_path}...")
                        tail_process.terminate() # Terminate tail
                        try:
                            tail_process.wait(timeout=2) # Wait a bit for it to close
                        except subprocess.TimeoutExpired:
                            tail_process.kill() # Force kill if it doesn't terminate
                    print(f"--- Finished {run_label} ---")
                
                pbar_overall.update(1)
    except ImportError:
        print("Error: tqdm library not found. Please install it using 'pip install tqdm'")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nSweep interrupted by user (KeyboardInterrupt).")
        # Ensure any running tail process is also stopped if the main loop is interrupted
        if tail_process and tail_process.poll() is None:
            print(f"Stopping active log tail due to sweep interruption...")
            tail_process.terminate()
            try: tail_process.wait(timeout=1)
            except subprocess.TimeoutExpired: tail_process.kill()
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during the sweep: {e}")
        sys.exit(1)

    print("\nAll simulation commands processed.")

if __name__ == "__main__":
    main() 