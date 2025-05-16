#!/bin/bash

# Script to run a parameter sweep for Kelvin-Helmholtz simulations
# Assumes 'athena' executable is in the PATH or ./bin/athena
# Assumes 'sweep_config.json' and 'inputs/kh_sweep_template.in' exist
# Uses 'execute_athena_with_progress.py' for running simulations with tqdm

ATHENA_EXEC="./bin/athena"
if ! [ -f "${ATHENA_EXEC}" ]; then
    ATHENA_EXEC=$(which athena) # Try to find athena in PATH
fi

if ! [ -x "${ATHENA_EXEC}" ]; then
    echo "Error: Athena executable not found or not executable at ${ATHENA_EXEC}" >&2
    echo "Please ensure Athena is compiled and accessible." >&2
    exit 1
fi

CONFIG_FILE="sweep_config.json"
TEMPLATE_INPUT_FILE="inputs/kh_sweep_template.in"
PYTHON_EXEC_SCRIPT="execute_athena_with_progress.py"

if ! [ -f "${CONFIG_FILE}" ]; then
    echo "Error: Sweep configuration file '${CONFIG_FILE}' not found." >&2
    exit 1
fi

if ! [ -f "${TEMPLATE_INPUT_FILE}" ]; then
    echo "Error: Template input file '${TEMPLATE_INPUT_FILE}' not found." >&2
    exit 1
fi

if ! [ -f "${PYTHON_EXEC_SCRIPT}" ]; then
    echo "Error: Python execution script '${PYTHON_EXEC_SCRIPT}' not found." >&2
    exit 1
fi

# Check for Python3
if ! command -v python3 &> /dev/null
then
    echo "Error: python3 is not installed. Please install Python 3." >&2
    exit 1
fi

# Check for jq for parsing JSON
if ! command -v jq &> /dev/null
then
    echo "Error: jq is not installed. Please install jq to parse the config file." >&2
    exit 1
fi

# Check for bc for calculations
if ! command -v bc &> /dev/null
then
    echo "Error: bc is not installed. Please install bc for calculations." >&2
    exit 1
fi

echo "Starting KHI Simulation Sweep..."

# Read the config file and loop through runs
num_runs=$(jq '.runs | length' ${CONFIG_FILE})
declare -a ATHENA_COMMANDS

for i in $(seq 0 $((${num_runs} - 1)))
do
    run_id=$(jq -r ".runs[${i}].id" ${CONFIG_FILE})
    athdf_dir=$(jq -r ".runs[${i}].athdf_dir" ${CONFIG_FILE})
    atwood_number=$(jq -r ".runs[${i}].atwood_number" ${CONFIG_FILE})

    echo ""
    echo "-----------------------------------------------------"
    echo "Preparing run: ${run_id}"
    echo "Output directory: ${athdf_dir}"
    echo "Atwood Number: ${atwood_number}"

    # Calculate drat from Atwood number
    # A = (rho_h - rho_l) / (rho_h + rho_l)
    # Let rho_l = 1.0. Then rho_h = drat.
    # A = (drat - 1) / (drat + 1)
    # A * (drat + 1) = drat - 1
    # A * drat + A = drat - 1
    # A + 1 = drat - A * drat
    # A + 1 = drat * (1 - A)
    # drat = (1 + A) / (1 - A)
    # Need bc for calculations with floating point
    drat=$(bc -l <<< "(1 + ${atwood_number}) / (1 - ${atwood_number})")
    echo "Calculated density ratio (drat): ${drat}"

    # Create output directory if it doesn't exist
    mkdir -p "${athdf_dir}"

    # Create a temporary input file for this run
    temp_input_file="${athdf_dir}/kh_${run_id}.in"
    cp "${TEMPLATE_INPUT_FILE}" "${temp_input_file}"

    # Modify the problem_id and drat in the temporary input file
    # Using | as sed delimiter because paths might contain /
    sed -i.bak "s|problem_id = KHSweep|problem_id = ${run_id}|" "${temp_input_file}"
    sed -i.bak "s|drat       = 1.0|drat       = ${drat}|" "${temp_input_file}"
    rm "${temp_input_file}.bak" # Remove sed backup file

    # Construct the Athena command
    # Output and error streams are redirected to a log file in the run's directory
    # The '-d' flag specifies the output directory for .athdf files.
    command_string=""${ATHENA_EXEC}" -i "${temp_input_file}" -d "${athdf_dir}/" > "${athdf_dir}/${run_id}.log" 2>&1"
    ATHENA_COMMANDS+=("${command_string}")
    
    echo "Command prepared: ${command_string}"
done

echo ""
echo "-----------------------------------------------------"
echo "All simulation configurations prepared."
echo "Executing simulations with progress bar..."

# Execute all commands using the Python script with tqdm
printf "%s\n" "${ATHENA_COMMANDS[@]}" | python3 "${PYTHON_EXEC_SCRIPT}"

# Check the exit status of the Python script (optional, as Python script handles its own errors)
# python_exit_status=$?
# if [ ${python_exit_status} -ne 0 ]; then
#     echo "Error: Python execution script failed with exit code ${python_exit_status}." >&2
#     # exit ${python_exit_status} # Decide if sweep script should also exit
# fi

echo ""
echo "-----------------------------------------------------"
echo "KHI Simulation Sweep Finished."
echo "Please check individual log files in each run's directory (${athdf_dir}/${run_id}.log) for details."

# Make the script executable
# chmod +x run_kh_sweep.sh 