#!/bin/bash

#SBATCH --job-name=MiniProject1-LL
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=6:0:0
#SBATCH --mem-per-cpu=100M
#SBATCH --account=PHYS030544

# Load anaconda environment
module add languages/anaconda3/2020-3.8.5

# Change to working directory, where job was submitted from
cd "${SLURM_SUBMIT_DIR}"

# Record some potentially useful details about the job: 
echo "Running on host $(hostname)"
echo "Started on $(date)"
echo "Directory is $(pwd)"
echo "Slurm job ID is ${SLURM_JOBID}"
echo "This jobs runs on the following machines:"
echo "${SLURM_JOB_NODELIST}"
echo "CPUs per task = ${SLURM_CPUS_PER_TASK}"
printf "\n\n"

# Submitting and timing code runs
# Recording start time
start_time=$(date +%s)

python setup_LebwohlLasher_mpi4cython.py build_ext --inplace

# Number of workers must be <= Lattice size
workers=16
echo "Workers = ${workers}"

x=10

for ((i=1; i<=$x; i++)); do
    mpiexec -n $workers python run_LebwohlLasher_mpi4cython.py 50 10 0.5
done

for ((i=1; i<=$x; i++)); do
    mpiexec -n $workers python run_LebwohlLasher_mpi4cython.py 50 50 0.5
done

for ((i=1; i<=$x; i++)); do
    mpiexec -n $workers python run_LebwohlLasher_mpi4cython.py 50 100 0.5
done

for ((i=1; i<=$x; i++)); do
    mpiexec -n $workers python run_LebwohlLasher_mpi4cython.py 50 250 0.5
done

for ((i=1; i<=$x; i++)); do
    mpiexec -n $workers python run_LebwohlLasher_mpi4cython.py 50 500 0.5
done

for ((i=1; i<=$x; i++)); do
    mpiexec -n $workers python run_LebwohlLasher_mpi4cython.py 50 1000 0.5
done

# End recording the end time
end_time=$(date +%s)

# Calculate and print the runtime
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"