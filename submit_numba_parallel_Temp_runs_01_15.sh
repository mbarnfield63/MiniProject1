#!/bin/bash

#SBATCH --job-name=MiniProject1-LL
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:0:0
#SBATCH --mem-per-cpu=100M
#SBATCH --account=PHYS030544

# Load anaconda environment
module add languages/anaconda3/2022.11-3.9.13 

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

# File run
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.1 0 > cProfile_Outputs/Numba/numba_parallel_01.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.2 0 > cProfile_Outputs/Numba/numba_parallel_02.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.3 0 > cProfile_Outputs/Numba/numba_parallel_03.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.4 0 > cProfile_Outputs/Numba/numba_parallel_04.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.5 0 > cProfile_Outputs/Numba/numba_parallel_05.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.6 0 > cProfile_Outputs/Numba/numba_parallel_06.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.7 0 > cProfile_Outputs/Numba/numba_parallel_07.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.8 0 > cProfile_Outputs/Numba/numba_parallel_08.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 0.9 0 > cProfile_Outputs/Numba/numba_parallel_09.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 1.0 0 > cProfile_Outputs/Numba/numba_parallel_10.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 1.1 0 > cProfile_Outputs/Numba/numba_parallel_11.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 1.2 0 > cProfile_Outputs/Numba/numba_parallel_12.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 1.3 0 > cProfile_Outputs/Numba/numba_parallel_13.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 1.4 0 > cProfile_Outputs/Numba/numba_parallel_14.txt
python -m cProfile LebwohlLasher_numba_parallel.py 50 50 1.5 0 > cProfile_Outputs/Numba/numba_parallel_15.txt

# End recording the end time
end_time=$(date +%s)

# Calculate and print the runtime
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"