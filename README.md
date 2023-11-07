# MiniProject1
Lebwohl-Lasher Model of Liquid Crystals for UoBristol SCIFM0004 - Accelerated Computing

## BC4 Module selection
Code containing MPI:
```
module add languages/anaconda3/2020-3.8.5
```
All other code:
```
module add languages/anaconda3/2022.11-3.9.13
```

## BC4 Submissions
Example submission files are given for all methods
- ```submit.sh``` can be used for all files within the 'code' directory, excluding MPI
- ```submit_multiple_runs_MPI.sh``` shows an example MPI run, along with looped runs.
- Cython and Cython w/ MPI have their own .sh files within their respective directories within code.

## Running code from command line
- cProfile outputs for runtimes:
  - python -m cProfile [FILENAME] [PARAMETERS] > cProfile_Outputs/OUTPUT_FILENAME.txt
- Cython build:
  - Move to cython or cython_mpi directories within code
  - python [SETUP_FILENAME] build_ext --inplace
- Cython run:
  - mpiexec -n [NUM_WORKERS] python [RUN_FILENAME] [PARAMETERS]

## Testing
Within the root directory ```pytest``` can be run to test all available functions, however Cython is excluded.

## Saving Timings
saving_time.py is a script for scraping all timings from slurm output files and saving them in a text file while also outputting the mean, max, min and sum of all timings within the run.
```
python saving_time.py <input_file.txt> <saving_filename>
```
