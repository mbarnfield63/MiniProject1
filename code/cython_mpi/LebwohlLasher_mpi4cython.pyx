"""
Basic Python Lebwohl-Lasher code.  Based on the paper 
P.A. Lebwohl and G. Lasher, Phys. Rev. A, 6, 426-429 (1972).
This version in 2D.

Run at the command line by typing:

python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE>

where:
  ITERATIONS = number of Monte Carlo steps, where 1MCS is when each cell
      has attempted a change once on average (i.e. SIZE*SIZE attempts)
  SIZE = side length of square lattice
  TEMPERATURE = reduced temperature in range 0.0 - 2.0.
  
The initial configuration is set at random. The boundaries
are periodic throughout the simulation.  During the
time-stepping, an array containing two domains is used; these
domains alternate between old data and new data.

SH 16-Oct-23
"""

# Build using: python setup_LebwohlLasher_cython.py build_ext --inplace
# Run using: 

#cython: language_level=3
import sys
import time
import datetime
import numpy as np

from mpi4py import MPI
cimport numpy as np
from libc.math cimport cos, exp, pi

#=======================================================================
def initdat(int nmax):
    """
    Arguments:
      nmax (int) = size of lattice to create (nmax,nmax).
    Description:
      Function to create and initialise the main data array that holds
      the lattice.  Will return a square lattice (size nmax x nmax)
	  initialised with random orientations in the range [0,2pi].
	Returns:
	  arr (float(nmax,nmax)) = array to hold lattice.
    """
    # cdef:
    #   double[:, :] arr = np.empty((nmax, nmax))

    arr = np.random.random_sample((nmax,nmax))*2.0*np.pi
    return arr
#=======================================================================

# Plotdat removed due to errors being too intrusive when using
# module add languages/anaconda3/2020-3.8.5
# which is required for mpi4py

#=======================================================================
def savedat(double[:, ::1] arr, int nsteps, double Ts, double[:] ratio, double[:] energy, double[:] order, int nmax):
    # Create filename based on the current date and time.
    current_datetime = datetime.datetime.now().strftime("%a-%d-%b-%Y-at-%I-%M-%S%p")
    filename = "/user/home/eh19374/Year4/MiniProject1/Outputs/LL-Output-{:s}.txt".format(current_datetime).encode('utf-8')
    
    # Declare variables
    cdef double runtime = 0.0  # Assign a value to runtime
    cdef int i

    # Write the data to the file
    with open(filename, "wb") as FileOut:
        # Write the header with run parameters
        FileOut.write("#=====================================================\n".encode('utf-8'))
        FileOut.write("# File created:        {:s}\n".format(current_datetime).encode('utf-8'))
        FileOut.write("# Size of lattice:     {:d}x{:d}\n".format(nmax, nmax).encode('utf-8'))
        FileOut.write("# Number of MC steps:  {:d}\n".format(nsteps).encode('utf-8'))
        FileOut.write("# Reduced temperature: {:5.3f}\n".format(Ts).encode('utf-8'))
        FileOut.write("# Run time (s):        {:8.6f}\n".format(runtime).encode('utf-8'))
        FileOut.write("#=====================================================\n".encode('utf-8'))
        FileOut.write("# MC step:  Ratio:     Energy:   Order:\n".encode('utf-8'))
        FileOut.write("#=====================================================\n".encode('utf-8'))
        
        # Write the columns of data
        for i in range(nsteps + 1):
            FileOut.write("\n   {:05d}    {:6.4f} {:12.4f}  {:6.4f} ".format(i, ratio[i], energy[i], order[i]).encode('utf-8'))
#=======================================================================
cdef double one_energy(double[:, :] arr, int ix, int iy, int nmax):
    """
    Arguments:
      arr (double[:, ::1]) = array that contains lattice data;
      ix (int) = x lattice coordinate of cell;
      iy (int) = y lattice coordinate of cell;
      nmax (int) = side length of square lattice.
    Description:
      Function that computes the energy of a single cell of the
      lattice taking into account periodic boundaries. Working with
      reduced energy (U/epsilon), equivalent to setting epsilon=1 in
      equation (1) in the project notes.
    Returns:
      en (double) = reduced energy of cell.
    """
    cdef:
      double en = 0.0
      int ixp, ixm, iyp, iym

    en = 0.0
    ixp = (ix+1)%nmax # These are the coordinates
    ixm = (ix-1)%nmax # of the neighbours
    iyp = (iy+1)%nmax # with wraparound
    iym = (iy-1)%nmax #
#
# Add together the 4 neighbour contributions
# to the energy
#
    ang = arr[ix,iy]-arr[ixp,iy]
    en += 0.5*(1.0 - 3.0*cos(ang)**2)
    ang = arr[ix,iy]-arr[ixm,iy]
    en += 0.5*(1.0 - 3.0*cos(ang)**2)
    ang = arr[ix,iy]-arr[ix,iyp]
    en += 0.5*(1.0 - 3.0*cos(ang)**2)
    ang = arr[ix,iy]-arr[ix,iym]
    en += 0.5*(1.0 - 3.0*cos(ang)**2)
    return en

#=======================================================================
cdef double all_energy(double[:, :] arr, int nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
      nmax (int) = side length of square lattice.
    Description:
      Function to compute the energy of the entire lattice. Output
      is in reduced units (U/epsilon).
	Returns:
	  enall (float) = reduced energy of lattice.
    """
    cdef: 
      double enall = 0.0

    # Calculate the energy of the entire lattice    
    for i in range(nmax):
        for j in range(nmax):
            enall += one_energy(arr, i, j, nmax)

    return enall
#=======================================================================
cdef get_order(double[:,:] arr, int nmax):
    cdef:
        double[:,:] Qab = np.zeros((3, 3))
        double[:,:] delta = np.eye(3)
        double[:,:,:] lab = np.empty((3, nmax, nmax))
        double[:] eigenvalues
        int scalar = (2*nmax*nmax)

    lab = np.vstack((np.cos(arr),np.sin(arr),np.zeros_like(arr))).reshape(3,nmax,nmax)
    for a in range(3):
        for b in range(3):
            for i in range(nmax):
                for j in range(nmax):
                    Qab[a,b] += 3*lab[a,i,j]*lab[b,i,j] - delta[a,b]
            Qab[a,b] /= scalar

  
    # Use your own function to compute eigenvalues
    eigenvalues = np.linalg.eigvals(Qab)
    # Return the maximum eigenvalue
    return np.max(eigenvalues)

#=======================================================================
cdef MC_step(double[:,:] arr, double Ts, int nmax):
    """
    Arguments:
	    arr (float(nmax,nmax)) = array that contains lattice data;
	    Ts (float) = reduced temperature (range 0 to 2);
      nmax (int) = side length of square lattice.
    Description:
      Function to perform one MC step, which consists of an average
      of 1 attempted change per lattice site.  Working with reduced
      temperature Ts = kT/epsilon.  Function returns the acceptance
      ratio for information.  This is the fraction of attempted changes
      that are successful.  Generally aim to keep this around 0.5 for
      efficient simulation.
	  Returns:
	    accept/(nmax**2) (float) = acceptance ratio for current MCS.
    """
    cdef:
      double scale = 0.1 + Ts
      int local_accept = 0
      int total_accept
      long[:,:] xran = np.random.randint(0, high=nmax, size=(nmax, nmax))
      long[:,:] yran = np.random.randint(0, high=nmax, size=(nmax, nmax))
      double[:,:] aran = np.random.normal(scale=scale, size=(nmax, nmax))
      long ix, iy
      double ang, en0, en1, boltz, ratio

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Determine the portion of the lattice to work on for each process
    rows_per_process = nmax // size
    start_row = rank * rows_per_process
    end_row = start_row + rows_per_process

    scale = 0.1 + Ts

    xran = np.random.randint(0, high=nmax, size=(nmax, nmax))
    yran = np.random.randint(0, high=nmax, size=(nmax, nmax))
    aran = np.random.normal(scale=scale, size=(nmax, nmax))

    for i in range(start_row, end_row):
        for j in range(nmax):
            ix = xran[i, j]
            iy = yran[i, j]
            ang = aran[i, j]
            en0 = one_energy(arr, ix, iy, nmax)
            arr[ix, iy] += ang
            en1 = one_energy(arr, ix, iy, nmax)
            if en1 <= en0:
                local_accept += 1
            else:
                boltz = exp(-(en1 - en0) / Ts)
                if boltz >= np.random.uniform(0.0, 1.0):
                    local_accept += 1
                else:
                    arr[ix, iy] -= ang

    # Sum the local_accept values from all processes
    total_accept = comm.allreduce(local_accept, op=MPI.SUM)
    ratio = total_accept / (nmax*nmax)
    return ratio

#=======================================================================
def main(program, int nsteps, int nmax, double temp):
    """
    Arguments:
	  program (string) = the name of the program;
	  nsteps (int) = number of Monte Carlo steps (MCS) to perform;
      nmax (int) = side length of square lattice to simulate;
	  temp (float) = reduced temperature (range 0 to 2);
	  pflag (int) = a flag to control plotting.
      threads (int) = number of threads.
    Description:
      This is the main function running the Lebwohl-Lasher simulation.
    Returns:
      NULL
    """
    #===== MPI =====#
    comm = MPI.COMM_WORLD
    taskid = comm.Get_rank()
    numtasks = comm.Get_size()
    numworkers = numtasks-1

    MAXWORKER = nmax # maximum number of worker tasks
    MINWORKER = 1 # minimum number of worker tasks
    MASTER = 0 # taskid of first process

    #===== MASTER =====#
    if taskid == MASTER:
      # Check if numworkers is within range - quit if not
      if (numworkers > MAXWORKER) or (numworkers < MINWORKER):
          print("ERROR: the number of tasks must be between %d and %d." % (MINWORKER+1,MAXWORKER+1))
          print("Quitting...")
          comm.Abort()

    # Create and initialise lattice
    lattice = initdat(nmax)

    # Create arrays to store energy, acceptance ratio and order parameter
    cdef:
      double[:] energy = np.zeros(nsteps+1)
      double[:] ratio = np.zeros(nsteps+1)
      double[:] order = np.zeros(nsteps+1)

    # Set initial values in arrays
    energy[0] = all_energy(lattice,nmax)
    ratio[0] = 0.5 # ideal value
    order[0] = get_order(lattice,nmax)

    # Begin doing and timing some MC steps.
    initial = time.time()
    for it in range(1,nsteps+1):
        ratio[it] = MC_step(lattice,temp,nmax)
        energy[it] = all_energy(lattice,nmax)
        order[it] = get_order(lattice,nmax)
    final = time.time()
    runtime = final-initial
    
    # Final outputs
    if taskid == MASTER:
        print("{}: Size: {:d}, Steps: {:d}, T*: {:5.3f}: Order: {:5.3f}, Time: {:8.6f} s".format(program, nmax,nsteps,temp,order[nsteps-1],runtime))
        # Plot final frame of lattice and generate output file
        savedat(lattice,nsteps,temp,runtime,ratio,energy,order,nmax)

#=======================================================================
# Main part of program, getting command line arguments and calling
if __name__ == '__main__':
    if int(len(sys.argv)) == 4:
        PROGNAME = sys.argv[0]
        ITERATIONS = int(sys.argv[1])
        SIZE = int(sys.argv[2])
        TEMPERATURE = float(sys.argv[3])
        main(PROGNAME, ITERATIONS, SIZE, TEMPERATURE)
    else:
        print("Usage: python {} <ITERATIONS> <SIZE> <TEMPERATURE>".format(sys.argv[0]))
#=======================================================================
