"""
Basic Python Lebwohl-Lasher code.  Based on the paper 
P.A. Lebwohl and G. Lasher, Phys. Rev. A, 6, 426-429 (1972).
This version in 2D.

Run at the command line by typing:

python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>

where:
  ITERATIONS = number of Monte Carlo steps, where 1MCS is when each cell
      has attempted a change once on average (i.e. SIZE*SIZE attempts)
  SIZE = side length of square lattice
  TEMPERATURE = reduced temperature in range 0.0 - 2.0.
  PLOTFLAG = 0 for no plot, 1 for energy plot and 2 for angle plot.
  
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
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

cimport numpy as np
from libc.math cimport cos, sin, exp
from libc.stdlib cimport rand, srand, RAND_MAX

#=======================================================================
cdef initdat(int nmax):
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
    cdef np.ndarray[np.float64_t, ndim=2] arr = np.random.random_sample((nmax, nmax)) * 2.0 * np.pi
    return arr
#=======================================================================
def plotdat(arr,pflag,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  pflag (int) = parameter to control plotting;
      nmax (int) = side length of square lattice.
    Description:
      Function to make a pretty plot of the data array.  Makes use of the
      quiver plot style in matplotlib.  Use pflag to control style:
        pflag = 0 for no plot (for scripted operation);
        pflag = 1 for energy plot;
        pflag = 2 for angles plot;
        pflag = 3 for black plot.
	  The angles plot uses a cyclic color map representing the range from
	  0 to pi.  The energy plot is normalised to the energy range of the
	  current frame.
	Returns:
      NULL
    """
    if pflag==0:
        return
    u = np.cos(arr)
    v = np.sin(arr)
    x = np.arange(nmax)
    y = np.arange(nmax)
    cols = np.zeros((nmax,nmax))
    if pflag==1: # colour the arrows according to energy
        mpl.rc('image', cmap='rainbow')
        for i in range(nmax):
            for j in range(nmax):
                cols[i,j] = one_energy(arr,i,j,nmax)
        norm = plt.Normalize(cols.min(), cols.max())
    elif pflag==2: # colour the arrows according to angle
        mpl.rc('image', cmap='hsv')
        cols = arr%np.pi
        norm = plt.Normalize(vmin=0, vmax=np.pi)
    else:
        mpl.rc('image', cmap='gist_gray')
        cols = np.zeros_like(arr)
        norm = plt.Normalize(vmin=0, vmax=1)

    quiveropts = dict(headlength=0,pivot='middle',headwidth=1,scale=1.1*nmax)
    fig, ax = plt.subplots()
    q = ax.quiver(x, y, u, v, cols,norm=norm, **quiveropts)
    ax.set_aspect('equal')
    plt.show()
#=======================================================================
def savedat(double[:, ::1] arr, int nsteps, double Ts, double[:] ratio, double[:] energy, double[:] order, int nmax):
    # Create filename based on the current date and time.
    current_datetime = datetime.datetime.now().strftime("%a-%d-%b-%Y-at-%I-%M-%S%p")
    filename = "./Outputs/LL-Output-{:s}.txt".format(current_datetime).encode('utf-8')
    
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
            FileOut.write("   {:05d}    {:6.4f} {:12.4f}  {:6.4f} ".format(i, ratio[i], energy[i], order[i]).encode('utf-8'))
#=======================================================================
cpdef double one_energy(np.ndarray[double, ndim=2] arr, int ix, int iy, int nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  ix (int) = x lattice coordinate of cell;
	  iy (int) = y lattice coordinate of cell;
      nmax (int) = side length of square lattice.
    Description:
      Function that computes the energy of a single cell of the
      lattice taking into account periodic boundaries.  Working with
      reduced energy (U/epsilon), equivalent to setting epsilon=1 in
      equation (1) in the project notes.
	Returns:
	  en (float) = reduced energy of cell.
    """
    cdef double en = 0.0
    cdef int ixp, ixm, iyp, iym
    cdef np.ndarray[double, ndim=1] neighbors

    ixp = (ix + 1) % nmax
    ixm = (ix - 1) % nmax
    iyp = (iy + 1) % nmax
    iym = (iy - 1) % nmax

    # Extract the values of the neighbors
    neighbors = np.array([arr[ixp, iy], arr[ixm, iy], arr[ix, iyp], arr[ix, iym]], dtype=float)

    # Compute the angular differences using a for loop
    cdef double ang_diff
    cdef double cos_squared = 0.0
    for i in range(4):
        ang_diff = arr[ix, iy] - neighbors[i]
        cos_squared += (cos(ang_diff) ** 2)

    # Calculate energy
    en = 0.5 * (1.0 - 3.0 * cos_squared)

    return en
#=======================================================================
cdef double all_energy(np.ndarray[double, ndim=2] arr, int nmax):
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
    cdef np.ndarray[long, ndim=2] ix, iy
    cdef double enall = 0.0

    ix, iy = np.meshgrid(np.arange(nmax, dtype=int), np.arange(nmax, dtype=int))
    # Calculate the energy of the entire lattice
    cdef int i, j
    for i in range(nmax):
        for j in range(nmax):
            enall += one_energy(arr, ix[i, j], iy[i, j], nmax)

    return enall
#=======================================================================
def get_order(np.ndarray[np.float64_t, ndim=2] arr, int nmax):
    cdef np.ndarray[np.float64_t, ndim=2] Qab = np.zeros((3, 3), dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=2] delta = np.eye(3, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=3] lab = np.empty((3, nmax, nmax), dtype=np.float64)

    # Calculate Qab using Cythonized loops
    for a in range(3):
        for b in range(3):
            for i in range(nmax):
                for j in range(nmax):
                    lab[a, i, j] = cos(arr[i, j]) if a == 0 else sin(arr[i, j]) if a == 1 else 0.0
                    Qab[a, b] += 3 * lab[a, i, j] * lab[b, i, j] - delta[a, b]

    Qab /= 2 * nmax * nmax

    cdef np.ndarray[np.float64_t] eigenvalues = np.empty(3, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=2] eigenvectors = np.empty((3, 3), dtype=np.float64)

    # Compute the eigenvalues and eigenvectors using Cythonized np.linalg.eig
    eigenvalues, eigenvectors = np.linalg.eig(Qab)
    
    # Return the maximum eigenvalue
    return eigenvalues.max()

#=======================================================================
def MC_step(np.ndarray[np.float64_t, ndim=2] arr, double Ts, int nmax):
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
    cdef double scale = 0.1 + Ts
    cdef int accept = 0
    cdef np.ndarray[long, ndim=2] xran = np.random.randint(0, high=nmax, size=(nmax, nmax))
    cdef np.ndarray[long, ndim=2] yran = np.random.randint(0, high=nmax, size=(nmax, nmax))
    cdef np.ndarray[np.float64_t, ndim=2] aran = np.random.normal(scale=scale, size=(nmax, nmax))
    cdef int i, j, ix, iy
    cdef double ang, en0, en1, boltz
    cdef double random = RAND_MAX
    cdef double ratio

    for i in range(nmax):
        for j in range(nmax):
            ix = xran[i, j]
            iy = yran[i, j]
            ang = aran[i, j]
            en0 = one_energy(arr, ix, iy, nmax)
            arr[ix, iy] += ang
            en1 = one_energy(arr, ix, iy, nmax)

            if en1 <= en0:
                accept += 1
            else:
                # Now apply the Monte Carlo test - compare
                # exp( -(E_new - E_old) / T* ) >= rand(0,1)
                boltz = exp(-(en1 - en0) / Ts)

                if boltz >= rand() / random:
                    accept += 1
                else:
                    arr[ix, iy] -= ang
    ratio = accept / (nmax * nmax) 
    return ratio

#=======================================================================
def main(program, int nsteps, int nmax, double temp, int pflag):
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
    # Create and initialise lattice
    lattice = initdat(nmax)
    # Plot initial frame of lattice
    plotdat(lattice,pflag,nmax)
    # Create arrays to store energy, acceptance ratio and order parameter
    energy = np.zeros(nsteps+1,dtype=np.dtype)
    ratio = np.zeros(nsteps+1,dtype=np.dtype)
    order = np.zeros(nsteps+1,dtype=np.dtype)
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
    print("{}: Size: {:d}, Steps: {:d}, T*: {:5.3f}: Order: {:5.3f}, Time: {:8.6f} s".format(program, nmax,nsteps,temp,order[nsteps-1],runtime))
    # Plot final frame of lattice and generate output file
    #savedat(lattice,nsteps,temp,ratio,energy,order,nmax)
    plotdat(lattice,pflag,nmax)
#=======================================================================
# Main part of program, getting command line arguments and calling
# main simulation function.
#
if __name__ == '__main__':
    if int(len(sys.argv)) == 5:
        PROGNAME = sys.argv[0]
        ITERATIONS = int(sys.argv[1])
        SIZE = int(sys.argv[2])
        TEMPERATURE = float(sys.argv[3])
        PLOTFLAG = int(sys.argv[4])
        main(PROGNAME, ITERATIONS, SIZE, TEMPERATURE, PLOTFLAG)
    else:
        print("Usage: python {} <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>".format(sys.argv[0]))
#=======================================================================
