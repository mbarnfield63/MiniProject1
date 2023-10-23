import pytest
import numpy as np
from ..code.LebwohlLasher_numba import *

def test_initdat():
    # testing values
    nmax = 5
    arr = initdat(nmax)
    
    # test shape and output is array
    assert arr.shape == (nmax, nmax)
    # test values between 0 and 2pi
    assert (0 <= arr).all() and (arr <= 2*np.pi).all()

def test_one_energy():
    # testing values
    nmax = 5
    np.random.seed(42)
    arr = np.random.random_sample((nmax, nmax))*2.0*np.pi
    ix, iy = 2, 2
    
    # test function
    energy = one_energy(arr, ix, iy, nmax)
    assert energy == -1.0512936529391477

def test_all_energy():
    # testing values
    nmax = 5
    np.random.seed(42)
    arr = np.random.random_sample((nmax, nmax))*2.0*np.pi
    
    # test function
    all = all_energy(arr, nmax)
    assert all == -22.709752842971472

def test_get_order():
    # testing values
    nmax = 5
    np.random.seed(42)
    arr = np.random.random_sample((nmax, nmax))*2.0*np.pi

    # test function
    order = get_order(arr, nmax)
    assert order == 0.3238593394231049

def test_MC_step():
    # testing values
    nmax = 5
    np.random.seed(42)
    arr = np.random.random_sample((nmax, nmax))*2.0*np.pi
    Ts = 0.5

    # test function
    MC = MC_step(arr, Ts, nmax)
    assert MC == 0.52