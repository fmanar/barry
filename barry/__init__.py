import numpy as np

from . import mvpoly

def fill(lam):
    """complete the coordinates

    based on sigma(lambda_i) = 1

    will someday combine all the arguments into a list or a tuple or something

    """
    total = 0
    ind = None
    for i, l in enumerate(lam):
        if l is None:
            ind = i
        else:
            total += l
    lam[ind] = 1 - total
    return lam

def get_pos(basis, lam):
    return np.matmul(basis, lam)

def get_lam(basis, pos):
    d = basis.shape[0]
    A = np.ones((d+1, d+1))
    A[0:d,:] = basis
    b = np.ones(((d+1), 1))
    b[0:d,:] = pos
    return np.linalg.solve(A, b)


    