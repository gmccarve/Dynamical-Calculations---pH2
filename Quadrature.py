import numpy as np
import sys

def QuadValues(A, B, G):

    AA = np.zeros((A))
    BB = np.zeros((B))
    GG = np.zeros((G))

    for i in range(1, A+1):
        AA[i-1] = (2 * np.pi * (i-1)) / A

    for i in range(1, G+1):
        GG[i-1] = (2 * np.pi * (i-1)) / G

    BetaV, BetaW = np.polynomial.legendre.leggauss(B)

    BB = np.arccos(BetaV)

    return AA, BB, BetaW, GG

