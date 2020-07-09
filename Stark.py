import numpy as np
from Kinetic import W_kin
from HCP import Shell as Lattice
from BraKet import LittleD
from Rotation import Rotate
from Quadrature import QuadValues
from INPUT import *
from PES import *
import sys

def sort_eig(A):
    eigVal, eigVec = np.linalg.eig(A)

    idx = eigVal.argsort()[::-1]
    eigVal = eigVal[idx]
    eigVec = eigVec[:,idx]

    return eigVal, eigVec


def main():

    # Calculate Kinetic energy of water molecule in the 
    # absence of an external field
    # First Parameter is the Jmax value
    # Second parameter is the vibrational excited state

    # Water excited states:
        # 0 - 0,0,0
        # 1 - 1,0,0
        # 2 - 0,1,0
        # 3 - 0,0,1

    KinE = W_kin(Jmax, State)

    # Initialize position of molecule in an experimentally
    # determined geometry with the center of mass set at (0,0,0)

    H2O = np.array([[-0.123085931651,  0.000000000000, 0.0000000000],
                    [ 0.984686107042,  1.430967820506, 0.0000000000],
                    [ 0.984686107042, -1.430967820506, 0.0000000000]])


    # Initialize quadrature values for alpha, beta, and gamma

    A_max, B_max, G_max = MaxQuad, MaxQuad, MaxQuad

    AA, BB, BW, GG = QuadValues(A_max, B_max, G_max)

    # Begin looping through the J, K, and M values

    col = 0
    row = 0

    dim = int((4.0/3.0) * Jmax**3 + 4*Jmax**2 + (11.0/3.0) * Jmax + 1.0000001)

    Hamil = np.zeros((dim, dim))

    for J in range(0, Jmax+1):
        BraNorm = np.sqrt((2 * J + 1) / (8 * np.pi**2))
        for K in range(-J, J+1):
            for M in range(-J, J+1):

                Bra_Nmin = max(0, K-M)
                Bra_Nmax = min(J-K, J+K)

                for JJ in range(0, Jmax+1):
                    KetNorm = np.sqrt((2 * JJ + 1) / (8 * np.pi**2))
                    for KK in range(-JJ, JJ+1):
                        for MM in range(-JJ, JJ+1):

                            Ket_Nmin = max(0, KK-MM)
                            Ket_Nmax = min(JJ-KK, JJ+KK)

                            H = 0.

                            angle_count = 0

                            for a in range(0, A_max):
                                alpha = (2 * np.pi * a) / A_max 
                                for b in range(0, B_max):
                                    beta = BB[b]
                                    beta_weight = BW[b]

                                    LilD_bra = LittleD(J, K, M, beta)
                                    LilD_ket = LittleD(JJ, KK, MM, beta)

                                    for g in range(0, G_max):
                                        gamma = (2 * np.pi * g) / G_max

                                        E = 0.

                                        Rot = Rotate(alpha, beta, gamma)

                                        BigD_bra = np.exp(-(M*alpha + K*gamma) * complex(0,1)) * LilD_bra
                                        BigD_ket = np.exp(-(MM*alpha + KK*gamma) * complex(0,1)) * LilD_ket

                                        BRA = BigD_bra.conjugate() * BraNorm
                                        KET = BigD_ket * KetNorm

                                        weight = ((2. * np.pi) / A_max) * ((2. * np.pi) / G_max) * beta_weight

                                        E += weight * BRA * KET * float(sys.argv[1])/100 * dipole * np.sin(beta) * np.cos(gamma)

                                        H += E

                            if abs(H.real) < thresh:
                                H = 0.0

                            Hamil[col, row] = H.real

                            col += 1

                col = 0
                row += 1

    val, vec = sort_eig(Hamil + KinE)

    #print (val.real)

    for j in val:
        print ('%.4f' % j.real)

    return

if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True, linewidth=20)
    main()

