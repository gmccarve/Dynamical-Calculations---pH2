import numpy as np
from sklearn.neighbors import KNeighborsRegressor as KNR
import sys
import time
import matplotlib.pyplot as plt

def IE(pos, XYZ, FIT, absc, weight, quad, sig):

    if np.linalg.norm(pos) > 3.5:
        idx = np.argpartition(np.linalg.norm(XYZ-pos, axis=1), 2)[:2]
        ReCalc = False
    else:
        ReCalc = True

    E = 0.

    for jx in range(quad):
        for jy in range(quad):
            for jz in range(quad):

                weight_tot = (weight[jx] * weight[jy] * weight[jz])

                move = sig * np.array(([absc[jx], absc[jy], absc[jz]]))

                new_pos = pos + move

                R = np.linalg.norm(new_pos)
                new_pos = new_pos / R

                if ReCalc == True:
                    idx = np.argpartition(np.linalg.norm(XYZ-new_pos, axis=1), 2)[:2]

                E += ((np.poly1d(FIT[idx[0],:])(R) + np.poly1d(FIT[idx[1],:])(R)) / 2.) * weight_tot

    return E / sum(weight)**3

def H2H():
    return 

def H2H2():
    return


if __name__ == '__main__':

    start = time.time()

    Sig = 0.72

    E = np.zeros((1000))

    SmearQuad = 4

    XYZ = np.loadtxt("XYZ")
    FIT = np.loadtxt("Water_H.polyfit")

    Smear_Absc, Smear_weight = np.polynomial.legendre.leggauss(SmearQuad)

    j = np.linspace(2.5, 8.0, 1000)
        
    count = 0
    for jj in j:
        #print (jj)
        E[count] = IE(np.array(([jj, 0., 0.])), XYZ, FIT, Smear_Absc, Smear_weight, SmearQuad, Sig)

        count += 1
    
    plt.scatter(j, E)

    for kk in E:
        print (kk)
    #plt.legend()
    #plt.show()

    #print (time.time() - start)
