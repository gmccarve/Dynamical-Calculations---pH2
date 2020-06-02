import numpy as np
from sklearn.neighbors import KNeighborsRegressor as KNR


def H2O_H(a):

    x, y, z = a[0], a[1], a[2]

    XYZ = np.loadtxt("XYZ")

    FIT = np.loadtxt("Water_H.polyfit")

    R = np.sqrt(x*x + y*y + z*z)

    neigh = KNR(n_neighbors=2, p=2)
    neigh.fit(XYZ, np.zeros((110)))

    dist, ID = neigh.kneighbors([[x, y, z]])

    E = (np.poly1d(FIT[ID[0][0],:])(R) + np.poly1d(FIT[ID[0][1],:])(R)) / 2.

    return E

def H2O_H2(a):

    x, y, z = a[0], a[1], a[2]

    XYZ = np.loadtxt("XYZ")

    FIT = np.loadtxt("Water_H2.polyfit")
    
    R = np.sqrt(x*x + y*y + z*z)
    
    neigh = KNR(n_neighbors=2, p=2)
    neigh.fit(XYZ, np.zeros((110)))
    
    dist, ID = neigh.kneighbors([[x, y, z]])

    E = (np.poly1d(FIT[ID[0][0],:])(R) + np.poly1d(FIT[ID[0][1],:])(R)) / 2.

    return E


def H2H():
    return 

def H2H2():
    return


