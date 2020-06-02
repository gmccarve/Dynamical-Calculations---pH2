import numpy as np

def Rotate(A, B, G):

    cosA = np.cos(A)
    cosB = np.cos(B)
    cosG = np.cos(G)

    sinA = np.sin(A)
    sinB = np.sin(B)
    sinG = np.sin(G)

    Rot = np.array([[ cosA * cosB * cosG - sinA * sinG, -cosA * cosB * sinG - sinA * cosG, cosA * sinB],
                    [ sinA * cosB * cosG + cosA * sinG, -sinA * cosB * sinG + cosA * cosG, sinA * sinB],
                    [-sinB * cosG,                       sinB * sinG,                      cosB       ]])


    return Rot

