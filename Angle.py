import numpy as np

def Angle(a, b, c):
    ba = a - b
    bc = c - b

    if np.linalg.norm(ba) == 0 or np.linalg.norm(bc) == 0 :
        return 0

    angle = np.arccos(np.dot(ba,bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))

    return angle
    
