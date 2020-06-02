import numpy as np

def W_kin(jmax, State):

    if State == 0 :
    # v1, v2, v3 = 0, 0, 0
        A = 27.8572554962
        B = 14.5144875043
        C = 9.27986250769

    elif State == 1:
    # v1, v2, v3 = 1, 0, 0
        A = 27.1388599924
        B = 14.2991000050
        C = 9.10149945389

    elif State == 2:
    # v1, v2, v3 = 0, 1, 0
        A = 31.0847300062
        B = 14.6748099834
        C = 9.13606001363

    elif State == 3:
    # v1, v2, v3 = 0, 0, 1
        A = 26.6303599946
        B = 14.4225400054
        C = 9.14182999037

    dim = int((4.0/3.0) * jmax**3 + 4*jmax**2 + (11.0/3.0) * jmax + 1.0000000001)

    IULC = 0

    RotMat = np.zeros((dim,dim))

    for J in range(0,jmax+1):
        for M in range(-J, J+1):
            for K in range(-J, J+1):
                for KK in range(-J, J+1):

                    val = 0.0

                    if K == KK:
                        val = ((0.50) * (B+C) * (J*(J+1)-K**2) + A*K**2)
                    elif K + 2 == KK:
                        val = ((0.25) * (B-C) * np.sqrt((J*(J+1) - K * (K+1))) * np.sqrt(J*(J+1) - (K+1) * (K+2)))
                    elif K - 2 == KK:
                        val = ((0.25) * (B-C) * np.sqrt((J*(J+1) - K * (K-1))) * np.sqrt(J*(J+1) - (K-1) * (K-2)))

                    try:
                        RotMat[IULC+KK+J, IULC+K+J] = val
                    except IndexError:
                        pass

            IULC = IULC + 2*J+1

    return RotMat
