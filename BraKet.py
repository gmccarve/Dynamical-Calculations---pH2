import numpy as np
import math

def fac(x):
    return math.factorial(x)

def cos(x):
    return math.cos(x)

def sin(x):
    return math.sin(x)

def sqrt(x):
    return math.sqrt(x)

def exp(x):
    try:
        return math.exp(x)

    except TypeError:
        return cmath.exp(x)

def acos(x):
    return math.acos(x)

def LittleD(J, K, M, beta):
    
    lilD = 0

    nmin = max(0, (K-M))
    nmax = min((J-M), (J+K))

    for N in range(nmin, nmax+1):

        lilw = sqrt(fac(J+M) * fac(J-M) * fac(J+K) * fac(J-K)) / (fac(J-M-N) * fac(J+K-N) * fac(N+M-K) * fac(N))
        BigW = lilw * ((cos(beta/2.0))**(2*J+K-M-2*N)) * ((-sin(beta/2.0))**(M-K+2*N))

        lilD = lilD + ((-1)**N) * BigW

    return lilD
