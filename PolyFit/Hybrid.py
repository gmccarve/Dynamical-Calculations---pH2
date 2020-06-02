import pandas as pd
import numpy as np
import math
import time
import sys
from pandas_ods_reader import read_ods
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import KNeighborsRegressor
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


R_vals = [2.5, 2.75, 3.0,  3.25, 3.5,  3.74, 3.79, 3.84, 4.0, 4.25,
          4.5, 4.75, 5.0,  5.25, 5.5,  5.75, 6.0,  6.25, 6.5, 6.75,
          7.0, 7.25, 7.53, 7.58, 7.63, 7.75, 8.0,  8.25, 8.5, 8.75,
          9.0]

KNN_Train_labels = np.loadtxt("H.txt")

FIT = np.zeros((110,16))
for j in range(110):
    FIT[j,:] = np.polyfit(R_vals, KNN_Train_labels[j,:], 15)

np.savetxt("Water_H.polyfit", FIT)

