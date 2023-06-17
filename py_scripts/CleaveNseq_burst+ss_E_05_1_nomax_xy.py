import numpy as np
import scipy.optimize
import sys
import matplotlib.pyplot as plt
def burst_ss(x, k2, k3, E):
    return E * ((np.square(k2/(k2+k3))*(1-np.exp(-1*(k2+k3)*x))) + ((k2 * k3 * x)/(k2 + k3)))
#load data
names = []
y_list = []
with open(sys.argv[1]) as infile:
    for line in infile:
        fields = line.split("\t")
        names.append(str (fields[0]))
        y_list.append([float(i) for i in fields[1:len(fields)]])
p0 = (0.01, .0001, 1) # start with values near those we expect
x = np.asarray (y_list[0])
for i in range(1,len(y_list)):
    try:
        params, cv = scipy.optimize.curve_fit(burst_ss, x, y_list[i], p0, bounds=( [0, 0, 0.5], [1000, 0.0001, 1]), maxfev=10000)
        k2, k3, E = params
        squaredDiffs = np.square(y_list[i] - burst_ss(x, k2,k3, E))
        squaredDiffsFromMean = np.square(y_list[i] - np.mean(y_list[i]))
        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
    except RuntimeError:
        k2, k3, E, rSquared = '0', '0', '0', '0'
    print (names[i],k2,k3,E,rSquared,sep='\t', file = sys.stdout)
