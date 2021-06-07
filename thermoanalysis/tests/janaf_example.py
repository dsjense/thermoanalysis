#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

from thermoanalysis.constants import AU2J, NA, C, J2AU, J2CAL
from thermoanalysis.thermo import thermochemistry
from thermoanalysis.QCData import QCData


def diff_arr(arr):
    return arr[1:] - arr[:-1]


def h2o_pbe0_orca():
    log = "logs/07_h2o_pbe0_orca_freq.out"
    qc = QCData(log, point_group="c2v")
    temperatures = [float(val) for val in (100, 200, 298.15, 300, 400, 500, 
                                           600, 700, 800, 900, 1000)]
    CPs_janaf = np.array([33.299, 33.349, 33.590, 33.596, 34.262, 35.226, 36.325, 
                          37.495, 38.721, 39.987, 41.268])
    thermos = [thermochemistry(qc, T) for T in temperatures]
    CPs = np.array([thermo.C_P for thermo in thermos])
    cols = np.array([CPs, CPs_janaf]).T
    print('Constant Pressure Heat Capacities in J / (mol * K)\n')
    print(tabulate(cols, ["ORCA", "JANAF"]))
    print("\nDiffs")
    print(tabulate(cols[1:] - cols[:-1]))
    
    Ss_janaf = np.array([152.388, 175.485, 188.834, 189.042, 198.788, 206.534, 
                         213.052, 218.739, 223.825, 228.459, 232.738])
    Ss = np.array([AU2J * NA * thermo.S_tot for thermo in thermos])
    cols = np.array([Ss, Ss_janaf]).T
    print(tabulate(cols, ["ORCA", "JANAF"]))


if __name__ == "__main__":
    h2o_pbe0_orca()
