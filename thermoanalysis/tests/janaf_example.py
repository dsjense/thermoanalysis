#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

from thermoanalysis.constants import AU2J, NA, C, J2AU, J2CAL
from thermoanalysis.thermo import thermochemistry
from thermoanalysis.QCData import QCData
from scipy.constants.constants import calorie, N_A


def print_thermo_comparison(tarr, arr1, arr2, header1, header2, title):
    len_title = len(title)
    print('-' * len_title)
    print(title)
    print('-' * len_title)
    diffs = [np.concatenate([[0], arr[1:] - arr[:-1]]) for arr in (arr1, arr2)]
    cols = np.array([tarr, arr1, arr2, diffs[0], diffs[1]]).T
    headers = ['T (K)', header1, header2, header1 + ' Diffs', header2+' Diffs']
    print(tabulate(cols, headers))
    print('')


def parse_janaf_str(janaf_str):
    janaf_str = janaf_str.replace('INFINITE', '10000000')
    # The first line is blank so we skip it with `[1:]`
    data_arr = np.array([[float(val) for val in line.split()]
                         for line in janaf_str.splitlines()[1:]])
    headers = ['T', 'C', 'S', 'gef', 'dH', 'dHf', 'dGf', 'logK']
    data = {}
    for i, header in enumerate(headers):
        data[header] = data_arr[:, i]
    return data


def find_index(arr, val, tol=1e-16):
    index = np.argwhere(np.abs(arr - val) < tol)
    if len(index) != 1:
        raise Exception("Index corresponding to val={} in arr={} "
                        "was not found.".format(val, arr))
    return index[0][0]
    

def h2o_pbe0_orca(tol=1e-16):
    janaf_str = """
0        0.          0.       INFINITE     -9.904     -238.921     -238.921     INFINITE
100     33.299     152.388     218.534     -6.615     -240.083     -236.584     123.579
200     33.349     175.485     191.896     -3.282     -240.900     -232.766     60.792
298.15  33.590     188.834     188.834      0.        -241.826     -228.582     40.047
300     33.596     189.042     188.835      0.062     -241.844     -228.500     39.785
400     34.262     198.788     190.159      3.452     -242.846     -223.901     29.238
500     35.226     206.534     192.685      6.925     -243.826     -219.051     22.884
600     36.325     213.052     195.550     10.501     -244.758     -214.007     18.631
700     37.495     218.739     198.465     14.192     -245.632     -208.812     15.582
800     38.721     223.825     201.322     18.002     -246.443     -203.496     13.287
900     39.987     228.459     204.084     21.938     -247.185     -198.083     11.496
1000    41.268     232.738     206.738     26.000     -247.857     -192.590     10.060
"""
    janaf_data = parse_janaf_str(janaf_str)
    log = "logs/07_h2o_pbe0_orca_freq.out"
    qc = QCData(log, point_group="c2v")
    temperatures = janaf_data['T']
    temperatures_nonzero = [tol] + [val for val in temperatures[1:]]
    tref = 298.15  # Reference temperature
    tref_index = find_index(temperatures, tref, tol)
    thermos = [thermochemistry(qc, T) for T in temperatures_nonzero]
    header1 = 'ORCA'
    header2 = 'JANAF'
    
    CPs_janaf = janaf_data['C']
    CPs = np.array([thermo.C_P for thermo in thermos])
    title = 'Constant Pressure Heat Capacities in J / (mol * K)'
    print_thermo_comparison(temperatures, CPs, CPs_janaf, header1, header2, title)
    
    Ss_janaf = janaf_data['S']
    Ss = np.array([AU2J * NA * thermo.S_tot for thermo in thermos])
    title = 'Entropy in J / (mol * K)'
    print_thermo_comparison(temperatures, Ss, Ss_janaf, header1, header2, title)
    
    dH_janaf = janaf_data['dH']
    H = np.array([AU2J * NA * thermo.H / 1000e0 for thermo in thermos])
    dH = H - H[tref_index]
    title = 'H(T) - H({}) in kJ / mol'.format(tref)
    print_thermo_comparison(temperatures, dH, dH_janaf, header1, header2, title)

    gef_janaf = janaf_data['gef']
    gef = Ss - (dH * 1000e0) / temperatures_nonzero
    title = '-[G(T) - H({})]/T in kJ / mol'.format(tref)
    print_thermo_comparison(temperatures, gef, gef_janaf, header1, header2, title)

    # All energies are converted to kJ
    AU2kJ = AU2J / 1000e0
    H2O_E0 = (thermos[0].U_el + thermos[0].ZPE) * AU2kJ * N_A
    print('thermos[0].U_el={}, thermos[0].ZPE={}'.format(thermos[0].U_el, thermos[0].ZPE))
    print(thermos[0].U_el)
    # H2O_E0 = -76.487766198369 * AU2kJ * N_A
    H2O_E0 = (-76.48776620 + 0.02087626) * AU2kJ * N_A
    # zpe = thermos[0].ZPE * AU2kJ * N_A
    # print('zpe={} AU/particle'.format(thermos[0].ZPE))
    # Electronic energy                ...    -76.48776620 Eh
    # Zero point energy                ...      0.02087626 Eh      13.10 kcal/mol
    # H2O_E0 += zpe
    H_E0 = -0.499631513332 * AU2kJ * N_A
    # H_E0 = -0.5 * AU2kJ * N_A
    print('AU2kJ * N_A={}'.format(AU2kJ * N_A))
    # H_E0 = (-0.501051104433 - 0.00000000) * AU2kJ * N_A
    O_E0 = (-75.12949136 - 0.00000003) * AU2kJ * N_A
    # O_E0 = -75.129491363507 * AU2kJ * N_A
    print('H2O_E0={}, H_E0={}, O_E0={}'.format(H2O_E0, H_E0, O_E0))
    dHf_H = 51.63 * calorie
    dHf_O = 58.99 * calorie
    dHf_H2O_0K = (2 * dHf_H + 1 * dHf_O) - (2 * H_E0 + 1 * O_E0 - H2O_E0)
    print('2 * dHf_H + 1 * dHf_O={}'.format(2 * dHf_H + 1 * dHf_O))
    print('dHf_H2O_0K=', dHf_H2O_0K)
    dH_H2O = 9.904
    dH_H = 8.467 / 2e0
    dH_O = 8.683 / 2e0
    # dHf_H2O_298K = dHf_H2O_0K + dH_H2O - thermos[0].ZPE * AU2kJ * N_A - (2 * dH_H + 1 * dH_O)
    dHf_H2O_298K = dHf_H2O_0K + dH_H2O - (2 * dH_H + 1 * dH_O)
    print('dH_H2O - (2 * dH_H + 1 * dH_O)={}'.format(dH_H2O - (2 * dH_H + 1 * dH_O)))
    print('dHf_H2O_298K={}'.format(dHf_H2O_298K))
    # -241.826
    
    # Example using the H2 + .5 O2 -> H2O reaction
    # H2 and O2 are both reference states with enthalpies of formation equal to zero
    dHf_H2 = 0e0
    dHf_O2 = 0e0
    H2_E0 = (-1.16842993 + 0.01006000) * AU2kJ * N_A
    O2_E0 = (-150.39167425 + 0.00386815) * AU2kJ * N_A
    H2O_E0 = (-76.48801045 + 0.02147735) * AU2kJ * N_A
    
    # Gaussian MO6
    H2_E0 = (-1.17221845130 + 0.009902) * AU2kJ * N_A
    O2_E0 = (-150.322374247+0.003897) * AU2kJ * N_A
    H2O_E0 = (-76.4321360023 + 0.021635) * AU2kJ * N_A
    
    # NWChem with aug-pc-2
    # H2_E0 = (  -1.180400347676 + 0.010054) * AU2kJ * N_A
    # O2_E0 = (-150.392526993215 + 0.003717) * AU2kJ * N_A
    # H2O_E0 = (-76.470735757716 + 0.021266) * AU2kJ * N_A
    
    # NWChem with aug-pc-3
    # H2_E0 = (-1.18070033 + 0.010056) * AU2kJ * N_A
    # O2_E0 = (-150.39936484 + 0.003730) * AU2kJ * N_A
    # H2O_E0 = (-76.47419695 + 0.021292) * AU2kJ * N_A
    
    # ORCA dispersion-corrected energies
    # H2_E0 = (  -1.16842993 + 0.01006000) * AU2kJ * N_A
    # O2_E0 = (-150.45700775 + 0.00389843) * AU2kJ * N_A
    # H2O_E0 = (-76.48776620 + 0.02087626) * AU2kJ * N_A
    vals = [0.01006000, 0.00389843, 0.02087626]
    vals = [0.0096868, 0.003308, 0.020845]
    vals = [0.010054, 0.003717, 0.021266]
    for val in vals: print('{:.6f}'.format(val * AU2kJ * N_A))
    # CFOUR CCSD(T) PVQZ
    # H2_E0 =   -1.173840027264172 * AU2kJ * N_A + 26.3347
    # O2_E0 = -150.204237625165092 * AU2kJ * N_A + 9.4961
    # H2O_E0 = -76.401162130790752 * AU2kJ * N_A + 56.9863
    
    dHf_H2O_0K = (dHf_H2 + 0.5 * dHf_O2) - (H2_E0 + 0.5 * O2_E0 - H2O_E0)
    print('Hello: dHf_H2O_0K={:.6g}'.format(dHf_H2O_0K))
    dH_H2 = 0.013206 * AU2kJ * N_A
    # dH_H2 = 8.467
    dH_O2 = 0.007203 * AU2kJ * N_A
    # dH_O2 = d8.683
    dH_H2O = 0.025415 * AU2kJ * N_A
    # dH_H2O = 9.904
    print('dH_H2O={}, dH_H2={}, dH_O2={}'.format(dH_H2O, dH_H2, dH_O2))
    dHf_H2O_298K = dHf_H2O_0K + dH_H2O - (dH_H2 + 0.5 * dH_O2)
    print('dHf_H2O_298K={}'.format(dHf_H2O_298K))
    print('AU2kJ*N_A={}'.format(AU2kJ*N_A/calorie))
    from scipy.constants import Boltzmann
    print('k_B*298.15 K = {} kJ/mol'.format(Boltzmann*298.15*N_A/1e3))

if __name__ == "__main__":
    h2o_pbe0_orca()
