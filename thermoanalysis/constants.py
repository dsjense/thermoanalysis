import scipy.constants as spc
from scipy.constants.codata import (_physical_constants_2002, 
    _physical_constants_2006, _physical_constants_2010)


physical_constants = {}
for val in (_physical_constants_2002, _physical_constants_2006, 
            _physical_constants_2010):
    physical_constants.update(val)

def value(key):
    return physical_constants[key][0]

# We force the use of CODATA 2010 values to agree with cclib.
 
# AU2EV = spc.value("Hartree energy in eV")  # Hartree to eV
AU2EV = value("Hartree energy in eV")  # Hartree to eV
# AU2J = spc.value("Hartree energy")  # Hartree to J
AU2J = value("Hartree energy")  # Hartree to J
AU2KJ = AU2J / 1000e0
CAL2J = spc.calorie
J2CAL = 1 / CAL2J  # Joule to calorie
# J2AU = 1 / spc.value("Hartree energy")  # Joule to Hartree
J2AU = 1 / AU2J  # Joule to Hartree
ANG2M = 1e-10  # Angstrom to meter
# ANG2AU = 1 / spc.value("Bohr radius") * ANG2M  # Angstrom to Bohr
ANG2AU = 1 / value("Bohr radius") * ANG2M  # Angstrom to Bohr
AMU2KG = value("unified atomic mass unit")  # Atomic mass units to kg
# R = spc.R  # J/(K*mol), ideal gas constant
R = value('molar gas constant')  # J/(K*mol), ideal gas constant
# C = spc.c  # Speed of light in m/s
C = value('speed of light in vacuum')
# PLANCK = spc.Planck  # J/s, Planck constant
PLANCK = value('Planck constant')
# KB = spc.Boltzmann  # J/K, Boltzmann constant
KB = value('Boltzmann constant')
KBAU = KB * J2AU
NA = spc.N_A  # 1/mol, Avogadro constant
AU2J_MOL = 1 / J2AU * NA
KCAL2J = CAL2J * 1000
CAL_MOL2AU = CAL2J * J2AU / NA
KCAL_MOL2AU = 1000 * CAL_MOL2AU
