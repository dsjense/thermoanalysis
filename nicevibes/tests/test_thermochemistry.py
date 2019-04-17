#!/usr/bin/env python3

import os
from pathlib import Path

from pytest import approx

from nicevibes.constants import AU2J_MOL, KCAL2J, CAL2J
from nicevibes.main import thermochemistry
from nicevibes.QCData import QCData


THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))


def test_g16_thermochemistry():
    log = str(THIS_DIR / "logs/04_dmso_hf_freq.log")
    qc = QCData(log, point_group="c1")
    T = 298.15
    thermo = thermochemistry(qc, T, kind="rrho")

    zpe_ref = 0.083950 * AU2J_MOL
    assert thermo.ZPE == approx(zpe_ref)

    u_trans_ref = 0.889 * KCAL2J
    assert thermo.U_trans  == approx(u_trans_ref, rel=1e-3)
    u_rot_ref = 0.889 * KCAL2J
    assert thermo.U_rot  == approx(u_trans_ref, rel=1e-3)
    u_vib_ref = 54.676 * KCAL2J
    assert thermo.U_vib  == approx(u_vib_ref, rel=1e-1)

    s_el_ref = 0.
    assert thermo.S_el == approx(s_el_ref)
    s_trans_ref = 38.978 * CAL2J
    assert (thermo.S_trans) == approx(s_trans_ref)
    s_rot_ref = 25.168 * CAL2J
    assert (thermo.S_rot) == approx(s_rot_ref, rel=1e-3)
    s_vib_ref = 11.519 * CAL2J
    assert (thermo.S_vib) == approx(s_vib_ref, rel=1e-3)


def test_orca_thermochemistry():
    log = str(THIS_DIR / "logs/05_dmso_hf_orca_freq.out")
    qc = QCData(log, point_group="c1")
    T = 298.15
    thermo = thermochemistry(qc, T, kind="qrrho")

    zpe_ref = 0.08393782 * AU2J_MOL
    assert thermo.ZPE == approx(zpe_ref, rel=1e-3)

    u_trans_ref = 0.00141627 * AU2J_MOL
    assert thermo.U_trans  == approx(u_trans_ref, rel=1e-3)
    u_rot_ref = 0.00141627 * AU2J_MOL
    assert thermo.U_rot  == approx(u_trans_ref, rel=1e-3)

    s_el_ref = 0.
    assert thermo.S_el == approx(s_el_ref)
    s_trans_ref = 0.01852169 * AU2J_MOL
    assert (thermo.S_trans*T) == approx(s_trans_ref, rel=1e-1)
    s_rot_ref = 0.01092172 * AU2J_MOL
    assert (thermo.S_rot*T) == approx(s_rot_ref, rel=1e-1)
    s_vib_ref = 0.00546508 * AU2J_MOL
    assert (thermo.S_vib*T) == approx(s_vib_ref, rel=1e-4)
