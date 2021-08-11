"""Tests for thermodynamics.py using pytest."""
from pathlib import Path

import pytest


file_path = Path(__file__)
file_dir = file_path.parent
sqlite_loc = file_dir / "../compounds.sqlite"

# Thermodynamics
loaded_db = False

if sqlite_loc.is_file():
    try:
        from equilibrator_api import Q_

        from minedatabase.thermodynamics import Thermodynamics

        thermo = Thermodynamics()

        try:
            thermo.load_thermo_from_postgres()
            loaded_db = True
            print("Thermo DB loaded from Postgres")
        except:
            try:
                thermo.load_thermo_from_sqlite()
                loaded_db = True
                print("Thermo DB loaded from SQLite")
            except:
                pass
    except:
        pass


@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_get_compound(pk_transformed):
    from_pk = thermo.get_eQ_compound_from_cid(
        "C75ec1165c59bd4ec72b0636a8ae4904a20b4c4f3", pk_transformed
    )
    cpd_smiles = pk_transformed.compounds["C75ec1165c59bd4ec72b0636a8ae4904a20b4c4f3"][
        "SMILES"
    ]
    from_eq = thermo.lc.get_compounds(cpd_smiles)

    assert from_pk == from_eq


@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_dgf_standard(pk_transformed):
    dgf = thermo.standard_dg_formation_from_cid(
        "C75ec1165c59bd4ec72b0636a8ae4904a20b4c4f3", pk_transformed
    )

    assert abs(dgf) == pytest.approx(2580, 2582)


@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_dgr_standard(pk_transformed):
    dgr = thermo.physiological_dg_prime_from_rid(
        "R6036f02cd619e4515f5180a4d650c27b7a8b40b59dc130a6c8de526d5adaeb0e",
        pk_transformed,
    )

    assert abs(dgr.value.magnitude) == pytest.approx(13, 15)


@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_dgr_prime(pk_transformed):
    low_p_h = Q_(5)
    high_p_h = Q_(9)
    low_ph_dgp = thermo.dg_prime_from_rid(
        "R6036f02cd619e4515f5180a4d650c27b7a8b40b59dc130a6c8de526d5adaeb0e",
        pk_transformed,
        p_h=low_p_h,
    )
    high_ph_dgp = thermo.dg_prime_from_rid(
        "R6036f02cd619e4515f5180a4d650c27b7a8b40b59dc130a6c8de526d5adaeb0e",
        pk_transformed,
        p_h=high_p_h,
    )

    assert abs(high_ph_dgp.value.magnitude) == pytest.approx(4, 7)
    assert abs(high_ph_dgp.value.magnitude) == pytest.approx(21, 23)
