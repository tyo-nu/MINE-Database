"""Tests for thermodynamics.py using pytest."""


import json
from pathlib import Path
import pytest
from minedatabase.thermodynamics import Thermodynamics

file_path = Path(__file__)
file_dir = file_path.parent

# Thermodynamics
loaded_db = False
thermo = Thermodynamics()

try:
    thermo.load_thermo_from_postgres()
    loaded_db = True
except:
    pass

if not loaded_db:
    try:
        thermo.load_thermo_from_sqlite()
        loaded_db = True
    except:
        pass

@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_get_compound(pk_transformed):
    from_pk = thermo.get_eQ_compound_from_cid(
        'C75ec1165c59bd4ec72b0636a8ae4904a20b4c4f3', pk_transformed)
    cpd_smiles = pk_transformed.compounds['C75ec1165c59bd4ec72b0636a8ae4904a20b4c4f3']["SMILES"]
    from_eq = thermo.pc.get_compounds(cpd_smiles)

    assert from_pk == from_eq

@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_dgf_standard(pk_transformed):
    dgf = thermo.standard_dg_formation_from_cid(
        'C75ec1165c59bd4ec72b0636a8ae4904a20b4c4f3', pk_transformed)

    assert abs(dgf[0]) == pytest.approx(2580, 2582) 

@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_dgr_standard(pk_transformed):
    dgr = thermo.physiological_dg_prime_from_rid(
        'R6036f02cd619e4515f5180a4d650c27b7a8b40b59dc130a6c8de526d5adaeb0e',
        pk_transformed
    )

    assert abs(dgr.value.magnitude) == pytest.approx(13, 15) 
