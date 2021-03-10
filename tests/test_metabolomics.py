"""Tests for metabolomics.py using pytest."""

import pytest

from minedatabase.metabolomics import MetabolomicsDataset, Peak

@pytest.fixture
def known_peaks():
    """Return a list of Peak instances with data (identified)."""
    peak1 = Peak(name='Test1Known', r_time=5.00, mz=867.1391, charge=True,
                 inchi_key='IRPOHFRNKHKIQA-UHFFFAOYSA-N')
    peak2 = Peak(name='Test2Known', r_time=8.00, mz=260.0297, charge=False,
                 inchi_key='HXXFSFRBOHSIMQ-FPRJBGLDSA-N')
    return [peak1, peak2]


@pytest.fixture
def unknown_peaks():
    """Return a list of Peak instances with data (unidentified)."""
    peak1 = Peak(name='Test1Unknown', r_time=3.00, mz=427.0294, charge=True)
    peak2 = Peak(name='Test2Unknown', r_time=1.50, mz=180.0000, charge=False)
    return [peak1, peak2]


@pytest.fixture
def native_cpd_set():
    """Return a set of native compound IDs."""
    return set('C03e0b10e6490ce79a7b88cb0c4e17c2bf6204352',
               'C189c6f78587772ccbd0c2b9e118cf88ff35c316a',
               'C90ba1f8d4d84f6305539f6d05a74497d4d5dfe06')


@pytest.fixture
def metabolomics_dataset(known_peaks, unknown_peaks, native_cpd_set):
    """Default instance of MetabolomicsDataset."""
    met_dataset = MetabolomicsDataset(name='Test Dataset',
                                      adducts=['[M-H]-', '[M+H]+'],
                                      known_peaks=known_peaks,
                                      unk_peaks=unknown_peaks,
                                      native_set=native_cpd_set)
                                      


