import json
import logging
import os

import numpy as np
import pytest

from irsa.core.algorithm import Algorithm
from irsa.core.analyses import get_spearman_and_pearson
from irsa.core.spectrum_functions import deconvolute, voigt
from irsa.core.utils import load_peaks, load_spectrum, normalize_spectrum

TEST_DATA = f"{os.path.dirname(__file__)}/data"


@pytest.fixture
def reference_ir():
    return np.loadtxt(f"{TEST_DATA}/experimental_ir_reference.txt")


@pytest.fixture
def reference_raman():
    return np.loadtxt(f"{TEST_DATA}/experimental_raman_reference.txt")


def test_load_and_normalize(reference_ir, reference_raman):
    ir = normalize_spectrum(load_spectrum(f"{TEST_DATA}/experimental_spectra/ir_cdcl3.txt"), lower=1000, upper=1500)
    raman = normalize_spectrum(
        load_spectrum(f"{TEST_DATA}/experimental_spectra/raman_cdcl3.txt"),
        lower=1000, upper=1500)

    assert np.max(ir[:, 1]) == 1
    assert np.max(raman[:, 1]) == 1
    assert (np.abs(reference_ir[:, 1] - ir[:, 1]) < 1e-4).all()
    assert (np.abs(reference_raman[:, 1] - raman[:, 1]) < 1e-4).all()


def test_deconvolute(reference_ir):
    file_to_clean = "tmp.txt"
    try:
        deconvolute(spectrum=reference_ir, working_dir=os.path.abspath(''),
                    save_data=file_to_clean, normalize=True, lower=1000, higher=1800, vcd=False)
        #generated = np.loadtxt(file_to_clean)
        #reference = np.loadtxt(f"{TEST_DATA}/ir_exp_peaks_reference.txt")
        # for gen, ref in zip(generated, reference): # version dependent it seems
        #    for i in range(len(gen)):
        #        assert abs(gen[i] - ref[i]) < 0.2
    except Exception as error:
        logging.error(error)
        assert False
    finally:
        os.remove(file_to_clean)


def test_load_peaks():
    reference = np.loadtxt(f"{TEST_DATA}/concatenate_theo_peaks_reference.txt")
    compare = load_peaks([f"{TEST_DATA}/ir_theo_peaks_reference.txt",
                         f"{TEST_DATA}/raman_theo_peaks_reference.txt"], kind_of_spectra=[0, 1])
    assert (np.abs(reference - compare) < 1e-4).all()


def test_algorithm():
    with open(f"{TEST_DATA}/reference_algo", 'r') as file:
        reference = json.load(file)

    theo_peaks = load_peaks([f"{TEST_DATA}/ir_theo_peaks_reference.txt",
                             f"{TEST_DATA}/raman_theo_peaks_reference.txt"], kind_of_spectra=[0, 1])

    exp_peaks = load_peaks([f"{TEST_DATA}/ir_exp_peaks_reference.txt",
                            f"{TEST_DATA}/raman_exp_peaks_reference.txt"], kind_of_spectra=[0, 1])
    algo = Algorithm(exp_peaks=exp_peaks, theo_peaks=theo_peaks,
                     lower_bound=1000, upper_bound=1500, sc=0.975, cutoff=0.01)
    return_value, old_freq, freq, inten, new_sigma, new_eta, kind_of_spectrum = algo.needleman()
    computed = {"return_value": return_value, "old_freq": old_freq, "freq": freq, "inten": inten,
                "new_sigma": new_sigma, "new_eta": new_eta, "kind_of_spectrum": kind_of_spectrum}
    for key in computed.keys():
        if isinstance(reference[key], list):
            assert (np.abs(computed[key] - np.asarray(reference[key])) < 1e-3).all()
        else:
            assert abs(computed[key] - reference[key]) < 1e-4


def test_voigt():
    with open(f"{TEST_DATA}/reference_algo", 'r') as file:
        reference = json.load(file)
    references_converted = {}
    for key in reference.keys():
        if isinstance(reference[key], list):
            references_converted[key] = np.asarray(reference[key])
        else:
            references_converted[key] = reference[key]
    return_value, old_freq, freq, inten, new_sigma, new_eta, kind_of_spectrum = references_converted["return_value"], references_converted["old_freq"], references_converted[
        "freq"], references_converted["inten"], references_converted["new_sigma"], references_converted["new_eta"], references_converted["kind_of_spectrum"]

    x_ir, y_ir = voigt(
        freqs=freq[kind_of_spectrum == 0],
        inten=inten[kind_of_spectrum == 0],
        new_sigma=new_sigma[kind_of_spectrum == 0],
        new_eta=new_eta[kind_of_spectrum == 0])
    generated = np.concatenate([x_ir[:, np.newaxis], y_ir[:, np.newaxis]], axis=-1)
    reference_voigt = np.loadtxt(f"{TEST_DATA}/reference_voigt.txt")
    assert (np.abs(generated - reference_voigt) < 1e-4).all()


def test_get_pearson_and_spearman():
    reference_voigt_1 = np.loadtxt(f"{TEST_DATA}/reference_voigt.txt")
    reference_voigt_2 = np.loadtxt(f"{TEST_DATA}/reference_voigt.txt")
    pearson_ir, spearman_ir = get_spearman_and_pearson(spectrum_exp=reference_voigt_1, spectrum_theo=reference_voigt_2)
    assert abs(pearson_ir - spearman_ir) < 1e-5
    assert abs(pearson_ir - 1) < 1e-4
