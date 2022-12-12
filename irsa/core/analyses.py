import os
from typing import Tuple

import numpy as np
import scipy
from scipy.interpolate import interp1d


def _interpolate(spectrum: np.array, lower: float, upper: float, kind: str = "cubic") -> np.array:
    """Helper function for generating interpolated spectra

    Args:
        spectrum (np.array):
        lower (float):
        upper (float):
        kind (str): Defaults to "cubic".
    Returns:
        np.array: The interpolated spectrum.
    """
    x_range = np.linspace(lower, upper, endpoint=True, num=upper - lower)[10:-10]
    idx = (spectrum[:, 0] > lower) & (spectrum[:, 0] < upper)
    spectrum = spectrum[idx]
    spectrum[:, 1] /= np.max(np.abs(spectrum[:, 1]))
    return interp1d(spectrum[:, 0], spectrum[:, 1], kind=kind)(x_range)


def get_spearman_and_pearson(
        spectrum_exp: np.array, spectrum_theo: np.array, lower: float = 1000, upper: float = 1500) -> Tuple[
        float, float]:
    """Returns the spearman and pearson coefficient of two spectra

    Args:
        spectrum_exp (np.array): The experimental spectrum.
        spectrum_theo (np.array): The theoretical spectrum
        lower (float, optional):  Defaults to 1000.
        upper (float, optional): Defaults to 1500.

    Returns:
        Tuple[float, float]: The pearson and spearmanr
    """
    spectrum_exp_interpolated = _interpolate(spectrum=spectrum_exp, lower=lower, upper=upper)
    spectrum_theo_interpolated = _interpolate(spectrum=spectrum_theo, lower=lower, upper=upper)
    pearsonr = scipy.stats.pearsonr(spectrum_exp_interpolated, spectrum_theo_interpolated)[0]
    spearmanr = scipy.stats.spearmanr(spectrum_exp_interpolated, spectrum_theo_interpolated)[0]
    return (pearsonr, spearmanr)
