import os
from typing import List

import numpy as np


def load_spectrum(path: str) -> np.array:
    """Loads a spectrum from a file with x, y coordinates (two columns, nothing more).

    Args:
        path (str): Path to file

    Returns:
        np.array: The spectrum loaded in a numpy array.
    """
    return np.loadtxt(path)


def normalize_spectrum(spectrum: np.array, lower: float = 1000, upper: float = 1500) -> np.array:
    """Returns the normalized and filtered spectrum

    Args:
        spectrum (np.array): The spectrum to normalize.
        lower (float, optional): Defaults to 1000.
        upper (float, optional): Defaults to 1500.

    Returns:
        np.array: The normalized spectrum.
    """
    idx = (spectrum[:, 0] > lower) & (spectrum[:, 0] < upper)
    spectrum[:, 1] = spectrum[:, 1] / np.max(np.abs(spectrum[idx, 1]))
    return spectrum[idx]


def load_peaks(paths: List[str], kind_of_spectra: List[int]) -> np.array:
    """Loads the deconvoluted peaks and concatenates them

    Args:
        path (str): Path to peaks file.
        kind_of_spectrum List[int]: kind of spectrum (in the order of paths) to be loaded.

    Returns:
        np.array: The peaks of the spectrum
    """
    concatenate = None
    if len(paths) != len(kind_of_spectra):
        raise ValueError(f"{len(paths)} must be equal to {len(kind_of_spectra)}")
    for idx, path in enumerate(paths):
        if os.path.exists(path):
            peaks = np.loadtxt(path)
            peaks[:, 0] /= np.max(np.abs(peaks[:, 0]))
            peaks[:, 3] = kind_of_spectra[idx]
            if idx == 0:
                concatenate = peaks
            else:
                concatenate = np.concatenate([concatenate, peaks], axis=0)
        else:
            raise FileNotFoundError(f"The file {path} does not exist")
    if concatenate is None:
        raise ValueError("No peaks loaded")
    return concatenate
