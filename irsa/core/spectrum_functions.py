import os
from pathlib import Path
from typing import Tuple

import numpy as np
import scipy
from lmfit import Model
from scipy import signal

WORKING_DIR = f"{os.path.dirname(__file__)}/../data/"


def voigt_fn(x: float, amp: float, cen: float, wid: float, eta: float) -> np.array:
    """Places a voigt function

    Args:
        x (float): x value on the graph.
        amp (float): The amplitute of the peak.
        cen (float): The center of the peak.
        wid (float): The width of the peak.
        eta (float): The mixing parameter of the gaussian the lorentzian.

    Returns:
        np.array: The voight function
    """
    tmp = ((x - cen) / (wid / 2))**2
    gauss = 1 * np.exp(-np.log(2) * tmp)
    lorentian = 1 / (1 + tmp)
    voigt_new = amp * (eta * lorentian + (1 - eta) * gauss)
    return voigt_new


def add_peak(prefix: str,
             center: float,
             amplitude: float = 0.5,
             sigma: float = 12,
             eta: float = 1,
             vcd: bool = False) -> Tuple[Model, dict]:
    """Adds a peak to the model.

    Args:
        prefix (str): The prefix name of the peak.
        center (float): The center of the peak.
        amplitude (float, optional): The amplitude of the peak. Defaults to 0.5.
        sigma (float, optional): The width of the peak. Defaults to 12.
        eta (float, optional): The mixing parameter of the peak. Defaults to 1.
        vcd (bool, optional): Wether it is a VCD peak or not. Defaults to False.

    Returns:
        Tuple[Model, dict]: The model and the parameters.
    """
    peak = Model(voigt_fn, prefix=prefix)
    pars = peak.make_params()

    pars[prefix + 'cen'].set(center, min=center - 2, max=center + 2, vary=True)
    if vcd is False:
        pars[prefix + 'amp'].set(amplitude, vary=True, min=0.00, max=1.1)
    else:
        pars[prefix + 'amp'].set(amplitude, vary=True, min=-1.1, max=1.1)

    pars[prefix + 'wid'].set(sigma, vary=True, min=1, max=20)
    pars[prefix + 'eta'].set(eta, vary=True, min=0, max=1)
    return peak, pars


def lorentzian_broadening(peaks: np.array, width=12) -> np.array:
    """Performs a Lorentzian broadening of the spectrum.

    Args:
        peaks (np.array): The peaks array.
        width (int, optional): The assumed peak width.. Defaults to 12.

    Returns:
        np.array: The broadened spectrum.
    """
    p = np.arange(500, 2000)
    x = (p[:, np.newaxis] - (peaks[:, 0])[np.newaxis, :]) / (width / 2)
    lorentzians = (peaks[:, 1])[np.newaxis, :] / (1 + x * x)
    y = np.sum(lorentzians, axis=-1)[:, np.newaxis]
    p = p[:, np.newaxis]
    spectrum = np.concatenate([p, y], axis=-1)
    return spectrum


def voigt(freqs: np.array,
          inten: np.array,
          new_sigma: np.array,
          new_eta: np.array,
          lower: float = 1000,
          upper: float = 1500) -> Tuple[np.array, np.array]:
    """Performs Voigt broadening.

    Args:
        freqs (np.array): The frequencies array (x axis).
        inten (np.array): The intensities of the peaks (y axis).
        new_sigma (np.array): The broadening array.
        new_eta (np.array): The mixing parameters
        lower (float, optional): Lower bound of the spectrum. Defaults to 1000.
        upper (float, optional): The upper bound of the spectrum. Defaults to 1500.

    Returns:
        Tuple[np.array, np.array]: The x axis and the y axis of the spectrum as numpy arrays.
    """
    x_axis = np.arange(lower, upper)
    list_append = []
    for i, freq in enumerate(freqs):
        tmp = ((x_axis - freq) / (new_sigma[i] / 2))**2
        lorentzian = inten[i] / (1 + tmp)
        gauss = inten[i] * np.exp(-np.log(2) * tmp)
        list_append.append(lorentzian * new_eta[i] + (1 - new_eta[i]) * gauss)
    list_append = np.asarray(list_append)
    y_axis = np.sum(list_append, axis=0)
    return x_axis, y_axis

def _get_index_for_peak(peaks_on_x_axis: np.array, 
                        x_dimension_of_spectrum: np.array) -> np.array:
    """Helper function to get peak indices.

    Args:
        peaks_on_x_axis (np.array): Provided x dimension of peaks.
        x_dimension_of_spectrum (np.array): x axis of the spectrum.

    Raises:
        ValueError: (min) If the provided peaks are out of scope.
        ValueError: (max) If the provided peaks are out of scope.

    Returns:
        np.array: indices of the peak.
    """

    if np.min(peaks_on_x_axis) < np.min(x_dimension_of_spectrum):
        raise ValueError("The lowest provided peak is out of scope")
    if np.max(peaks_on_x_axis) > np.max(x_dimension_of_spectrum):
        raise ValueError("The highest provided peak is out of scope")
    distances_between_peaks = np.abs(peaks_on_x_axis[:, np.newaxis] - x_dimension_of_spectrum[np.newaxis, :])

    return np.argmin(distances_between_peaks, axis=1)


def deconvolute(spectrum: np.array,
                working_dir: str = WORKING_DIR,
                save_data: str = 'ir_exp_peaks.txt',
                normalize: bool = False,
                lower: float = 1000,
                higher: float = 1800,
                vcd: bool = False,
                peaks_on_x_axis: np.array = None) -> None:
    """Deconvolutes the spectrum, and saves the information in the text file {save_data}.

    Args:
        spectrum (np.array): The spectrum which should be deconvoluted.
        working_dir (str, optional): The working directory. Defaults to WORKING_DIR.
        save_data (str, optional): Name of the file to be saved. Defaults to 'ir_exp_peaks.txt'.
        normalize (bool, optional): Whether the data should be normalized. Defaults to False.
        lower (float, optional): Lower bound of the spectrum to be considered. Defaults to 1000.
        higher (float, optional): Upper bound of the spectrum to be considered. Defaults to 1800.
        vcd (bool, optional): Whether the data is a VCD spectrum. Defaults to False.
        peaks_on_x_axis (np.array, optional): Provide manuall selected x peaks.
    """
    if normalize:
        idx = (spectrum[:, 0] > lower) & (spectrum[:, 0] < higher)
        spectrum[:, 1] = spectrum[:, 1] / np.max(np.abs(spectrum[idx, 1]))
        spectrum = spectrum[idx]
    if peaks_on_x_axis is None:
        ind_ir, _ = scipy.signal.find_peaks(spectrum[:, 1])
        if vcd:
            ind_ir = []
            for i in range(1, len(spectrum) - 1):
                if abs(spectrum[i - 1, 1]) < abs(spectrum[i, 1]) > abs(spectrum[i + 1, 1]) and abs(spectrum[i, 1]) > 0.03:
                    ind_ir.append(i)
            ind_ir = np.asarray(ind_ir)
    else:
        ind_ir = _get_index_for_peak(peaks_on_x_axis=peaks_on_x_axis,
                                     x_dimension_of_spectrum=spectrum[:, 0])

    peaks = spectrum[ind_ir]
    params, model, write_state = None, None, []
    for i in range(0, len(peaks)):
        peak, pars = add_peak('lz%d_' % (i + 1), center=peaks[i, 0],
                              amplitude=peaks[i, 1], vcd=vcd)
        if i == 0:
            model = peak
            params = model.make_params()
        else:
            model = model + peak
        params.update(pars)

    model.eval(params, x=spectrum[:, 0])
    result = model.fit(spectrum[:, 1], params, x=spectrum[:, 0])
    result.eval_components()

    for _, par in result.params.items():
        write_state.append(par.value)

    write_state = np.asarray(write_state)
    write_state = write_state.reshape(-1, 4)
    with open(Path(f"{working_dir}/{save_data}"), "w") as file:
        for i in write_state:
            file.write(str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + " 0 " + str(i[3]) + "\n")
