from cmath import inf
from typing import Tuple

import numpy as np


def backtrace(p_mat: np.array,
              n_theo_peaks: int,
              n_exp_peaks: int,
              freq_i: np.array,
              inten_i: np.array,
              exp_freq_j: np.array,
              sigma: np.array,
              vcd: np.array,
              eta: np.array) -> Tuple[np.array]:
    """Backtraces the ir alignment. Assigns unassigned peaks to new position.

    Args:
        p_mat (np.array):
        n_theo_peaks (int):
        n_exp_peaks (int):
        freq_i (np.array):
        inten_i (np.array):
        exp_freq_j (np.array):
        sigma (np.array):
        vcd (np.array):
        eta (np.array):

    Returns:
        Tuple[np.array]: Tuple of np.arrays with information about the aligned spectrum.
    """
    new_freq, old_freq, new_inten, new_sigma = [], [], [], []
    new_eta, non_matched_sigma, new_inten_vcd, non_matched_freq = [], [], [], []
    matched_freq, vcd_ir_array, non_matched_inten, non_matched_inten_vcd = [], [], [], []
    n_theo_peaks = n_theo_peaks - 1
    n_exp_peaks = n_exp_peaks - 1
    current_scaling_factor = 1
    factors = []
    while True:
        if p_mat[n_theo_peaks, n_exp_peaks] == "D":
            new_freq.append(exp_freq_j[n_exp_peaks - 1])
            old_freq.append(freq_i[n_theo_peaks - 1])
            new_inten.append(inten_i[n_theo_peaks - 1])
            new_sigma.append(sigma[n_theo_peaks - 1])
            new_eta.append(eta[n_theo_peaks - 1])
            vcd_ir_array.append(vcd[n_theo_peaks - 1])
            current_scaling_factor = exp_freq_j[n_exp_peaks - 1] / freq_i[n_theo_peaks - 1]
            matched_freq.append(n_theo_peaks - 1)
            factors.append(current_scaling_factor)
            n_theo_peaks = n_theo_peaks - 1
            n_exp_peaks = n_exp_peaks - 1
        elif p_mat[n_theo_peaks, n_exp_peaks] == "V":
            non_matched_inten.append(n_theo_peaks - 1)
            non_matched_sigma.append(n_theo_peaks - 1)
            non_matched_inten_vcd.append(n_theo_peaks - 1)
            non_matched_freq.append(n_theo_peaks - 1)
            n_theo_peaks = n_theo_peaks - 1
        elif p_mat[n_theo_peaks, n_exp_peaks] == "H":
            n_exp_peaks = n_exp_peaks - 1
        else:
            break

    for i, non_matched_f in enumerate(non_matched_freq):
        closest_distance = inf
        sf = 1
        for j, matched_f in enumerate(matched_freq):
            dis = abs(freq_i[non_matched_f] - freq_i[matched_f])
            if dis < closest_distance:
                closest_distance = dis
                sf = factors[j]
        new_freq.append(freq_i[non_matched_f] * sf)
        new_sigma.append(sigma[non_matched_sigma[i]])
        new_eta.append(eta[non_matched_sigma[i]])
        vcd_ir_array.append(vcd[non_matched_f])
        old_freq.append(freq_i[non_matched_f])
        new_inten.append(inten_i[non_matched_f])
        new_inten_vcd.append(0)
    return (np.asarray(new_freq),
            np.asarray(new_inten),
            np.asarray(old_freq),
            np.asarray(new_sigma),
            np.asarray(new_eta),
            np.asarray(vcd_ir_array))


def pointer(di: float, ho: float, ve: float) -> str:
    """returns a point to the cell with the minimum path.

    Args:
        di (float): diagonal.
        ho (float): horizontal.
        ve (float): vertical.

    Returns:
        str: Minimum path encoded as str.
    """
    point = min(di, min(ho, ve))
    if di == point:
        return "D"
    if ho == point:
        return "H"
    return "V"


class Algorithm:
    """Needleman Wunsch Algorithm Class

    Returns:
        _type_: Object.
    """

    def __init__(self, theo_peaks: np.array, exp_peaks: np.array,
                 cutoff: float = 0.01, lower_bound: float = 1000, upper_bound: float = 1800,
                 sc: float = 0.98, algo: int = 0) -> None:
        """Init function.

        Args:
            theo_peaks (np.array): The theoretical, deconvoluted peaks.
            exp_peaks (np.array): The experimental, deconvoluted peaks.
            cutoff (float, optional): The cutoff applied for the function. Defaults to 0.01.
            lower_bound (float, optional): The lower bound of the spectrum to be considered. Defaults to 1000.
            upper_bound (float, optional): The upper bound of the spectrum to be considered . Defaults to 1800.
            sc (float, optional): The scaling fcator. Defaults to 0.98.
            algo (int, optional):  The algorithm number (scoring function). Defaults to 0.
        """

        self.cutoff, self.theo_peaks, self.exp_peaks = cutoff, theo_peaks, exp_peaks
        self.algo, self.lower_bound, self.upper_bound, self.sc = algo, lower_bound, upper_bound, sc

    # this is the scoring function
    def diagonal(self, freq_i: float, inten_i: float, exp_freq_j: float, exp_inten_j: float, 
                 exp_vcd: float = 0, inten_vcd: float = 0, width_j: float = 0, width_i: float = 0) -> float:
        """The scoring fuction of the algorithm. Implement here new stuff if you want to, and provide a new algo number

        Args:
            freq_i (float): The frequency of peak i.
            inten_i (float): The intensity of peak i.
            exp_freq_j (float): The frequency of the experimental peak j.
            exp_inten_j (float): The intensity of the experimental peak j.
            exp_vcd (float, optional): The frequency of the experimental vcd peak. Defaults to 0.
            inten_vcd (float, optional): The intensity of the vcd peak. Defaults to 0.
            width_j (float, optional): the width of peak j. Defaults to 0.
            width_i (float, optional): The width of peak i. Defaults to 0.

        Returns:
            float: _description_
        """
        value = np.inf
        if self.algo == 0:
            if inten_vcd == exp_vcd:
                tmp_bool_1 = min(abs(1 - exp_freq_j / freq_i),
                                 abs(1 - freq_i / exp_freq_j)
                                 ) < self.cutoff
                tmp_bool_2 = exp_freq_j > self.lower_bound
                tmp_bool_3 = exp_freq_j < self.upper_bound
                if tmp_bool_1 and tmp_bool_2 and tmp_bool_3:
                    x_dummy = min(abs(1 - exp_freq_j / freq_i),
                                  abs(1 - freq_i / exp_freq_j))
                    width_dummy = min(abs(1 - width_j / width_i),
                                      abs(1 - width_i / width_j))
                    freq_contrib = np.exp(-1 / (1 - abs(x_dummy / self.cutoff)**2))
                    y_dummy = min(abs(1 - inten_i / exp_inten_j),
                                  abs(1 - exp_inten_j / inten_i))
                    inten_contrib = np.exp(-1 / (1 - abs(y_dummy / 1)**2))
                    sigma_contrib = np.exp(-1 / (1 - abs(width_dummy / 8)**2))
                    if min(abs(1 - width_i / width_j), abs(1 - width_j / width_i)) < 8:
                        if abs(1 - inten_i / exp_inten_j) < 1 or abs(1 - exp_inten_j / inten_i) < 1:
                            value = -inten_contrib * freq_contrib * sigma_contrib

        return value

    def needleman(self):
        """Performs the needleman Wunsch algorithm for the provided data.

        Returns:
            _type_: A tuple of objects with information one might be interested in.
        """
        idx = np.argsort(self.theo_peaks[:, 1])
        self.theo_peaks = self.theo_peaks[idx]
        freq = self.theo_peaks[:, 1] * self.sc
        inten = self.theo_peaks[:, 0]
        sigma = self.theo_peaks[:, 2]
        vcd = self.theo_peaks[:, 3]
        try:
            eta = self.theo_peaks[:, 4]
        # pylint: disable=broad-except
        except Exception:
            eta = np.ones((len(sigma)))

        idx = np.argsort(self.exp_peaks[:, 1])
        self.exp_peaks = self.exp_peaks[idx]
        exp_freq = self.exp_peaks[:, 1]

        exp_inten = self.exp_peaks[:, 0]
        exp_sigma = self.exp_peaks[:, 2]
        exp_inten_vcd = self.exp_peaks[:, 3]

        n_theo_peaks = len(freq) + 1
        n_exp_peaks = len(exp_freq) + 1
        al_mat = np.zeros((n_theo_peaks, n_exp_peaks))
        p_mat = np.zeros((n_theo_peaks, n_exp_peaks), dtype='U25')  # string
        for i in range(1, n_theo_peaks):
            al_mat[i, 0] = al_mat[i - 1, 0]
            p_mat[i, 0] = 'V'
        for i in range(1, n_exp_peaks):
            al_mat[0, i] = al_mat[0, i - 1]
            p_mat[0, i] = 'H'
        p_mat[0, 0] = "S"
        for i in range(1, n_theo_peaks):  # theoretical
            for j in range(1, n_exp_peaks):  # experimental
                di = self.diagonal(freq[i - 1],
                                   inten[i - 1],
                                   exp_freq[j - 1],
                                   exp_inten[j - 1],
                                   exp_vcd=exp_inten_vcd[j - 1],
                                   inten_vcd=vcd[i - 1],
                                   width_j=exp_sigma[j - 1],
                                   width_i=sigma[i - 1])
                di = al_mat[i - 1, j - 1] + di
                ho = al_mat[i, j - 1]
                ve = al_mat[i - 1, j]
                al_mat[i, j] = min(di, min(ho, ve))
                p_mat[i, j] = pointer(di, ho, ve)
        freq, inten, old_freq, new_sigma, eta_new, vcd = backtrace(p_mat=p_mat,
                                                                   n_theo_peaks=n_theo_peaks,
                                                                   n_exp_peaks=n_exp_peaks,
                                                                   freq_i=freq,
                                                                   inten_i=inten,
                                                                   exp_freq_j=exp_freq,
                                                                   sigma=sigma,
                                                                   vcd=vcd,
                                                                   eta=eta)
        returnvalue = al_mat[n_theo_peaks - 1, n_exp_peaks - 1]
        return -returnvalue, old_freq, freq, inten, new_sigma, np.asarray(eta_new), np.asarray(vcd)
