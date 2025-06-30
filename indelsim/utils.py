import math
from operator import indexOf

import numpy as np
import matplotlib.pyplot as plt


def random_zipf(a: float, max_int: int) -> int:
    am1: float = a - 1.0
    b: float = pow(2.0, am1)
    while 1:
        U: float = 1.0 - np.random.random(1)
        V: float = np.random.random(1)
        X: float = math.floor(pow(U, -1 / am1))
        if X > max_int or X < 1:
            continue
        T: float = pow(1 + 1 / X, am1)
        if V * X * (T - 1) / (b - 1) <= T / b:
            return int(X)


def calc_trunc_zipf(alpha: float, max_val: int) -> int:
    while True:
        z: int = np.random.zipf(alpha)
        if z <= max_val:
            return z


def plot_distribution(distribution_list: list[float], bins: int, density: bool, file_name: str,
                      measure: str):  # task 1 / st 5
    fig, ax = plt.subplots()
    ax.hist(distribution_list, bins=bins, density=density)
    ax.set_title(f'{measure} Distribution')
    plt.xlabel(f'{measure} Value')
    plt.ylabel("counts")
    plt.savefig(file_name)
    plt.show()


def calc_msa_from_naive_nodes(sequences: list[list[int]], ancestors: list[int]) -> list[list[int]]:
    msa: list[list[int]] = [sequences[0]]
    for seq_inx in range(1, len(sequences)):
        msa.append([-1] * len(sequences[0]))
        ancestor_inx: int = ancestors[seq_inx]
        last_msa_inx: int = -1
        for char_inx in sequences[seq_inx]:
            try:
                char_inx_on_msa: int = msa[ancestor_inx].index(char_inx)
                last_msa_inx = char_inx_on_msa
            except ValueError:
                char_inx_on_msa: int = last_msa_inx + 1
                for s in msa:
                    s[char_inx_on_msa:char_inx_on_msa] = [-1]
                last_msa_inx = char_inx_on_msa
            msa[seq_inx][char_inx_on_msa] = char_inx
    return msa


def get_msa_as_str_list(msa: list[list[int]], padding: int) -> list[str]:
    res: list[str] = []
    for seq in msa:
        res.append(','.join([str(i).rjust(padding, ' ') for i in seq]))
    return res



