from collections import Counter

# from sklearn.metrics import precision_recall_curve
# from sklearn.metrics import auc

import math
import numpy as np
# import pandas as pd
import random as rnd
import matplotlib.pyplot as plt

from classes.indel_event import IndelEvent


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


def calc_sequence_log_prob(sequence: list[int], prob_array: list[float]) -> float:  # task 2 / st 1
    prob = 0
    for digit in sequence:
        prob += math.log10(prob_array[digit - 1])
    return round(prob, 4)




