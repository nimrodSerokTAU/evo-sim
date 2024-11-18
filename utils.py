from collections import Counter

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

import math
import numpy as np
import pandas as pd
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


def present_prob_as_log(probability: float) -> float:  # task 2 / st 1
    return round(math.log10(probability), 4)


def is_sequence_fair(runs, fair_prob: list[float], biased_prob: list[float]):  # task 2 / st 3
    fair_ll = calc_sequence_log_prob(runs, fair_prob)
    biased_ll = calc_sequence_log_prob(runs, biased_prob)
    return 1 if fair_ll > biased_ll else 0


def is_counter_fair(counter: list[int] or np.array, fair_prob: list[float], biased_prob: list[float]):  # task 2 / st 3
    fair_ll = calc_counter_log_prob(counter, fair_prob)
    biased_ll = calc_counter_log_prob(counter, biased_prob)
    return 1 if fair_ll > biased_ll else 0


def calc_counter_log_prob(counter: list[int], prob_array: list[float]) -> float:  # task 2 / st 2
    prob = 0
    for i in range(1, 7):
        prob += math.log10(prob_array[i - 1]) * counter[i - 1]
    return round(prob, 4)


def read_matching_matrix(file_path: str) -> tuple[list[list[int]], dict[str, int]]:
    codes_dict_to_inx: dict[str, int] = {}
    match_matrix: list[list[int]] = []
    with open(file_path) as infile:
        for l_inx, line in enumerate(infile):
            line = line.strip()
            is_first_col = True
            if l_inx == 0:
                for c in line.split(' '):
                    if c != '':
                        codes_dict_to_inx[c] = len(codes_dict_to_inx)
            else:
                match_matrix.append([])
                for c in line.split(' '):
                    if c != '':
                        if not is_first_col:
                            match_matrix[l_inx - 1].append(int(c))
                        is_first_col = False
    return match_matrix, codes_dict_to_inx


def calc_inserted_seq(original_seq_len: int, alpha: float) -> tuple[list[str], int]:
    codons: list[str] = ['a', 't', 'g', 'c']
    place: int = rnd.randint(0, original_seq_len)
    insertion_size: int = calc_trunc_zipf(alpha, 50)
    insertion: list[str] = []
    for i in range(insertion_size):
        insertion.append(codons[rnd.randint(0, len(codons) - 1)])
    return insertion, place


def calc_inserted_seq_len_and_place(original_seq_len: int, alpha: float) -> IndelEvent:
    place: int = rnd.randint(0, original_seq_len)
    insertion_size: int = calc_trunc_zipf(alpha, 50)
    return IndelEvent(length=insertion_size, place=place, is_insertion=True)


def insertion_event(seq: list[str], alpha: float) -> list[str]:
    insertion, place = calc_inserted_seq(len(seq), alpha)
    seq[place:place] = insertion
    return seq


