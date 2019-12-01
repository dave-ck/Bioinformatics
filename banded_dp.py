import time

import numpy as np
from banded_matrix import MatrixBand


# IMPORTANT: refactor and verify function signatures

def banded_dynprog(alphabet, scoring_matrix, seq_s, seq_t, max_diff, min_diff, show_alignment=True):
    scoring_matrix = np.array(scoring_matrix)
    len_s, len_t = len(seq_s), len(seq_t)
    backtrack = MatrixBand(max_diff, min_diff, len_t + 1, dtype="str",
                           default="E")  # points U(p), L(eft), D(iagonal), or E(nd)
    matrix_v = MatrixBand(max_diff, min_diff, len_t + 1, dtype="int32",
                          default=0)  # because doing local alignment, no negative values arise naturally
    # todo: change -1 default to 0 default once done with debug
    for i in range(1, len_s + 1):
        for j in range(max(i - max_diff, 1), min(i - min_diff, len_t)+1):
            # every (i, j) to reach this point satisfies min_diff <= i-j <= max_diff
            best = 0
            reset = True
            # # always exists
            from_diagonal = matrix_v.get(i - 1, j - 1) + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            if i - j + 1 <= max_diff:
                from_up = matrix_v.get(i, j - 1) + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            else:
                from_up = False
            if i - 1 - j >= min_diff:
                from_left = matrix_v.get(i - 1, j) + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            else:
                from_left = False
            if from_diagonal and best <= from_diagonal:
                best = from_diagonal
                backtrack.set(i, j, "D")
                reset = False
            if from_left and best <= from_left:
                best = from_left
                backtrack.set(i, j, "U")
                reset = False
            if from_up and best <= from_up:
                best = from_up
                backtrack.set(i, j, "L")
                reset = False
            if reset:
                backtrack.set(i, j, "E")
            matrix_v.set(i, j, best)
    alignment_s, alignment_t = [], []
    show_s, show_t = "", ""
    i, j = matrix_v.max_value_coords()
    best_score = matrix_v.get(i, j)
    pointer = backtrack.get(i, j)
    turn = 0
    while pointer != "E":
        turn += 1
        if pointer == "D":
            j -= 1
            i -= 1
            alignment_s.insert(0, i)
            alignment_t.insert(0, j)
            show_s = seq_s[i] + show_s
            show_t = seq_t[j] + show_t
            pointer = backtrack.get(i, j)
        elif pointer == "L":
            j -= 1
            show_t = seq_t[j] + show_t
            show_s = " " + show_s
            pointer = backtrack.get(i, j)
        elif pointer == "U":
            i -= 1
            show_s = seq_s[i] + show_s
            show_t = " " + show_t
            pointer = backtrack.get(i, j)
    if show_alignment:
        display_alignment(show_s, show_t)
    return [best_score, alignment_s, alignment_t]


# one-line function, keep for cleaner code
def score(alphabet, scoring_matrix, char1, char2):
    return scoring_matrix[(char1 == "_") * -1 or alphabet.index(char1)][(char2 == "_") * -1 or alphabet.index(char2)]


def display_alignment(string1, string2):    # todo: refactor to check alignment value and take lists of indices as input
    string3 = ''
    for i in range(min(len(string1), len(string2))):
        if string1[i] == string2[i]:
            string3 = string3 + "|"
        else:
            string3 = string3 + " "
    print('Alignment ')
    print('String1: ' + string1)
    print('         ' + string3)
    print('String2: ' + string2 + '\n\n')


# test code below
#
#
# sigma = "ABCD"
# score_matrix = [[1, -5, -5, -5, -1],
#                 [-5, 1, -5, -5, -1],
#                 [-5, -5, 5, -5, -4],
#                 [-5, -5, -5, 6, -4],
#                 [-1, -1, -4, -4, -9]]
#
# diff = 6
# string1 = "AAAAAABBADADBBAD"
# string2 =       "BBDADABBAD" # todo: investigate why fucky stuff happens when switched out for "BCDBBAD"
# a = banded_dynprog(sigma, score_matrix, string1, string2, diff + 1, diff - 1)
# print("String1: {}, String2: {}".format(string1, string2))
