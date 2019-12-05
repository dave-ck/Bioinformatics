import math

import numpy as np


# IMPORTANT: refactor and verify function signatures

def dynproglin(alphabet, scoring_matrix, seq_s, seq_t, show_alignment=False):
    # run through, find end points

    # initialize 2 columns (or rows) each for matrix_v and backtrack; these alternate at each step
    # backtrack contains tuples pointing to sequence start
    # once run-through is complete, read endpoints for the optimal local alignment from pointers
    # run linear space global alignment on the cropped (spliced) sequences

    # re-adjust the output of global alignment by adding the starting coordinates to each index
    # so i.e. "ABB" and "BBA" would produce [2, [0,1], [0,1]]
    # for global alignment of "BB" and "BB" - make it [2, [1,2], [0,1]]
    pass


# one-line function, keep for cleaner code
def score(alphabet, scoring_matrix, char1, char2):
    return scoring_matrix[(char1 == "_") * -1 or alphabet.index(char1)][(char2 == "_") * -1 or alphabet.index(char2)]


def needleman_wunsch_char_v_seq(alphabet, scoring_matrix, seq_s, seq_t):
    if len(seq_t) != 1:
        seq_s, seq_t = seq_t, seq_s
        swapped = True  # keep track of whether s and t were swapped around to allow wlog assumption len_t == 1
    else:
        swapped = False
    len_s, len_t = len(seq_s), len(seq_t)
    if len_s != 1 and len_t != 1:
        raise ValueError("At least one of the sequences passed must comprise exactly one character.")
    char_t = seq_t  # semantics
    indel_char_t_score = score(alphabet, scoring_matrix, char_t, "_")
    all_indel_score = indel_char_t_score
    for i in range(len_s):
        all_indel_score += score(alphabet, scoring_matrix, seq_s[i], "_")
    best_score = all_indel_score
    print("All indel:", all_indel_score)
    best_alignment_s = []
    best_alignment_t = []
    for i in range(len_s):
        # score = score of full indel - score for indel of s_i - score for indel of t + score for sub(s_i, t)
        current_score = all_indel_score \
                        - indel_char_t_score \
                        - score(alphabet, scoring_matrix, seq_s[i], "_") \
                        + score(alphabet, scoring_matrix, char_t, seq_s[i])
        if current_score > best_score:
            best_score = current_score
            best_alignment_s = [i]
            best_alignment_t = [0]
    if swapped:
        best_alignment_s, best_alignment_t = best_alignment_t, best_alignment_s
    return best_alignment_s, best_alignment_t


def hirschberg(alphabet, scoring_matrix, seq_s, seq_t):
    if len(seq_t) != 1:
        seq_s, seq_t = seq_t, seq_s
        swapped = True  # keep track of whether s and t were swapped around to allow wlog assumption len_t == 1
    else:
        swapped = False
    alignment_s = []
    alignment_t = []
    if len(seq_s) == 1 or len(seq_t) == 1:
        alignment_s, alignment_t = needleman_wunsch_char_v_seq(alphabet, scoring_matrix, seq_s, seq_t)
    else:
        xlen = len(X)
        xmid = math.floor(len(X) / 2)
        ylen = len(Y)
        scoreL = get_global_score(alphabet, scoring_matrix, X[:xmid], Y)
        revY = Y[::-1]
        scoreR = get_global_score(alphabet, scoring_matrix, X[xmid:], revY)
        sum_scores = scoreL + scoreR
        print(np.where(sum_scores == np.amax(sum_scores))[0][0])
        # ymid = arg max ScoreL + rev(ScoreR)
        ymid = np.where(sum_scores == np.amax(sum_scores))[0][0]
        Z1, W1 = hirschberg(alphabet, scoring_matrix, X[:xmid], Y[:ymid])
        Z2, W2 = hirschberg(alphabet, scoring_matrix, X[xmid:], Y[ymid:])
        Z, W = Z1 + Z2, W1 + W2
    return alignment_s, alignment_t


def get_local_score_and_endpoints(alphabet, scoring_matrix, seq_s, seq_t):
    scoring_matrix = np.array(scoring_matrix)
    len_s, len_t = len(seq_s), len(seq_t)
    matrix_v = np.zeros([len_s + 1, 2], dtype=int)  # because doing local alignment, no negative values
    backtrack = np.empty([len_s + 1, 2], dtype=(int, 2))  # pointer to last "reset"
    backtrack[:, :] = (0, 0)
    for i in range(1, len_s + 1):
        backtrack[i, 0] = (i, 0)
    best_score = 0
    best_ptr_start, best_ptr_end = None, None
    for j in range(1, len_t + 1):
        backtrack[0, 1] = (0, j)
        for i in range(1, len_s + 1):
            # fill column j (held at 1) using column j-1 (held at 0)
            print()
            print("i, j: {}, {}".format(i, j))
            print("seq_s, seq_t: {}, {}".format(seq_s, seq_t))
            print("seq_s[i-1:i+1], seq_t[j-1:j+1]: {}, {}".format(seq_s[i - 1:i + 1], seq_t[j - 1:j + 1]))
            best_for_cell = 0
            reset = True
            from_diagonal = matrix_v[i - 1, 0] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            from_left = matrix_v[i, 0] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_up = matrix_v[i - 1, 1] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            if best_for_cell <= from_diagonal:
                print("From diagonal.")
                best_for_cell = from_diagonal
                print(backtrack)
                backtrack[i, 1] = np.copy(backtrack[i - 1, 0])
                matrix_v[i, 1] = best_for_cell
                reset = False
            if best_for_cell <= from_left:
                print("From left.")
                best_for_cell = from_left
                backtrack[i, 1] = np.copy(backtrack[i, 0])
                matrix_v[i, 1] = best_for_cell
                reset = False
            if best_for_cell <= from_up:
                print("From up.")
                best_for_cell = from_up
                backtrack[i, 1] = np.copy(backtrack[i - 1, 1])
                matrix_v[i, 1] = best_for_cell
                reset = False
            if reset:
                print("Reset.")
                matrix_v[i, 1] = 0
                backtrack[i, 1] = (i, j)
            if best_for_cell > best_score:
                # update pointers to start & end for best local alignment
                best_score = best_for_cell
                best_ptr_end = (i, j)
                # print("Backtrack:\n", backtrack)
                best_ptr_start = tuple(np.copy(backtrack[i, 1]))  # need to apply tuple to avoid pass-by-reference
            print()
        # switch out column 1 for column 0
        matrix_v[:, [0, 1]] = matrix_v[:, [1, 0]]
        backtrack[:, [0, 1]] = backtrack[:, [1, 0]]
    return [best_score, best_ptr_start, best_ptr_end]


def get_global_score(alphabet, scoring_matrix, seq_s, seq_t):
    scoring_matrix = np.array(scoring_matrix)
    len_s, len_t = len(seq_s), len(seq_t)
    matrix_v = np.zeros([len_s + 1, 2], dtype=int)
    matrix_v[0, 0] = 0
    for i in range(1, len_s + 1):
        matrix_v[i, 0] = matrix_v[i - 1, 0] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
    for j in range(1, len_t + 1):
        matrix_v[0, 1] = matrix_v[0, 0] + score(alphabet, scoring_matrix, seq_t[j - 1], "_")
        for i in range(1, len_s + 1):
            # fill column j (held at 1) using column j-1 (held at 0)
            from_diagonal = matrix_v[i - 1, 0] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            from_up = matrix_v[i, 0] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_left = matrix_v[i - 1, 1] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            best_for_cell = max(from_diagonal, from_up, from_left)
            matrix_v[i, 1] = best_for_cell
        # switch out column 1 for column 0
        matrix_v[:, [0, 1]] = matrix_v[:, [1, 0]]
    return matrix_v[:, 0]  # return the last column computed, now at index 0


def display_alignment(string1, string2):
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

# From DUO example
sigma = "ABCD"
score_matrix = [[1, -5, -5, -5, -1],
                [-5, 1, -5, -5, -1],
                [-5, -5, 5, -5, -4],
                [-5, -5, -5, 6, -4],
                [-1, -1, -4, -4, -9]]

seq_s = "CC"
seq_t = "ACC"

# from wikipedia example: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
# sigma = "ACGT"
# score_matrix = [[2, -1, -1, -1, -2],
#                 [-1, 2, -1, -1, -2],
#                 [-1, -1, 2, -1, -2],
#                 [-1, -1, -1, 2, -2],
#                 [-2, -2, -2, -2, -2]]
#
# seq_s = "AGTACGCA"
# seq_t = "TATGC"

# todo: check that endpoints are correctly identified: issue example from DUO not consistent with indices returned
# High-score segments of perfect matches separated by short segments of cheap indels.
a = get_local_score_and_endpoints(sigma, score_matrix, seq_s, seq_t)
print("Yeet")
# a = hirschberg(sigma, score_matrix, seq_s, seq_t)
# a = needleman_wunsch_char_v_seq(sigma, score_matrix, "C", "C")

print(a)
