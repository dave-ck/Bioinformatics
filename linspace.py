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


def get_score_and_endpoints(alphabet, scoring_matrix, seq_s, seq_t):
    scoring_matrix = np.array(scoring_matrix)
    len_s, len_t = len(seq_s), len(seq_t)
    matrix_v = np.zeros([len_s + 1, 2], dtype=int)  # because doing local alignment, no negative values
    backtrack = np.empty([len_s + 1, 2], dtype=(int, 2))  # pointer to last "reset"
    backtrack[:, :] = (0, 0)
    best_score = 0
    best_ptr_start, best_ptr_end = None, None
    for j in range(1, len_t + 1):
        for i in range(1, len_s + 1):
            # fill column j (held at 1) using column j-1 (held at 0)
            best_for_cell = 0
            reset = True
            from_diagonal = matrix_v[i - 1, 0] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            from_up = matrix_v[i, 0] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_left = matrix_v[i - 1, 1] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            if best_for_cell <= from_diagonal:
                best_for_cell = from_diagonal
                backtrack[i, 1] = backtrack[i - 1, 0]
                matrix_v[i, 1] = best_for_cell
                reset = False
            if best_for_cell <= from_left:
                best_for_cell = from_left
                backtrack[i, 1] = backtrack[i, 0]
                matrix_v[i, 1] = best_for_cell
                reset = False
            if best_for_cell <= from_up:
                best_for_cell = from_up
                backtrack[i, 1] = backtrack[i - 1, 1]
                matrix_v[i, 1] = best_for_cell
                reset = False
            if reset:
                matrix_v[i, 1] = 0
                backtrack[i, 1] = (i, j)
            if best_for_cell > best_score:
                # update pointers to start & end for best local alignment
                best_score = best_for_cell
                best_ptr_end = (i, j)
                best_ptr_start = tuple(backtrack[i, 1])  # need to apply tuple to avoid pass-by-reference
                if best_ptr_start[0] == 0 and best_ptr_start[1] == 0:
                    print("Hi")
                else:
                    print("Howdy", best_ptr_start)
        print(backtrack)
        # switch out column 1 for column 0
        matrix_v[:, [0, 1]] = matrix_v[:, [1, 0]]
        backtrack[:, [0, 1]] = backtrack[:, [1, 0]]
    print("Howdyado", best_ptr_start)
    return [best_score, best_ptr_start, best_ptr_end]


# one-line function, keep for cleaner code
def score(alphabet, scoring_matrix, char1, char2):
    return scoring_matrix[(char1 == "_") * -1 or alphabet.index(char1)][(char2 == "_") * -1 or alphabet.index(char2)]


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


sigma = "ABCD"
score_matrix = [[1, -5, -5, -5, -1],
                [-5, 1, -5, -5, -1],
                [-5, -5, 5, -5, -4],
                [-5, -5, -5, 6, -4],
                [-1, -1, -4, -4, -9]]

a = get_score_and_endpoints(sigma, score_matrix, "BBBBAACAA", "AAACBB")
print(a)
