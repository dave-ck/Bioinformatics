import math
import submit_main
import numpy as np


# IMPORTANT: refactor and verify function signatures

def dynproglin(alphabet, scoring_matrix, seq_s, seq_t, show_alignment=False):
    # todo: chop off to make local
    score, start_ptr, end_ptr = get_local_score_and_endpoints(alphabet, scoring_matrix, seq_s, seq_t)
    start_s, start_t = start_ptr
    end_s, end_t = end_ptr
    print("Calling hirschberg({}, {}, {}, {})".format(alphabet, score_matrix, seq_s[start_s:end_s], seq_t[start_t:end_t]))
    indices_s, indices_t = hirschberg(alphabet, score_matrix, seq_s[start_s:end_s], seq_t[start_t:end_t])
    # res = submit_main.dynprog(alphabet, score_matrix, seq_s[start_s:end_s], seq_t[start_t:end_t])
    # indices_s, indices_t = res[1], res[2]
    indices_s = list(map(lambda x: x + start_s, indices_s))
    indices_t = list(map(lambda x: x + start_t, indices_t))
    return score, indices_s, indices_t


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
    if len(seq_s) > len(seq_t):
        seq_s, seq_t = seq_t, seq_s
        swapped = True  # keep track of whether s and t were swapped around to allow wlog assumption len_t > len_s
    else:
        swapped = False
    print("Called on s={}, t={}".format(seq_s, seq_t))
    if len(seq_s) == 1 or len(seq_t) == 1:
        alignment_s, alignment_t = needleman_wunsch_char_v_seq(alphabet, scoring_matrix, seq_s, seq_t)
    else:
        len_t = len(seq_t)
        t_mid = math.floor(len_t / 2)
        print("Running global on: {}, {}".format(seq_s, seq_t[t_mid:]))
        score_left = get_global_score(alphabet, scoring_matrix, seq_s, seq_t[:t_mid])
        rev_s = seq_s[::-1]
        seq_t_pt2 = seq_t[t_mid:]
        rev_t_pt2 = seq_t_pt2[::-1]
        print("Running global on: {}, {}".format(rev_s, seq_t[t_mid:]))
        score_right = get_global_score(alphabet, scoring_matrix, rev_s, rev_t_pt2)
        rev_right = list(reversed(score_right))
        sum_scores = score_left + rev_right
        print(score_left)
        print(rev_right)
        print(sum_scores)
        s_mid = np.where(sum_scores == np.amax(sum_scores))[0][0]
        print(s_mid)
        # todo: fix relative addressing
        alignment_s_1, alignment_t_1 = hirschberg(alphabet, scoring_matrix, seq_s[:s_mid], seq_t[:t_mid])
        alignment_s_2, alignment_t_2 = hirschberg(alphabet, scoring_matrix, seq_s[s_mid:], seq_t[t_mid:])
        alignment_s_2 = list(map(lambda x: x + s_mid, alignment_s_2))
        alignment_t_2 = list(map(lambda x: x + t_mid, alignment_t_2))
        alignment_s, alignment_t = alignment_s_1 + alignment_s_2, alignment_t_1 + alignment_t_2
    if swapped:
        alignment_t, alignment_s = alignment_s, alignment_t
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
            best_for_cell = 0
            reset = True
            from_diagonal = matrix_v[i - 1, 0] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            from_left = matrix_v[i, 0] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_up = matrix_v[i - 1, 1] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            if best_for_cell <= from_diagonal:
                best_for_cell = from_diagonal
                backtrack[i, 1] = np.copy(backtrack[i - 1, 0])
                matrix_v[i, 1] = best_for_cell
                reset = False
            if best_for_cell <= from_left:
                best_for_cell = from_left
                backtrack[i, 1] = np.copy(backtrack[i, 0])
                matrix_v[i, 1] = best_for_cell
                reset = False
            if best_for_cell <= from_up:
                best_for_cell = from_up
                backtrack[i, 1] = np.copy(backtrack[i - 1, 1])
                matrix_v[i, 1] = best_for_cell
                reset = False
            if reset:
                matrix_v[i, 1] = 0
                backtrack[i, 1] = (i, j)
            if best_for_cell > best_score:
                # update pointers to start & end for best local alignment
                best_score = best_for_cell
                best_ptr_end = (i, j)
                # print("Backtrack:\n", backtrack)
                best_ptr_start = tuple(np.copy(backtrack[i, 1]))  # need to apply tuple to avoid pass-by-reference
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


def score_alignment(alphabet, scoring_matrix, alignment_s, alignment_t, seq_s, seq_t):
    if not seq_t and not seq_s:
        return 0
    current_score = 0
    for i, j in zip(alignment_s, alignment_t):
        current_score += score(alphabet, scoring_matrix, seq_s[i], seq_t[j])  # score of all subs
    for index in range(len(alignment_s) - 1):
        i = alignment_s[index]
        i_p = alignment_s[index + 1]
        j = alignment_t[index]
        j_p = alignment_t[index + 1]
        for char in seq_s[i + 1:i_p]:
            current_score += score(alphabet, scoring_matrix, char, "_")  # score of all s_dels
        for char in seq_t[j + 1:j_p]:
            current_score += score(alphabet, scoring_matrix, char, "_")  # score of all t_dels
    return current_score


# test code below

# From DUO example
sigma = "ABCD"
score_matrix = [[1, -5, -5, -5, -1],
                [-5, 1, -5, -5, -1],
                [-5, -5, 5, -5, -4],
                [-5, -5, -5, 6, -4],
                [-1, -1, -4, -4, -9]]

seq_s = "CDC"
seq_t = "CAAADC"
# adding   ^ "A" causes issues
# todo: check if issues still arise when dynprog (simple dynamic) is used instead of Hirschberg; shouldn't...


# from wikipedia example: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
# sigma = "ACGT"
# score_matrix = [[2, -1, -1, -1, -2],
#                 [-1, 2, -1, -1, -2],
#                 [-1, -1, 2, -1, -2],
#                 [-1, -1, -1, 2, -2],
#                 [-2, -2, -2, -2, -2]]
#
# seq_s = "AACCCC"
# seq_t = "CAACC"

# High-score segments of perfect matches separated by short segments of cheap indels.
# a = get_local_score_and_endpoints(sigma, score_matrix, seq_s, seq_t)
a = dynproglin(sigma, score_matrix, seq_s, seq_t)
# a = needleman_wunsch_char_v_seq(sigma, score_matrix, "C", "C")
print("Yeet")
print(a)
print(seq_s)
print(seq_t)
alignment_s = a[1]
alignment_t = a[2]

b = score_alignment(sigma, score_matrix, alignment_s, alignment_t, seq_s, seq_t)
print(b)
