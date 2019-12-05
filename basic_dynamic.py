import numpy as np


# IMPORTANT: refactor and verify function signatures

def dynprog(alphabet, scoring_matrix, seq_s, seq_t, show_alignment=False):
    scoring_matrix = np.array(scoring_matrix)
    len_s, len_t = len(seq_s), len(seq_t)
    backtrack = np.empty([len_s + 1, len_t + 1], dtype=str)  # points U(p), L(eft), D(iagonal), or E(nd)
    backtrack[:, :] = "E"
    matrix_v = np.zeros([len_s + 1, len_t + 1], dtype=int)  # because doing local alignment, no negative values
    for i in range(1, len_s + 1):
        for j in range(1, len_t + 1):
            best = 0
            reset = True
            from_diagonal = matrix_v[i - 1, j - 1] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            from_left = matrix_v[i, j - 1] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_up = matrix_v[i - 1, j] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            if best <= from_diagonal:
                best = from_diagonal
                backtrack[i, j] = "D"
                matrix_v[i, j] = best
                reset = False
            if best <= from_left:
                best = from_left
                backtrack[i, j] = "L"
                matrix_v[i, j] = best
                reset = False
            if best <= from_up:
                best = from_up
                backtrack[i, j] = "U"
                matrix_v[i, j] = best
                reset = False
            if reset:
                matrix_v[i, j] = 0  # explicit but redundant - overwrite 0 from np.zeros init
                backtrack[i, j] = "E"
    alignment_s, alignment_t = [], []
    show_s, show_t = "", ""
    best_score = np.amax(matrix_v)
    best_mat = np.where(matrix_v == best_score)
    i, j = best_mat[0][0], best_mat[1][0]
    pointer = backtrack[i, j]
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
            pointer = backtrack[i, j]
        elif pointer == "L":
            j -= 1
            show_s = seq_s[i] + show_s
            show_t = " " + show_t
            pointer = backtrack[i, j]
        elif pointer == "U":
            i -= 1
            show_t = seq_t[j] + show_t
            show_s = " " + show_s
            pointer = backtrack[i, j]
    if show_alignment:
        display_alignment(show_t, show_s)
    return [best_score, alignment_s, alignment_t]


# one-line function, keep for cleaner code
def score(alphabet, scoring_matrix, char1, char2):
    return scoring_matrix[(char1 == "_") * -1 or alphabet.index(char1)][(char2 == "_") * -1 or alphabet.index(char2)]


def display_alignment(show_s, show_t):
    string3 = ''
    for i in range(min(len(show_s), len(show_t))):
        if show_s[i] == show_t[i]:
            string3 = string3 + "|"
        else:
            string3 = string3 + " "
    print('Alignment ')
    print('String s: ' + show_s)
    print('         ' + string3)
    print('String t: ' + show_t + '\n\n')


# test code below


sigma = "ABCD"
score_matrix = [[1, -5, -5, -5, -1],
                [-5, 1, -5, -5, -1],
                [-5, -5, 5, -5, -4],
                [-5, -5, -5, 6, -4],
                [-1, -1, -4, -4, -9]]

a = dynprog(sigma, score_matrix, "ACCDCA", "ACCDCA", show_alignment=True)
print(a)

