import numpy as np


# IMPORTANT: refactor and verify function signatures


def dynprog_np(alphabet, scoring_matrix, seq_s, seq_t, debug=False):
    scoring_matrix = np.array(scoring_matrix)
    m, n = len(seq_s), len(seq_t)
    backtrack = np.empty([m + 1, n + 1], dtype=str)  # points U(p), L(eft), D(iagonal), or E(nd)
    backtrack[:, :] = "E"
    matrix_v = np.zeros([m + 1, n + 1], dtype=int)  # because doing local alignment, no negative values
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            best = 0
            reset = True
            from_diagonal = matrix_v[i - 1, j - 1] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1])
            from_up = matrix_v[i, j - 1] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_left = matrix_v[i - 1, j] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
            if from_left >= matrix_v[i - 1, j]:
                print(from_left)
                print(matrix_v[i - 1, j])
                print("Wat.")
            if best <= from_diagonal:
                best = from_diagonal
                backtrack[i, j] = "D"
                matrix_v[i, j] = best
                reset = False
            if best <= from_left:
                best = from_left
                backtrack[i, j] = "U"
                matrix_v[i, j] = best
                reset = False
            if best <= from_up:
                best = from_up
                backtrack[i, j] = "L"
                matrix_v[i, j] = best
                reset = False
            if reset:
                matrix_v[i, j] = 0  # explicit but redundant - overwrite 0 from np.zeros init
                backtrack[i, j] = "E"
    alignment_s, alignment_t = [], []
    show_s, show_t = "", ""
    best_score = np.amax(matrix_v)
    print("Best score (in-function, at compute):", best_score)
    best_mat = np.where(matrix_v == best_score)
    i, j = best_mat[0][0], best_mat[1][0]
    pointer = backtrack[i, j]
    print("Max value is {}, found at {}, {}".format(matrix_v[i, j], i, j))
    print(matrix_v)
    print(backtrack)
    turn = 0
    while pointer != "E":
        print("Step {}, pointer {}".format(turn, pointer))
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
    display_alignment(show_s, show_t)
    print("Best score (in-function, at return):", best_score)
    return [show_s, show_t, best_score, alignment_s, alignment_t]


def dynprog(alphabet, scoring_matrix, seq_s, seq_t):
    # m rows, n columns - bottom corner is mat[n+1][m+1]
    m, n = len(seq_s), len(seq_t)
    matrix_v = init_v(m, n)
    backtrack = init_backtrack(m, n)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            matrix_v[j][i] = max(matrix_v[j - 1][i - 1]
                                 + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1]),
                                 matrix_v[j - 1][i]
                                 + score(alphabet, scoring_matrix, "_", seq_t[j - 1]),
                                 matrix_v[j][i - 1]
                                 + score(alphabet, scoring_matrix, seq_s[i - 1], "_"))

            if matrix_v[j][i] == matrix_v[j - 1][i - 1] + score(alphabet, scoring_matrix, seq_s[i - 1], seq_t[j - 1]):
                backtrack[j][i] = "D"
            elif matrix_v[j][i] == matrix_v[j - 1][i] + score(alphabet, scoring_matrix, "_", seq_t[j - 1]):
                backtrack[j][i] = "U"
            elif matrix_v[j][i] == matrix_v[j][i - 1] + score(alphabet, scoring_matrix, seq_s[i - 1], "_"):
                backtrack[j][i] = "L"
            if matrix_v[j][i] < 0:
                matrix_v[j][i] = 0
                backtrack[j][i] = "E"
    i, j = m, n
    aligned_s = ""
    aligned_t = ""
    alignment_s = []
    alignment_t = []
    pointer = backtrack[j][i]
    while pointer != "E":
        if pointer == "D":
            j -= 1
            i -= 1
            alignment_s.insert(0, i)
            alignment_t.insert(0, j)
            aligned_s = seq_s[i] + aligned_s
            aligned_t = seq_t[j] + aligned_t
            pointer = backtrack[j][i]
        elif pointer == "L":
            i -= 1
            aligned_s = seq_s[i] + aligned_s
            aligned_t = " " + aligned_t
            pointer = backtrack[j][i]
        elif pointer == "U":
            j -= 1
            aligned_t = seq_t[j] + aligned_t
            aligned_s = " " + aligned_s
            pointer = backtrack[j][i]
    return [aligned_s, aligned_t, matrix_v[n][m], alignment_s, alignment_t]


def init_backtrack(length1, length2):
    mat = [["" for i in range(length1 + 1)] for j in range(length2 + 1)]
    for i in range(1, length1 + 1):
        mat[0][i] = "L"
    for j in range(1, length2 + 1):
        mat[j][0] = "U"
    mat[0][0] = "E"
    return mat


def init_v(length1, length2):
    mat = [[0 for i in range(length1 + 1)] for j in range(length2 + 1)]
    for i in range(1, length1 + 1):
        mat[0][i] = mat[0][i - 1] - 2
    for j in range(1, length2 + 1):
        mat[j][0] = mat[j - 1][0] - 2
    return mat


# one-line function, leave to keep clean (?)
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

a = dynprog_np(sigma, score_matrix, "AAAAACCDDCCDDAAAAACC", "CCAAADDAAAACCAAADDCCAAAA")
print("Score:   ", a[2])
print("Indices: ", a[3], a[4])
