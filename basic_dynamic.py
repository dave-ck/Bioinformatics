import numpy

# IMPORTANT: refactor and verify function signatures

sigma = "ABC"
matrix = [[1, -1, -2, -1], [-1, 2, -4. - 1], [-2, -4, 3, -2], [-1, -1, -2, 0]]
seqA = "ABAACACA"
seqB = "CABBCACA"


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


def pp(o):  # pretty print function, prints each line of 2D array separately for easier reading
    for i in o:
        print(i)


def display_alignment(alignment):
    string1 = alignment[0]
    string2 = alignment[1]
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


# display_alignment(dynprog(sigma, matrix, seqA, seqB))
a = dynprog("ABCD",[
[ 1,-5,-5,-5,-1],
[-5, 1,-5,-5,-1],
[-5,-5, 5,-5,-4],
[-5,-5,-5, 6,-4],
[-1,-1,-4,-4,-9]], "AAAAACCDDCCDDAAAAACC", "CCAAADDAAAACCAAADDCCAAAA")
display_alignment(a)
print("Score:   ", a[2])
print("Indices: ", a[3], a[4])

print("\n\n")
a = dynprog("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
print("Score:   ", a[2])
print("Indices: ", a[3],a[4])
display_alignment(a)
