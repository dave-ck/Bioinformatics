import numpy as np
import json  # todo: remove import - used only for debug prints


def heuralign(alphabet, scoring_matrix, seq_s, seq_t, ktup=2, show_alignment=False):
    if len(seq_s) < len(seq_t):
        seq_t, seq_s = seq_s, seq_t
    # initialize index table as dictionary of every sequence of letters of length ktup appearing in s
    scoring_matrix = np.array(scoring_matrix)
    len_s, len_t = len(seq_s), len(seq_t)
    index_table = {}
    for i in range(len_s - ktup):
        x = seq_s[i:i + ktup]
        if x not in index_table:
            index_table.update({x: [i]})  # initialize to empty list
        else:
            index_table[x].append(i)
    # sort according to distance i-j
    ij_diff_dict = {}
    for j in range(len_t):
        current_subseq = seq_t[j:j + ktup]
        if current_subseq in index_table:
            for i in index_table[current_subseq]:
                if i - j in ij_diff_dict:
                    ij_diff_dict[i - j].append(i)  # instead of storing {diff:(i,j)} store {diff:i} and find j=i+diff
                else:
                    ij_diff_dict[i - j] = [i]
    # tag diagonals (frequence analysis for each? naive count of matches along diagonal?)

    # for high i-j offset frequency (diagonals with many pieces) combine pieces into regions by extending pieces
    # greedily in a straight line

    # identify those diagonals which contain the most matches (?)

    # use banded DP on these diagonals (i.e. with width 12 - why??)
    print(json.dumps(index_table, indent=2))
    print(json.dumps(ij_diff_dict, indent=2))


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


sigma = "ABCD"
score_matrix = [[1, -5, -5, -5, -1],
                [-5, 1, -5, -5, -1],
                [-5, -5, 5, -5, -4],
                [-5, -5, -5, 6, -4],
                [-1, -1, -4, -4, -9]]

s = "CCAAADDAAAACCAAADDCCAAAA"
t = "AAAAACCDDCCDDAAAAACC"

heuralign(sigma, score_matrix, s, t, 3)
