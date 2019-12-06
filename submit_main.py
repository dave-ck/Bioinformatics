import numpy as np


# TODO: IMPORTANT: refactor and verify function signatures

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
            from_up = matrix_v[i, j - 1] + score(alphabet, scoring_matrix, "_", seq_t[j - 1])
            from_left = matrix_v[i - 1, j] + score(alphabet, scoring_matrix, seq_s[i - 1], "_")
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
            show_t = seq_t[j] + show_t
            show_s = " " + show_s
            pointer = backtrack[i, j]
        elif pointer == "U":
            i -= 1
            show_s = seq_s[i] + show_s
            show_t = " " + show_t
            pointer = backtrack[i, j]
    if show_alignment:
        display_alignment(show_s, show_t)
    return [best_score, alignment_s, alignment_t]


# one-line function, keep for cleaner code
def score(alphabet, scoring_matrix, char1, char2):
    return scoring_matrix[(char1 == "_") * -1 or alphabet.index(char1)][(char2 == "_") * -1 or alphabet.index(char2)]


# todo: remove (along with references) from submission .py
def display_alignment(show_s, show_t):
    string3 = ''
    for i in range(min(len(show_s), len(show_t))):
        if show_s[i] == show_t[i]:
            string3 = string3 + "|"
        else:
            string3 = string3 + " "
    print('Alignment ')
    print('String s: ' + show_s)
    print('          ' + string3)
    print('String t: ' + show_t + '\n\n')

def dynproglin(alphabet, scoring_matrix, seq_s, seq_t, show_alignment=False):
    score, start_ptr, end_ptr = get_local_score_and_endpoints(alphabet, scoring_matrix, seq_s, seq_t)
    start_s, start_t = start_ptr
    end_s, end_t = end_ptr
    indices_s, indices_t = hirschberg(alphabet, scoring_matrix, seq_s[start_s:end_s], seq_t[start_t:end_t])
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
    if len(seq_s) == 1 or len(seq_t) == 1:
        alignment_s, alignment_t = needleman_wunsch_char_v_seq(alphabet, scoring_matrix, seq_s, seq_t)
    else:
        len_t = len(seq_t)
        t_mid = len_t // 2
        score_left = get_global_score(alphabet, scoring_matrix, seq_s, seq_t[:t_mid])
        rev_s = seq_s[::-1]
        seq_t_pt2 = seq_t[t_mid:]
        rev_t_pt2 = seq_t_pt2[::-1]
        score_right = get_global_score(alphabet, scoring_matrix, rev_s, rev_t_pt2)
        rev_right = list(reversed(score_right))
        sum_scores = score_left + rev_right
        s_mid = np.where(sum_scores == np.amax(sum_scores))[0][0]
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

def heuralign(alphabet, scoring_matrix, seq_s, seq_t, show_alignment=False):
    ktup = 2
    band_radius = 15
    min_score_to_extend = -3
    min_seeds_to_extend = 2
    total_diags_to_extend = 10
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
    max_score_for_diff = {}
    all_scores = []
    for j in range(len_t):
        current_subseq = seq_t[j:j + ktup]
        if current_subseq in index_table:
            for i in index_table[current_subseq]:
                if i - j in ij_diff_dict:
                    ij_diff_dict[i - j].append(i)  # instead of storing {diff:(i,j)} store {diff:i} and find j=i+diff
                else:
                    ij_diff_dict[i - j] = [i]
    for diff in ij_diff_dict:
        if len(ij_diff_dict[diff]) > min_seeds_to_extend:
            max_score_for_diff[diff] = 0
            seed_list = ij_diff_dict[diff]
            # display_diagonal(seq_s, seq_t, diff)
            for seed_i in seed_list:
                # j = diff - i
                start_i, end_i = seed_i, seed_i + ktup  # non-inclusive of end index
                start_j, end_j = start_i - diff, start_i + ktup - diff  # tested j-calculation [gud]
                best_start_i, best_end_i = start_i, end_i
                best_start_j, best_end_j = start_j, end_j
                s_sub = seq_s[start_i:end_i]
                t_sub = seq_t[start_j:end_j]
                current_score = score_seqs(alphabet, scoring_matrix, s_sub, t_sub)  # maintain current score
                best_score = current_score
                updated = True
                while updated:
                    updated = False
                    while current_score > min_score_to_extend:
                        # extend left by decreasing start_i and end_j
                        if start_i < 1 or start_j < 1:
                            break  # cannot decrease index beyond zero - done with this loop
                        start_i -= 1
                        start_j -= 1
                        current_score += score(alphabet, scoring_matrix, seq_s[start_i], seq_t[start_j])
                        if start_i in seed_list:  # ensure only one pass
                            seed_list.remove(start_i)
                        if current_score > best_score:
                            updated = True
                            best_score = current_score
                            best_start_i, best_start_j = start_i, start_j
                    # revert to best start found so far
                    start_i, start_j = best_start_i, best_start_j
                    current_score = best_score
                    while current_score > min_score_to_extend:
                        # extend right by increasing start_i and end_j
                        if end_i >= len_s or end_j >= len_t:
                            break  # cannot increase index beyond sequence size
                        end_i += 1
                        end_j += 1
                        if end_i in seed_list:  # ensure only one pass
                            seed_list.remove(end_i)
                        current_score += score(alphabet, scoring_matrix, seq_s[end_i - 1], seq_t[end_j - 1])
                        if current_score > best_score:
                            updated = True
                            best_score = current_score
                            best_end_i, best_end_j = end_i, end_j
                    end_i, end_j = best_end_i, best_end_j
                    current_score = best_score
                    # reset/backtrack to highest score yet
                max_score_for_diff[diff] = max(max_score_for_diff[diff], best_score)
                all_scores.append(best_score)
    if not all_scores:
        print("Could not find {} seeds on the same (i-j) diagonal with ktup = {}; trying banded DP on i-j=0".format(min_seeds_to_extend, ktup))
        return banded_dynprog(alphabet, scoring_matrix, seq_s, seq_t, band_radius*2, band_radius*2)
    cutoff = sorted(all_scores)[-1 * min(len(all_scores), total_diags_to_extend)]
    good_diffs = []
    for diff in max_score_for_diff:
        if max_score_for_diff[diff] > cutoff:
            good_diffs.append(diff)
    # run banded DP on everything in good_diffs
    best = [0, [], []]
    diff_bandwidth = {i: (i - band_radius, i + band_radius) for i in good_diffs}
    for midpoint in range(-1*len_t, len_s):  # combine diagonals
        for diag_1 in good_diffs:
            if midpoint in range(*diff_bandwidth[diag_1]):
                for diag_2 in good_diffs:
                    if midpoint in range(*diff_bandwidth[diag_2]) and diag_2 != diag_1:
                        good_diffs.remove(diag_2)
                        diff_bandwidth[diag_1] = (min(diff_bandwidth[diag_1][0], diff_bandwidth[diag_2][0]),
                                                  max(diff_bandwidth[diag_1][1], diff_bandwidth[diag_2][1]))
                        del diff_bandwidth[diag_2]
    for diff in good_diffs:
        min_diff, max_diff = diff_bandwidth[diff]
        answer = banded_dynprog(alphabet, scoring_matrix, seq_s, seq_t, max_diff, min_diff, show_alignment)
        if answer[0] > best[0]:
            best = answer
    return best

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


def score_seqs(alphabet, scoring_matrix, seq1, seq2):
    return sum(score(alphabet, scoring_matrix, char1, char2) for (char1, char2) in zip(seq1, seq2))


class MatrixBand:
    def __init__(self, max_diff, min_diff, len_t, dtype, default):
        self.max_diff = max_diff
        self.min_diff = min_diff
        self.array = np.empty([max_diff - min_diff + 1, len_t], dtype=dtype)
        self.array[:] = default

    def readdress(self, i, j):
        # j remains the same
        if i < 0 or j < 0 or i - j < self.min_diff or i - j > self.max_diff:
            raise ValueError("Coordinates which were either negative or outside the band were passed.")
        i = i - j - self.min_diff
        return i, j

    def revert_address(self, i, j):  # inverse of readdress function
        return j + i + self.min_diff, j

    def get(self, i, j):
        return self.array[self.readdress(i, j)]

    def set(self, i, j, value):
        self.array[self.readdress(i, j)] = value

    def print(self):
        print("Array is now")
        print(self.array)

    def max_value_coords(self):
        best_score = np.amax(self.array)
        best_mat = np.where(self.array == best_score)
        i, j = best_mat[0][0], best_mat[1][0]
        coords_out = self.revert_address(i, j)
        if self.array[i, j] == 0:
            coords_out = (0, 0)
        return coords_out