import banded_dp
import numpy as np
import json  # todo: remove import - used only for debug prints


def heuralign(alphabet, scoring_matrix, seq_s, seq_t, show_alignment=False):
    ktup = 2
    band_radius = 15
    min_score_to_extend = -3
    min_seeds_to_extend = 2
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
    print("Diagonal seeds dict")
    print(json.dumps(ij_diff_dict, indent=2))
    # select and print maximum diagonals
    for diff in ij_diff_dict:
        if len(ij_diff_dict[diff]) > min_seeds_to_extend:
            max_score_for_diff[diff] = 0
            seed_list = ij_diff_dict[diff]
            # display_diagonal(seq_s, seq_t, diff)
            for seed_i in seed_list:
                # todo: properly cull seeds from seed_list during search
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
                # print("Best alignment from seed_i = {} and diff = {}: score = {}:".format(seed_i, diff, best_score))
                # print("Best start i = {}; Best start j = {}".format(best_start_i, best_start_j))
                # print("Best end i = {}; Best end j = {}".format(best_end_i, best_end_j))
                # s_sub = seq_s[best_start_i:best_end_i]
                # t_sub = seq_t[best_start_j:best_end_j]
                # print("Score calculated from scratch on:\n{}\n{}\n = {}:".format(s_sub, t_sub,
                #                                                                  score_seqs(alphabet,
                #                                                                             scoring_matrix,
                #                                                                             t_sub, s_sub)))
                # display_alignment(s_sub, t_sub)
    print(max_score_for_diff)
    if not all_scores:
        return "Could not find {} seeds on the same (i-j) diagonal with ktup = {}".format(min_seeds_to_extend, ktup)
    cutoff = sorted(all_scores)[-1 * min(len(all_scores), 10)]  # take 10th biggest as cutoff, if it exists
    good_diffs = []
    for diff in max_score_for_diff:
        if max_score_for_diff[diff] > cutoff:
            good_diffs.append(diff)
    # run banded DP on everything in good_diffs
    best = [0, [], []]
    for diff in good_diffs:
        max_diff = diff + band_radius
        min_diff = diff - band_radius
        answer = banded_dp.banded_dynprog(alphabet, scoring_matrix, seq_s, seq_t, max_diff, min_diff, show_alignment)
        if answer[0] > best[0]:
            best = answer
    return best


# have obtained all diagonals
# select highest-scoring ones
# identify ones whose bands would overlap - these get wider bands
# run banded DP on their bands with base width 12


# one-line function, keep for cleaner code
def score(alphabet, scoring_matrix, char1, char2):
    return scoring_matrix[(char1 == "_") * -1 or alphabet.index(char1)][(char2 == "_") * -1 or alphabet.index(char2)]


def score_seqs(alphabet, scoring_matrix, seq1, seq2):
    return sum(score(alphabet, scoring_matrix, char1, char2) for (char1, char2) in zip(seq1, seq2))


def display_diagonal(string1, string2, diff):
    print("With ijdiff == {}:".format(diff))
    if diff < 0:
        display_alignment(-1 * diff * "-" + string1, string2)
    else:
        display_alignment(string1, diff * "-" + string2)


def display_alignment(string1, string2):
    string3 = ""
    for i in range(min(len(string1), len(string2))):
        if string1[i] == string2[i]:
            string3 = string3 + "|"
        else:
            string3 = string3 + " "
    print("String1: " + string1)
    print("         " + string3)
    print("String2: " + string2 + "\n\n")


sigma = "ABCD"
score_matrix = [[1, -5, -5, -5, -1],
                [-5, 1, -5, -5, -1],
                [-5, -5, 5, -5, -4],
                [-5, -5, -5, 6, -4],
                [-1, -1, -4, -4, -9]]

s = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
t = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

print(heuralign(sigma, score_matrix, s, t, False))


print("Hello, world!")
