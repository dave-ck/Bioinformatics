from submit_main import dynprog, score, dynproglin, heuralign


def score_alignment(alphabet, scoring_matrix, seq_s, seq_t, alignment_s, alignment_t):
    if not alignment_s and not alignment_t:  # todo. catch case where len(ali_s)==1
        return 0
    total_score = 0
    for i in range(alignment_s[0], alignment_s[-1] + 1):
        if i in alignment_s:
            index = alignment_s.index(i)
            j = alignment_t[index]
            subst = score(alphabet, scoring_matrix, seq_s[i], seq_t[j])
            total_score += subst
        else:
            delete = score(alphabet, scoring_matrix, seq_s[i], "_")
            total_score += delete
    for j in range(alignment_t[0], alignment_t[-1] + 1):
        if j in alignment_t:
            index = alignment_t.index(j)
            i = alignment_s[index]
            pass  # already accounted for in loop above; don't want to double up
        else:
            delete = score(alphabet, scoring_matrix, seq_t[j], "_")
            total_score += delete
    return total_score


test_cases = {"test.py": {"test_params": ("ABC",
                                          [[1, -1, -2, -1],
                                           [-1, 2, -4, -1],
                                           [-2, -4, 3, -2],
                                           [-1, -1, -2, 0]],
                                          "AABBAACA", "CBACCCBA"),
                          "correct_indices": ([3, 5, 6], [1, 2, 3]),
                          "correct_score": 5},
              "interspersed": {"test_params": ("ABCD",
                                               [[1, -5, -5, -5, -1],
                                                [-5, 1, -5, -5, -1],
                                                [-5, -5, 5, -5, -4],
                                                [-5, -5, -5, 6, -4],
                                                [-1, -1, -4, -4, -9]],
                                               "AAAAACCDDCCDDAAAAACC",
                                               "CCAAADDAAAACCAAADDCCAAAA"),
                               "correct_indices": ([5, 6, 7, 8, 9, 10, 11, 12, 18, 19],
                                                   [0, 1, 5, 6, 11, 12, 16, 17, 18, 19]),
                               "correct_score": 39},
              "single matches": {"test_params": ("ABCD",
                                                 [[1, -5, -5, -5, -1],
                                                  [-5, 1, -5, -5, -1],
                                                  [-5, -5, 5, -5, -4],
                                                  [-5, -5, -5, 6, -4],
                                                  [-1, -1, -4, -4, -9]],
                                                 "AACAAADAAAACAADAADAAA",
                                                 "CDCDDD"),
                                 "correct_indices": ([2, 6, 11, 14, 17], [0, 1, 2, 3, 4]),
                                 "correct_score": 17},
              "lonk": {"test_params": ("ABCD",
                                       [[1, -5, -5, -5, -1],
                                        [-5, 1, -5, -5, -1],
                                        [-5, -5, 5, -5, -4],
                                        [-5, -5, -5, 6, -4],
                                        [-1, -1, -4, -4, -9]],
                                       "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
                                       "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"),
                       "correct_indices": (
                           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55,
                            56,
                            57, 58],
                           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66,
                            67,
                            68, 69]),
                       "correct_score": 81},
              "custom": {"test_params": ("ABC",
                                         [[1, -1, -2, -1],
                                          [-1, 2, -4, -1],
                                          [-2, -4, 3, -2],
                                          [-1, -1, -2, 0]],
                                         "AABBAACA", "CBACCCBA"),
                         "correct_indices": ([3, 5, 6], [1, 2, 3]),
                         "correct_score": 5},
              }

for test_choice in ["interspersed", "lonk"]:
    a = heuralign(*test_cases[test_choice]["test_params"], False)
    print(a)
    if a[0] == test_cases[test_choice]["correct_score"]:
        print("Score correctly found:", a[0])
    else:
        print("Correct score:", test_cases[test_choice]["correct_score"])
        print("Score found:", a[0])
    if (a[1], a[2]) == test_cases[test_choice]["correct_indices"]:
        print("Indices correctly found: ", a[1], a[2])
    else:
        print("Correct indices:", *test_cases[test_choice]["correct_indices"])
        print("Indices found: ", a[1], a[2])
        print("Indices found give a score of:", score_alignment(*test_cases[test_choice]["test_params"], a[1], a[2]))
