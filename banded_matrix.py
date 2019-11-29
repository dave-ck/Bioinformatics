import time
import numpy as np


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


# TEST CODE BELOW - COMPARISON WITH "NAIVE" QUADRATIC-SPACE IMPLEMENTATION
k = 3
len_t = 15  # j indexed
len_s = 10  # i indexed
diff = 0

m_good = MatrixBand(diff + k, diff - k, len_t, dtype='int', default=0)
m_good.set(0, 0, 3)
m_good.set(1, 0, 2)  # good
m_good.set(3, 1, 3)  # good
m_good.set(0, 1, 5)  # good
m_good.max_value_coords()
m_good.print()
