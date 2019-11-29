import time
import numpy as np


class MatrixBand:
    def __init__(self, max_diff, min_diff, len_t, dtype, default):
        self.max_diff = max_diff
        self.min_diff = min_diff
        self.array = np.empty([max_diff - min_diff, len_t], dtype=dtype)
        self.array[:] = default

    def readdress(self, i, j):
        # j remains the same
        if i < 0 or j < 0 or i - j < self.min_diff or i - j > self.max_diff:
            raise ValueError("Coordinates which were either negative or outside the band were passed.")
        i = i - j - self.min_diff
        return i, j

    def get(self, i, j):
        return self.array[self.readdress(i, j)]

    def set(self, i, j, value):
        self.array[self.readdress(i, j)] = value

    def print(self):
        print("Array is now")
        print(self.array)



"""
# TEST CODE BELOW - COMPARISON WITH "NAIVE" QUADRATIC-SPACE IMPLEMENTATION
k = 3
len_t = 8  # j indexed
len_s = 10  # i indexed
diff = 2

m_good = MatrixBand(diff + k + 1, diff - k, len_t, dtype='str', default='-')
m_good.set(0, 0, 'C')
m_good.set(1, 0, 'D')    # good
m_good.set(3, 1, 'A')    # good
m_good.set(8, 4, 'A')    # good
m_good.set(0, 1, 'U')    # good
m_good.print()

m_bad = np.array([[i - j - diff for i in range(len_t)] for j in range(len_s)])
m_bad = abs(m_bad)
m_bad[m_bad > abs(k)] = -1
m_bad[m_bad > 0] = 1

m_bad[0, 0] = 15
m_bad[1, 1] = 15
m_bad[6, 3] = 15
print(m_bad)

start = time.time()
m_good = MatrixBand(diff + k, diff - k -1 , len_t)
print("Generated linear-space band in: {}".format(time.time()-start))

start = time.time()
m_bad = np.zeros([len_t, len_s], dtype='int8')
print("Generated quadratic-space matrix in: {}".format(time.time()-start))


"""
