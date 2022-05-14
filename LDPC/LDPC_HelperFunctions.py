import numpy.random
from sage.all import *
from numpy.random import default_rng
from numpy.random import uniform
from numpy.random import normal


# 2 is representing the erasure symbol
# returns the modulation of the codeword with added noise
def channel_noise(s, channel, p):
    F = RealField(7)
    if channel == 'BSC':
        #noise = vector(F, [1 if x <= p else 0 for x in list(uniform(0, 1, size=len(s)))])
        noisepos = sample(range(0, len(s)), floor(len(s) * p))
        noise = vector(F, [1 if a in noisepos else 0 for a in range(len(s))])
        r = vector(F, list(map(lambda y: (2 * y) - 1, (vector(F, s)+noise) % 2)))
    elif channel == 'AWGN':
        #noise = vector(F, list(default_rng().normal(0, p, len(s))))
        inoise = vector(RealField(10), (numpy.random.normal(0, p, size=len(s))))
        #r = 2*vector(F, s) - vector(F, [1]*len(s)) + noise
        r = (2*s)-vector([1]*len(s))+inoise
        #ri = vector(GF(2), [1 if a > 0 else 0 for a in r])

    else: # channel == 'BEC'
        noise = vector(F, [1 if x <= p else 0 for x in list(uniform(0, 1, size=len(s)))])
        #noisepos = sample(range(0, len(s)), floor(len(s)*p))
        #noise = vector(F, [1 if a in noisepos else 0 for a in range(len(s))])
        #noise = list(map(lambda lis, i: lis[i] = 1, [0]*len(s)) list(random.sample(range(0, len(s)), floor(len(s)*p))))
        s_mod = vector(F, list(map(lambda y: (2 * y) - 1, vector(F, s))))
        r = vector(F, [2 if noise[i] == 1 else s_mod[i] for i, e in enumerate(s_mod)])
    return r


def P(perm, z):
    if perm == -1:
        return zero_matrix(GF(2), z, z)
    if perm == 0:
        return identity_matrix(GF(2), z)
    else:
        res_matrix = []
        for i in range(z):
            row = [0]*z
            row[(perm+i)%z]= 1
            res_matrix.append(row)
        return matrix(GF(2), res_matrix)


def Protograph(base_matrix, z):
    protograph = Matrix(GF(2), len(base_matrix.rows())*z, len(base_matrix.columns())*z)
    for i_row, row in enumerate(base_matrix):
        for i_col, col in enumerate(row):
            protograph.set_block(i_row * z, i_col * z, P(row[i_col], z))
    return protograph


def get_base_matrix(bg, ils, zc):
    matrix = []
    f = open("base_matrices\\" + f"NR_{bg}_{ils}_{zc}.txt", "r")
    txt_m = f.read().split('\n')
    for row in txt_m:
        row = list(map(int, [c for c in row.split(' ') if c!= '']))
        matrix.append(row)
    if matrix[-1] == []:
        matrix.pop(-1)          # matrix contains an empty list at the end
    return Matrix(matrix)


def print_H_square(H, i_row, i_col, zc):
    print(H.matrix_from_rows_and_columns(list(range(i_row*zc, (i_row+1)*zc)), list(range(i_col*zc, (i_col+1)*zc))))
