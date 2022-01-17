import random

import numpy.random
from sage.all import *
from numpy.random import default_rng
from numpy.random import uniform


# 2 is representing the erasure symbol
def channel_noise(s, channel, p):
    if channel == 'BSC':
        noise = [1 if x <= p else 0 for x in list(numpy.random.uniform(0, 1, size=len(s)))]
        r = vector(GF(2), s) + vector(GF(2), noise)

    elif channel == 'AWGN':
        noise = vector(RealField(10), list(default_rng().normal(0, p, len(s))))
        r = 2*vector(RealField(10), s) - vector(RealField(10), [1]*len(s)) + noise

    else: # channel == 'BEC'
        r = vector(GF(2), [2 if list(numpy.random.uniform(0, 1, size=len(s)))[i] <= p else s[i] for i, e in enumerate(s)])
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
    matrix.pop(-1)          # matrix contains an empty list at the end
    return Matrix(matrix)


def print_H_square(H, i_row, i_col, zc):
    print(H.matrix_from_rows_and_columns(list(range(i_row*zc, (i_row+1)*zc)), list(range(i_col*zc, (i_col+1)*zc))))
