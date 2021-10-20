from sage.all import *

P = [0, 1, 2, 4, 3, 5, 6, 7, 8, 16, 9, 17, 10, 18, 11, 19, 12, 20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31]


def main_sub_block_interleaver(d, N):
    y = []
    for n in range(N):
        i = floor((32*n)/N)
        jn = P[i]*(N//32) + Integer(mod(n, N//32))
        y.append(d[jn])
    return y

# TODO: må sikkert ha en for-løkke for hele blokken