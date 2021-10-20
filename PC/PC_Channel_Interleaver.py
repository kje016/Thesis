from sage.all import *


def main_channel_interleaver(e, E, channel):
    I_BIL = channel in ['PDCCH', 'PBCH']
    if I_BIL:
        f = []
        T = ceil( (sqrt(8*E+1)-1)/2 )
        V = matrix(T)
        k = 0
        for i in range(T):
            for j in range(T-i):
                if k < E:
                    V[i, j] = e[k]
                else:
                    V[i, j] = 3     # Matrix can't store None values
                k += 1

        k = 0
        for j in range(T):
            for i in range(T-j):
                if V[i, j] in [0, 1]:
                    f.append(V[i, j])
                    k += 1
        return f
    return e[:]






