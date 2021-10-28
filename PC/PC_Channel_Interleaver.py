from sage.all import *


def r_fun(i, T):
    return (i*(2*T-i+1))/2


def main_channel_interleaver(e, E, channel):
    I_BIL = channel in ['PDCCH', 'PBCH']
    if I_BIL:
        f = []
        T = ceil( (sqrt(8*E+1)-1)/2 )
        V = matrix(T)
        for i in range(T):
            for j in range(T):
                if i+j >= T or r_fun(i, T)+j >= E:
                    V[i, j] = 3  # Matrix can't store None values
                else:
                    V[i, j] = e[r_fun(i, T)+j]
        k = 0
        for j in range(T):
            for i in range(T-j):
                if V[i, j] in [0, 1]:
                    f.append(V[i, j])
                    k += 1
        return f, I_BIL
    return e[:], I_BIL


def inv_channel_interleaver(f, E, I_BIL):
    if not I_BIL:
        return f
    else:
        T = ceil((sqrt(8*E+1)-1)/2)
        V = Matrix(T)
        k = 0
        for i in range(T):
            for j in range(T):
                if j+i >= T or r_fun(j, T)+i >= E:
                    V[j, i] = 3
                else:
                    V[j, i] = f[k]
                    k += 1
        e = []
        for i in range(T):
            for j in range(T):
                if V[i, j] == 3:
                    break
                else:
                    e.append(V[i, j])
        return e
















