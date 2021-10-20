# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
P = [0, 1, 2, 4, 3, 5, 6, 7, 8, 16, 9, 17, 10, 18, 11, 19, 12, 20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31]

with open('Reliability_Sequence.txt') as f:
    reliability_sequence = f.read()
reliability_sequence = reliability_sequence.replace(' ', '')
reliability_sequence = reliability_sequence.split(",")
reliability_sequence = [int(x) for x in reliability_sequence]


def get_n_pc_bits(channel, A, E):
    n_pc = 0
    n_wm_pc = 0
    if channel in ["PUCCH", "PUSCH"] and 12 <= A <=19:
        n_pc = 3
        if E-A > 175:
            n_wm_pc = 1
    return n_pc, n_wm_pc


def J_fun(elem, N):
    i = floor((32*elem)/N)
    result = P[i] * (N//32) + Integer(mod(elem, N/32))
    return result



def get_Q_N0(N):
    Q_N_0 = []
    for elem in reliability_sequence:
        if elem-1 <= N:
            Q_N_0.append(elem-1)
        if len(Q_N_0) >= N:
            break
    return Q_N_0


def freeze(N, K, E, npc):
    Q_N_F = []
    if E < N:
        if K/E <= 7/16: # puncturing
            for n in range(N-E):
                Q_N_F = Q_N_F.union({J_fun(n, N)})
            if E >= (3*N)/4:
                Q_N_F.append([a for a in range(ceil((3*N)/4 - E/2))])
            else:
                Q_N_F.append([a for a in range(ceil((9*N)/16 - E/4))])
        else: # shortening
            for n in range(E, N):
                Q_N_F.append(J_fun(n, N))

    Q_N_0 = get_Q_N0(N)
    Q_N_I = set(Q_N_0[:K+npc])
    Q_N_F = set(Q_N_0)-Q_N_I

    return Q_N_I, Q_N_F


def main(N, c_ap, A, E, channel):
    K = len(c_ap)
    npc, n_wm_pc = get_n_pc_bits(channel, A, E)
    QNI, QNF = freeze(N, K, E, npc)
    u = []
    c_get = 0
    for i in range(N):
        if i in QNF:
            u.append(0)
        else:
            u.append(c_ap[c_get])
            c_get += 1
    return u, npc, n_wm_pc

