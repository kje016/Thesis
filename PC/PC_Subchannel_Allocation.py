# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
from collections import Counter

import HelperFunctions
import PC_Rate_Matching as RM
P = [0, 1, 2, 4, 3, 5, 6, 7, 8, 16, 9, 17, 10, 18, 11, 19, 12, 20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31]

reliability_sequence = HelperFunctions.get_realiability_sequence()

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
        if elem < N:
            Q_N_0.append(elem)
        if len(Q_N_0) >= N:
            break
    return Q_N_0


def freeze(N, K, E, npc):
    """QNI = []
    if E < N:
        if K/E <= 7/16: # puncturing
            for n in range(N-E):
                QNI.append(J_fun(n, N))
            if E >= (3*N)/4:
                QNI.extend([a for a in range(ceil((3*N)/4 - E/2))])
            else:
                QNI.extend([a for a in range(ceil((9*N)/16 - E/4))])
        else: # shortening
            for n in range(E, N):
                QNI.append(J_fun(n, N))"""
    QN0 = get_Q_N0(N)
    matching_scheme = RM.matching_selection(E, N, K)
    MS = RM.get_rm_set(U=N-E, matching_scheme=matching_scheme, QN0=QN0)
    QNI = set(list(set(QN0)-MS)[-(K+npc):])
    QNF = set(QN0) - QNI
    return QNF, QNI, MS, matching_scheme


def main(N, c_ap, A, E, channel):
    npc, n_wm_pc = get_n_pc_bits(channel, A, E)
    QNF, QNI, MS, matching_scheme = freeze(N, len(c_ap), E, npc)
    u = []
    c_get = 0
    for i in range(N):
        if i in QNF:
            u.append(0)
        else:
            u.append(c_ap[c_get])
            c_get += 1
    return u, npc, n_wm_pc, QNF, QNI, MS, matching_scheme

"""
def main(N, c_ap, A, E, channel):
    K = len(c_ap)
    npc, n_wm_pc = get_n_pc_bits(channel, A, E)
    QNF, QNI = freeze(N, K, E, npc)
    u = []
    c_get = 0
    for i in range(N):
        if i in QNF:
            u.append(0)
        else:
            u.append(c_ap[c_get])
            c_get += 1
    return u, npc, n_wm_pc, QNF, QNI
"""