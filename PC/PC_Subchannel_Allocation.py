# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
from collections import Counter
import HelperFunctions
import PC_Rate_Matching as RM

reliability_sequence = HelperFunctions.get_realiability_sequence()


def get_n_pc_bits(channel, A, E):
    n_pc, n_wm_pc = 0, 0
    if channel in ["PUCCH", "PUSCH"] and 12 <= A <=19:
        n_pc = 3
        if E-A > 175:
            n_wm_pc = 1
    return n_pc, n_wm_pc


def get_Q_N0(N):
    Q_N_0 = []
    for elem in reliability_sequence:
        if elem < N:
            Q_N_0.append(elem)
        if len(Q_N_0) >= N:
            break
    return Q_N_0


def freeze(N, K, E, npc):
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
