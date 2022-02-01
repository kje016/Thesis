# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions
import PC_Rate_Matching as RM

reliability_sequence = HelperFunctions.get_realiability_sequence()


def get_n_pc_bits(A, E, I_IL):
    n_pc, n_wm_pc = 0, 0
    if I_IL == 0 and 12 <= A <= 19:
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


def freeze(N, K, E, npc, R):
    QN0 = get_Q_N0(N)
    matching_scheme = RM.matching_selection(E=E, N=N, K=K)
    # tess = RM.get_rm_set(8, 'shortening', [0,1,2,4,3,5,6,7])
    MS = RM.get_rm_set(U=N-E, matching_scheme=matching_scheme, QN0=QN0)
    if matching_scheme == 'repetition':
        QNI = set(QN0[-K:])
        QNF = set(QN0) - QNI
    else:
        R = [a for a in QN0 if a not in MS]
        QNF = set(R[:-K]).union(MS)
        QNI = set(QN0) - QNF
    return QNF, QNI, MS, matching_scheme


def main(N, c_ap, A, E, I_IL, R):
    npc, n_wm_pc = get_n_pc_bits(A, E, I_IL)
    npc, n_wm_pc = 0, 0     # TODO: Implement pc bits
    QNF, QNI, MS, matching_scheme = freeze(N, len(c_ap), E, npc, R)

    u, c_get = [], 0
    for i in range(N):
        if i in QNF:
            u.append(0)
        else:
            u.append(c_ap[c_get])
            c_get += 1
    return u, npc, n_wm_pc, QNF, QNI, MS, matching_scheme
