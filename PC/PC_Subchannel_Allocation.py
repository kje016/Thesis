# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions
import PC_Rate_Matching as RM

reliability_sequence = HelperFunctions.get_realiability_sequence()


def get_n_pc_bits(A, E, I_IL):
    n_pc, n_wm_pc = 0, 0
    return n_pc, n_wm_pc # TODO: remove this line
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


def freeze(N, K, E):
    QN0 = get_Q_N0(N)
    matching_scheme = RM.matching_selection(E=E, N=N, K=K)
    MS = RM.get_rm_set(U=N-E, matching_scheme=matching_scheme, QN0=QN0)
    if len(MS) == 0:
        QNI = set(QN0[-K:])
        QNF = set(QN0)-QNI
    elif matching_scheme == 'repetition':
        QNI = set(QN0[-K:])
        QNF = set(QN0) - QNI
    else:
        R = [a for a in QN0 if a not in MS]
        QNF = set(R[:N-K-len(MS)]).union(MS)
        QNI = set(QN0) - QNF
    return QNF, QNI, MS, matching_scheme


def main(N, c_ap, K, E, I_IL, R):
    npc, n_wm_pc = get_n_pc_bits(K, E, I_IL)
    QNF, QNI, MS, matching_scheme = freeze(N, K, E)
    QNPC = QNI[-(npc-n_wm_pc):]

    u, c_get = [], 0
    for i in range(N):
        if i in QNI:
            u.append(c_ap[c_get])
            c_get += 1
        else:
            u.append(0)
    return u, npc, n_wm_pc, QNF, QNI, MS, matching_scheme


def get_n_wm_pc(GN, n_wm_pc, QNI, npc):
    if n_wm_pc <= 0:
        return None

    pc1 = QNI[-(npc - n_wm_pc):]
    Q_tilde = list(filter(lambda a: a not in pc1, QNI))
    GN_I = list(map(lambda a: a.hamming_weight(), GN.matrix_from_rows(Q_tilde)))
    return pc1 + [Q_tilde[GN_I.index(max(GN_I))]]


def calc_u(N, QNI, c_ap, QNPC):
    k, u = 0, []
    if QNPC:
        y0, y1, y2, y3, y4 = 0, 0, 0, 0, 0
        for n in range(N):
            yt = y0
            y0, y1, y2, y3, y4 = y1, y2, y3, y4, yt
            if n in QNI:
                if n in QNPC:
                    u.append(y0)
                else:
                    u.append(c_ap[k])
                    k = k+1
                    y0 = mod(y0+u[n], 2)
            else:
                u.append(0)
    else:
        for i in range(N):
            if i in QNI:
                u.append(c_ap[k])
                k = k+1
            else:
                u.append(0)
    return u
