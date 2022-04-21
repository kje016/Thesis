# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


def matching_selection(E, N, K):
    if N > E:
        if K / E <= 7 / 16:
            rate_matching = "puncturing"
        else:
            rate_matching = "shortening"
    else:
        rate_matching = "repetition"
    return rate_matching


def get_rm_set(U, matching_scheme, QN0):
    bit_len = len("{0:b}".format(len(QN0)-1))
    bm = [int(("0" * (bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1], 2) for a in list(range(len(QN0)))]
    if matching_scheme == "puncturing":
        ms = set([bm.index(a) for a in QN0[:U]])
    elif matching_scheme == "shortening":
        ms = set([bm.index(a) for a in QN0[-U:]])
    else: # matching_scheme == 'repetition'
        ms = set([QN0.index(a) for a in bm[:abs(U)]])
    return ms


def circular_buffer(y, matching_set, matching_scheme):
    if matching_scheme == 'repetition':
        return y[:] + [y[a] for a in matching_set]
    e = []
    for a in range(len(y)):
        if a not in matching_set:
            e.append(y[a])
    return e


def bec_inv_circbuf(N, ee, matching_scheme, MS):
    y, counter = [], 0
    ms_operand = oo if matching_scheme == 'shortening' else 0

    for a in range(N):
        if a in MS:
            yi = (ee[N+counter] if ee[counter] == 2 else ee[counter])*-oo if matching_scheme == 'repetition' else\
                ms_operand
            y.append(yi)
        else:
            y.append(2 if ee[counter] == 2 else ee[counter]*-oo)
            counter += 1
    return vector(RealField(10), y)


def inv_circular_buffer(N, ee, matching_scheme, MS, p_cross, channel, N0):
    channel = channel.split('_')[0]
    if not matching_scheme:
        return ee
    if channel == 'BEC':
        return bec_inv_circbuf(N, ee, matching_scheme, MS)
    y, counter = [], 0
    llr = -(4/N0) if channel == 'AWGN' else log(p_cross / (1 - p_cross))
    ms_operand = oo if matching_scheme=='shortening' else 0

    for a in range(N):
        if a in MS:
            yi = ee[N+counter]*llr + ee[counter]*llr if matching_scheme == 'repetition' else\
                ms_operand
            y.append(yi)
        else:
            y.append(ee[counter]*llr)
            counter += 1

    return vector(RealField(10), y)
