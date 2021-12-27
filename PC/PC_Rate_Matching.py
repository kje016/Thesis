# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions as HF


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
    if matching_scheme == "repetition":
        return None
    else:
        bit_len = len("{0:b}".format(len(QN0)-1))
        bnm = [int(("0"*(bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1], 2) for a in QN0]
        if matching_scheme == "puncturing":
            matching_set = set(bnm[:U])
        else: # matching_scheme == "shortening"
            matching_set = set(bnm[-U::])
        return matching_set


def circular_buffer(y, matching_scheme, matching_set):
    e, counter = [], 0
    for a in range(len(y)):
        if a not in matching_set:
            e.append(y[counter])
        counter += 1
    return e


def main_circular_buffer(y, E, N, U, K, QNF, QNI):
    matchin_scheme = "puncturing"
    matching_scheme = matching_selection(E, N, K)
    N, K, U = 8, 4, 2
    y = [0, 1, 0, 1, 1, 0, 1, 0]

    e = []
    if matching_scheme == "repetition":
        e = y + y[0:U]  # first U bits are transmitted twice
    else:
        rnm = [a-1 for a in HF.get_inf_pos(8)]
        bit_len = len("{0:b}".format(rnm[-1]))
        bnm = [int(("0"*(bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1], 2) for a in rnm]
        print(f"Bnm := {bnm}")
        if matchin_scheme == "puncturing":
            matching_set = set(bnm[:U])
        else: # matching_scheme == "shortening"
            matching_set = set(bnm[-U::])

        QNI, QNF = QNI-set(matching_set), QNF + set(matching_set)
        y_get = 0
        for a in range(N):
            if a in QNI:
                e.append(y[y_get])
                y_get += 1
            if len(e) == N-U:
                break
        return e, (matching_set, matching_scheme)


def inv_circular_buffer(N, ee, matching_scheme, MS, p_cross):
    y, counter = [], 0
    F, llr1 = RealField(10), log(p_cross / (1 - p_cross))
    if matching_scheme == "shortening":
        for a in range(N):
            if a in MS:
                y.append(-oo)   # -oo since the llr1 will be negative
            else:
                y.append(ee[counter])
                counter += 1
    elif matching_scheme == "puncturing":
            for a in range(N):
                if a in MS:
                    y.append(0.5)  # 0.5 so the llr equals 0 in the lambda function below
                else:
                    y.append(ee[counter])
                    counter += 1
    return vector(F, list(map(lambda x: x - 1, 2 * vector(F, y)))) * llr1


# p_cross gjelder egt bare for BSC. evnt bytte til channel_metric/channel_characteristic
def LLR_fun(e, channel, p_cross):
    F = RealField(10)
    if channel == 'BSC':
        llr1 = log(p_cross/(1 - p_cross))
        return vector(F, list(map(lambda x: x - 1, 2 * vector(F, e)))) * llr1
