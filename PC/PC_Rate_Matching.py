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


# bit = "0" * (bit_len - len("{0:b}".format(a))) + "{0:b}".format(a)
# br = ("0" * (bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1]
def get_rm_set(U, matching_scheme, QN0):
    bit_len = len("{0:b}".format(len(QN0)-1))
    # bnm = [int(("0" * (bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1], 2) for a in list(range(len(QN0)))]
    #for a in QN0:
        #bnm.append(int(("0" * (bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1], 2))
    bnm = [int(("0" * (bit_len - len("{0:b}".format(a))) + "{0:b}".format(a))[::-1], 2) for a in QN0]
    if matching_scheme == "puncturing":
        matching_set = set(bnm[:U])
    elif matching_scheme == "shortening":
        matching_set = set(bnm[-U:])
    else: # matching_scheme == 'repetition'
        matching_set = set(bnm[:abs(U)])
    return matching_set


def circular_buffer(y, matching_set, matching_scheme):
    if matching_scheme == 'repetition':
        return y[:] + [y[a] for a in matching_set]
    e, counter = [], 0
    for a in range(len(y)):
        if a not in matching_set:
            e.append(y[counter])
        counter += 1
    return e


def bec_inv_circbuf(N, ee, matching_scheme, MS, p_cross):
    y, counter = [], 0
    F, llr1 = RealField(10), log(p_cross / (1 - p_cross))
    if matching_scheme == "shortening":
        for a in range(N):
            if a in MS:
                y.append(oo)
            else:
                y.append(2 if ee[counter] == 2 else ee[counter]*llr1*oo)
            counter += 1
    elif matching_scheme == "puncturing":
        for a in range(N):
            if a in MS:
                y.append(2 if ee[counter] == 2 else ee[counter] * -oo)
            else:
                y.append(ee[counter]*llr1)
            counter += 1
    elif matching_scheme == 'repetition':
        getter = 0
        for a in range(N):
            if a in MS:
                y.append( (ee[N+getter] if ee[counter] == 2 else ee[counter])*llr1 )
                getter += 1
            else:
                y.append(ee[counter]*llr1)
            counter += 1
    return vector(F, y)


def inv_circular_buffer(N, ee, matching_scheme, MS, p_cross, channel):
    if channel == 'BEC':
        return bec_inv_circbuf(N, ee, matching_scheme, MS, p_cross)
    y, counter = [], 0
    F, llr1 = RealField(10), log(p_cross / (1 - p_cross))
    if matching_scheme == "shortening":
        for a in range(N):
            if a in MS:
                y.append(oo)
            else:
                if channel == 'AWGN':
                    y.append(ee[counter]*-1)
                else:
                    y.append(ee[counter]*llr1)
                counter += 1
    elif matching_scheme == "puncturing":
        for a in range(N):
            if a in MS:
                y.append(0)  # 0.5 so the llr equals 0 in the lambda function below
            else:
                if channel == 'AWGN':
                    y.append(ee[counter] * -1)
                else:
                    y.append(ee[counter]*llr1)
                counter += 1
    elif matching_scheme == 'repetition':
        getter = 0
        for a in range(N):
            if a in MS:
                y.append( (ee[N+getter]*llr1 + ee[counter]*llr1))
                print(ee[N+getter], ee[counter], a, getter)
                getter += 1
            else:
                y.append(ee[counter]*llr1)
                counter += 1

    return vector(F, y)
    # return vector(F, list(map(lambda x: x - 1, 2 * vector(F, y)))) * llr1

