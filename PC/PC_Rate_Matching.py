# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


def bit_selection(E, N, K):
    if N > E:
        if K / E <= 7 / 16:
            rate_matching = "puncturing"
        else:
            rate_matching = "shortening"
    else:
        rate_matching = "repetition"
    return rate_matching


def main_circular_buffer(y, E, N, U, K):
    matching_scheme = bit_selection(E, N, K)
    if matching_scheme == "puncturing":
        e = y[U:]    # first U bits are not transmitted
    elif matching_scheme == "shortening":
        e = y[0:len(y)-U]    # last U bits are not transmitted
    else:
        e = y + y[0:U]       # first U bits are transmitted twice
    return e, matching_scheme


# p_cross gjelder egt bare for BSC. evnt bytte til channel_metric/channel_characteristic
def LLR_fun(e, channel, p_cross):
    F = RealField(10)
    if channel == 'BSC':
        llr1 = log(p_cross/(1 - p_cross))
        return vector(F, list(map(lambda x: x - 1, 2 * vector(F, e)))) * llr1


def inv_circular_buffer(e, matching_scheme, U, channel, p_cross):
    e = list(LLR_fun(e, channel, p_cross))
    if matching_scheme == 'shortening':
        # the corresponding LLRs set to infinity and appended to the set of E received LLRs
        y = e + [+Infinity]*U
    elif matching_scheme == "puncturing":
        # the corresponding LLRs can be set to infinity
        # and prepended to the set of E received LLRs
        y = [0]*U + e
    elif matching_scheme == 'repetition':   # repetition
        # the LLRs pertaining to the replicas of each bit in sequence y
        # may be accumulated, in order to obtain a corresponding sequence
        # of N LLRs.
        e1, e2 = e[:U], e[U:]
        res = list(map(lambda x, y: x+y, e1, e2))
        y = res + e1[len(res):]
    else:
        return e
    return y
