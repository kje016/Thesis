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
    return e
