from sage.all import *

import PC_Rate_Matching
import BEC_SCL
import BSC_SCL
import PC_CRC


def PC_Decoding(r, QNF, ms, N, MS, p_cross, channel, N0, pol):
    """Channel de-Interleaver"""
    yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=r, matching_scheme=ms, MS=MS, p_cross=p_cross, channel=channel, N0=N0)

    if channel == 'BEC':
        uu = BEC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)
    else:
        uu = BSC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)

    for dec in uu:
        if vector(GF(2), PC_CRC.bit_long_division(list(dec.inf_bits), pol)) == 0:
            return dec.inf_bits[:-pol.degree()], True
    return uu[0].inf_bits[:-pol.degree()], False

