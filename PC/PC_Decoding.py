from sage.all import *

import BEC_SC
import LLR_SC
import PC_Rate_Matching
import BEC_SCL
import LLR_SCL
import PC_CRC


def PC_Decoding(r, QNF, ms, N, MS, p_cross, channel, N0, I_IL):
    """Channel de-Interleaver"""
    if channel[-2:] == 'SC':
        yy = r
        del r
    else:
        yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=r, matching_scheme=ms, MS=MS, p_cross=p_cross, channel=channel, N0=N0)
        del r

    if channel == 'BEC_SC':
        uu = BEC_SC.decoder(d=yy, N=N, frozen_set=QNF)
    elif channel == 'BEC_SCL':
        uu = BEC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)

    if channel[-2:] == 'SC':
        uu = LLR_SC.decoder(d=yy, N=N, frozen_set=QNF)
    else:
        if channel.split('_')[0] == 'BSC':
            uu = LLR_SCL.decoder(d=yy,N=N, frozen_set=QNF,p_cross=p_cross,I_IL=I_IL,PI=None, C=None, pol=None)
    del yy
        #uu = LLR_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)

    """
    for dec in uu:
        if vector(GF(2), PC_CRC.bit_long_division(list(dec.inf_bits), pol)) == 0:
            return dec.inf_bits[:-pol.degree()], True
    """
    #return uu[0].inf_bits[:-pol.degree()], False
    return uu

