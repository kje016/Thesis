from sage.all import *

import BEC_SC
import LLR_SC
import PC_Rate_Matching
import BEC_SCL
import LLR_SCL


def PC_Decoding(r, QNF, ms, N, MS, p_cross, channel, N0, I_IL):
    """Channel de-Interleaver"""
    if not MS:
        yy = r
    else:
        yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=r, matching_scheme=ms, MS=MS, p_cross=p_cross,
                                                  channel=channel, N0=N0)

    if channel == 'BEC_SC':
        uu = BEC_SC.decoder(d=yy, N=N, frozen_set=QNF)
    elif channel == 'BEC_SCL':
        uu = BEC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)

    if channel[-2:] == 'SC':
        uu = LLR_SC.decoder(d=yy, N=N, frozen_set=QNF)
    else:
        if channel.split('_')[0] == 'BSC':
            uu = LLR_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross, I_IL=I_IL, PI=None, C=None, pol=None)
    del yy
        #uu = LLR_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)

    return uu
