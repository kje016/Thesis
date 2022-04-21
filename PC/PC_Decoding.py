from sage.all import *

import BEC_SC
import LLR_SC
import PC_Rate_Matching
import BEC_SCL
import LLR_SCL


def PC_Decoding(r, QNF, ms, N, MS, p_cross, channel, N0, I_IL):
    """Channel de-Interleaver"""
    yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=r, matching_scheme=ms, MS=MS, p_cross=p_cross,
                                              channel=channel, N0=N0)
    if channel[-2:] == 'SC':
        uu = BEC_SC.decoder(d=yy, N=N, frozen_set=QNF) if channel[:3]=='BEC' else\
            LLR_SC.decoder(d=yy, N=N, frozen_set=QNF)
    else:
        uu = BEC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross) if channel[:3]== 'BEC' else\
            LLR_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross, I_IL=I_IL, PI=None, C=None)
    return uu
