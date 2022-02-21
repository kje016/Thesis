from sage.all import *

import PC_CRC
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Rate_Matching
import PC_Subchannel_Allocation


def PC_Encoding(a, A, R, I_IL):
    n_min, n_max = 5, 10 - I_IL  # n_max = 10 for uplink, 9 for downlink.
    # a_ap, C = PC_Code_Block_Segmentation.main_block_segmentation(I_IL, a, A, G) # TODO: repeats bits?
    for a_elem in [a]:
        pol = PC_CRC.get_pol(A, I_IL)
        c = PC_CRC.CRC_calc(a_elem, pol)
        K = len(c)
        E = ceil(K/R)
        n = min(ceil(log(E, 2)), n_max)
        N = 2 ** n
        """ Interleaving the CRC bits """
        c_ap = PC_Input_Bits_Interleaver.interleaver(I_IL, c, A)  # TODO: ctess
        """ Subchannel allocation """
        u, n_pc, n_wm_pc, QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.main(N=N, c_ap=c_ap, K=K, E=E,
                                                                                        I_IL=I_IL, R=R)
        """ Polar Code Encoding """
        d = PC_Encoder.main_encoder(u=u, n=n, n_pc=n_pc, n_wm_pc=n_wm_pc)
        """ Rate Matching by Circular Buffer    """
        e = PC_Rate_Matching.circular_buffer(y=d, matching_set=MS, matching_scheme=matching_scheme)
    return N, e, QNF, matching_scheme, MS, pol
