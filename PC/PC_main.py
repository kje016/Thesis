# cd Desktop/Thesis/PySageMath/PC
import numpy as np
from sage.all import *

import sys
import BEC_SCL
import LLR_SCL
import HelperFunctions as HF
import PC_Code_Block_Segmentation
import PC_CRC
import PC_Input_Bits_Interleaver
import PC_Subchannel_Allocation
import PC_Encoder
import PC_Rate_Matching
import test_CRC

import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt
from scipy.stats import norm


# sage PC_main.py 10001000000110100110101111100100 PDCCH BSC
# sage PC_main.py '#information bits' 'Bool: Bit interleaver' 'Rate' 'Channel'
# sage PC_main.py 12 0 1/2 awgn

# sage PC_main.py 12 1 1/2 awgn

SNR = vector(RealField(10), [1, 1.5, 2, 2.5, 3, 3.5, 5, 4.5, 5, 5.5, 6])

if __name__ == "__main__":
    cntr, runs = 0, 250
    A = int(sys.argv[1])
    I_IL, channel = int(sys.argv[2]), sys.argv[4].upper()
    R = [int(x) for x in sys.argv[3].split('/')]
    R = R[0] / R[1]

    sigi = 4
    sigma = vector(RealField(10), map(lambda z: sqrt(1 / (2 * R * 10 ** (z / 10))), SNR))

    N0 = 2 * sigma[sigi] ** 2
    p_cross = sigma[sigi]
    if channel == 'BSC':
        p_cross = 0.1
    n_min, n_max = 5, 10 - I_IL  # n_max = 10 for uplink, 9 for downlink.
    " Mother polar code length and rate matching selection    "
    for p in range(runs):
        a = list(random_vector(GF(2), A))
        """ Code-Block Segmentation """
        # a_ap, C = PC_Code_Block_Segmentation.main_block_segmentation(I_IL, a, A, G) # TODO: repeats bits?
        for a_elem in [a]:
            pol = PC_CRC.get_pol(A, I_IL)
            #c = PC_CRC.CRC_calc(a_elem, pol)
            #cc = PC_CRC.CRC_checksum(list(c), pol)
            c, C = test_CRC.CRC(a, A, pol)
            #print(f"c := {c}")
            tc = test_CRC.CRC_check(c, len(c), pol)
            K = len(c)
            E = ceil(K / R)
            n = min(ceil(log(E, 2)), n_max)
            N = 2 ** n
            ################################################################################
            """                           Encoding                                      """
            ################################################################################
            """ Interleaving the CRC bits """
            c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL, c, A, C)   # TODO: ctess
            """ Subchannel allocation """
            u, n_pc, n_wm_pc, QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.main(N=N, c_ap=c_ap, K=K, E=E, I_IL=I_IL, R=R)

            """ Polar Code Encoding """
            d = PC_Encoder.main_encoder(u=u, n=n, n_pc=n_pc, n_wm_pc=n_wm_pc)
            """ Rate Matching by Circular Buffer    """
            e = PC_Rate_Matching.circular_buffer(y=d, matching_set=MS, matching_scheme=matching_scheme)
            ################################################################################
            """                           Decoding                                       """
            ################################################################################
            r = list(HF.channel_noise(s=e, channel=channel, p=p_cross))

            """ Rate de-Matching Circular Buffer    """
            yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=r, matching_scheme=matching_scheme, MS=MS, p_cross=p_cross, channel=channel, N0=N0)

            ee = (vector(RealField(10), e) * 2).apply_map(lambda a: a - 1)
            ty = PC_Rate_Matching.inv_circular_buffer(N=N, ee=ee, matching_scheme=matching_scheme, MS=MS, p_cross=p_cross, channel=channel, N0=N0)
            dd = vector(RealField(10), [1 if a == 0 else -1 for a in d])
            """ SC Decoder  """
            # print(f"MS := {MS}")
            if channel == 'BEC':
                uu = BEC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)
            else:
                uu = LLR_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross, I_IL=I_IL, PI=PI, pol=pol, C=C) # TODO: testing for rate-matched non-noise
            if I_IL:
                if vector(GF(2), a) == uu[0].inf_bits:
                    cntr += 1
            else:
                for dec in uu:
                    if dec*C == 0:
                        cntr += 1
                        break
    print(f"len QNF := {len(QNF)}")
    print(f" N-E = {N-E}, {matching_scheme}")
    print(f"successful decodings := {(cntr/runs)*100}%")

