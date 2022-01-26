# cd Desktop/Thesis/PySageMath/PC
from sage.all import *

import sys

import BSC_SCL
import HelperFunctions as HF
import PC_Code_Block_Segmentation
import PC_CRC
import PC_CRC
import PC_Input_Bits_Interleaver
import PC_Subchannel_Allocation
import PC_Encoder
import PC_Rate_Matching
import test_CRC


# sage PC_main.py 10001000000110100110101111100100 PDCCH BSC
# sage PC_main.py '#information bits' 'Bool: Bit interleaver' 'Rate' 'Channel'
# sage PC_main.py 12 0 1/2 BSC
if __name__ == "__main__":
    cntr, runs = 0, 20
    A = int(sys.argv[1])
    I_IL, channel, p_cross = int(sys.argv[2]), sys.argv[4].upper(), 0.1
    R = [int(x) for x in sys.argv[3].split('/')]
    R = R[0] / R[1]

    n_min, n_max = 5, 10 - I_IL  # n_max = 10 for uplink, 9 for downlink.

    " Mother polar code length and rate matching selection    "
    for p in range(runs):
        E = ceil(A / R)
        n = max(min(ceil(log(E, 2)), n_max), n_min)
        print(f" n := {n}")
        N = 2 ** n
        G = N
        a = list(random_vector(GF(2), A))
        """ Code-Block Segmentation """
        a_ap, C = PC_Code_Block_Segmentation.main_block_segmentation(I_IL, a, A, G) # TODO: repeats bits?
        for a_elem in a_ap:
            pol = PC_CRC.get_pol(len(a_elem), I_IL)
            c = PC_CRC.main_CRC(a_elem, pol)
            ctess = test_CRC.CRC(a, A, pol)
            ################################################################################
            """                           Encoding                                      """
            ################################################################################
            """ Interleaving the CRC bits """
            c_ap = PC_Input_Bits_Interleaver.main_bit_interleaver(I_IL, c, A)
            """ Subchannel allocation """
            u, n_pc, n_wm_pc, QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.main(N=N, c_ap=c_ap, A=A, E=E, I_IL=I_IL)
            print(matching_scheme)
            """ Polar Code Encoding """
            d = PC_Encoder.main_encoder(u=u, N=N, n_pc=n_pc, n_wm_pc=n_wm_pc)

            """ Rate Matching by Circular Buffer    """
            e = PC_Rate_Matching.circular_buffer(y=d, matching_set=MS)

            ################################################################################
            """                           Decoding                                       """
            ################################################################################
            breakpoint()
            r = list(HF.channel_noise(s=e, channel=channel, p=p_cross))
            """Channel de-Interleaver"""
            # ee = PC_Channel_Interleaver.inv_channel_interleaver(f=f, E=E, I_BIL=I_BIL)
            ee = (vector(RealField(10), e)*2).apply_map(lambda a: a-1)

            """ Rate de-Matching Circular Buffer    """
            yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=r, matching_scheme=matching_scheme, MS=MS, p_cross=p_cross)

            """ SC Decoder  """
            if channel == 'BEC':
                pass
            else:
                uu = BSC_SCL.decoder(d=yy, N=N, frozen_set=QNF, p_cross=p_cross)
            for dec in uu:
                if dec.inf_bits == c:
                    print(f"CRC_check := {PC_CRC.bit_long_division(list(dec.inf_bits), pol), matching_scheme}")
                    cntr += 1
                    break
            #print(f"uu := {uu}")
    print(f"correct := {(cntr/runs)*100}%")

