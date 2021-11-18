# cd Desktop/Thesis/PySageMath/PC
from sage.all import *

import sys
import PC_Code_Block_Segmentation
import PC_CRC
import PC_Input_Bits_Interleaver
import PC_Subchannel_Allocation
import PC_Encoder
import PC_Sub_Block_Interleaver
import PC_Rate_Matching
import PC_Channel_Interleaver
import LL_PC_DEC
import LL_PC_BSC

"""calling PC_main and its arguments 
* a = sys.argv[1] are the original information bits before they are treated in the polar code process
* physical_channel = sys.argv[2] is the phsyical channel used for encryption
"""

# sage PC_main.py 10001000000110100110101111100100 PDCCH
if __name__ == "__main__":
    G = 200  # TODO: have G to not be hard-coded
    physical_channel = sys.argv[2]  # PUCCH

    a, A = [int(x) for x in sys.argv[1]], len(sys.argv[1])
    A = len(a)      # a := 10001000000110100110101111100100

    " Mother polar code length and rate matching selection    "
    E = ceil(G/2)           # code length (after rate matching)
    K = len(a) + 24         #  '+ 11' since using the crc11 polnomial TODO: how to decide K?

    " computing n"
    n_min, n_max = 5, 10    # n_max = 10 for uplink, 9 for downlink.
    if E <= (9/8)*(2**(ceil(log(E, 2)-1))) and K/E < 9/16:      # from etsi ts 138 212 V15.6.0 (2019-07)
        n1 = ceil(log(E, 2))-1          # from etsi ts 138 212 V15.6.0 (2019-07)
    else:
        n1 = ceil(log(E, 2))
        # from etsi ts 138 212 V15.6.0 (2019-07)
    R_min = 1/8                         # from etsi ts 138 212 V15.6.0 (2019-07)
    n2 = ceil(log(K/R_min, 2))
    n = max(min(n1, n2, n_max), n_min)
    N = 2**n
    if physical_channel == "PBCH":
        N = 512

    """ Code-Block Segmentation """
    a_ap, C = PC_Code_Block_Segmentation.main_block_segmentation(physical_channel, a, A, G)

    g = []
    for a_elem in a_ap:
        """ adding CRC bits"""
        c, polynomial = PC_CRC.main_CRC(a_elem, physical_channel)
        #print(f"CRC_codeword := {c} \n polynomial used := {polynomial.degree()}")

        """ Interleaving the CRC bits """
        I_IL, c_ap = PC_Input_Bits_Interleaver.main_bit_interleaver(physical_channel, c)
        #print(f" interleaver CRC := {c_ap}")

        """ Subchannel allocation """
        #print(f"N={N}, K={K}, E={E}, U={N-E} ")
        u, n_pc, n_wm_pc, frozen_set = PC_Subchannel_Allocation.main(N=N, c_ap=c_ap, A=A, E=E, channel=physical_channel)
        #print(f"u := {u}")

        """ Polar Code Encoding """
        d = PC_Encoder.main_encoder(u=u, N=N, n_pc=n_pc, n_wm_pc=n_wm_pc)
        #print(f"codeword: := {d}")
        """ Sub-Block Interleaver"""
        y = PC_Sub_Block_Interleaver.main_sub_block_interleaver(d=d, N=N)
        #print(f"d interleaved := {y}")

        """ Rate Matching by Circular Buffer    """
        e, matching_scheme = PC_Rate_Matching.main_circular_buffer(y, E, N, N - E, K)
        print(f"rate-matched codeword := {e}")


        """ Channel Interleaver """
        f, I_BIL = PC_Channel_Interleaver.main_channel_interleaver(E=E, e=e, channel=physical_channel)
        print(f"f := {f}")
        g.append(e)

        ################################################################################
        """                           Inverse operations                             """
        ################################################################################
        """Channel de-Interleaver"""
        ee = PC_Channel_Interleaver.inv_channel_interleaver(f=f, E=E, I_BIL=I_BIL)

        """ Rate de-Matching Circular Buffer    """
        yy = PC_Rate_Matching.inv_circular_buffer(ee, matching_scheme, N-E)

        """ Sub-Block de-Interleaver    """
        dd = PC_Sub_Block_Interleaver.inv_sub_block_interleaver(yy, len(yy))
        #breakpoint()
        """ SC Decoder  """
        uu = LL_PC_BSC.SC_decoder(dd, N, frozen_set, 0.2)
        print(f"uu := {uu}")

    g = g if C==1 else [elem for sublist in g for elem in sublist]
    if G%2 != 0:
        g.append(0)
