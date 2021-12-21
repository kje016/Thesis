# cd Desktop/Thesis/PySageMath/PC
from sage.all import *

import sys

import BSC_SCL
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
* a = sys.argv[1] length of the original information bits before they are treated in the polar code process
* physical_channel = sys.argv[2] is the phsyical channel used for encryption
"""

# sage PC_main.py 10001000000110100110101111100100 PDCCH BSC
# sage PC_main.py 4 PDCCH BSC
if __name__ == "__main__":
    G = 200  # TODO: have G to not be hard-coded
    physical_channel = sys.argv[2]  # PUCCH
    soft_channel = sys.argv[3]
    p_cross = 0.2   # TODO: for now, only BSC
    # a, A = list(random_vector(GF(2), int(sys.argv[1]))), int(sys.argv[1])
    a = [1, 0, 1, 0]
    A = len(a)
    print(f"a:= {a}")

    " Mother polar code length and rate matching selection    "
    E = ceil(G/2)           # code length (after rate matching)
    K = A + 24         # '+ 24' since using the crc24 polnomial TODO: how to decide K?

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
        ################################################################################
        """                           Encoding                                      """
        ################################################################################
        """ adding CRC bits"""
        c, polynomial = PC_CRC.main_CRC(a_elem, physical_channel)
        #print(f"CRC_codeword := {c} \n polynomial used := {polynomial.degree()}")
        """ Interleaving the CRC bits """
        I_IL, c_ap = PC_Input_Bits_Interleaver.main_bit_interleaver(physical_channel, c)
        #print(f" interleaver CRC := {c_ap}")
        """ Subchannel allocation """
        #print(f"N={N}, K={K}, E={E}, U={N-E} ")
        u, n_pc, n_wm_pc, QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.main(N=N, c_ap=c_ap, A=A, E=E, channel=physical_channel)
        #print(f"u := {u}")

        """ Polar Code Encoding """
        d = PC_Encoder.main_encoder(u=u, N=N, n_pc=n_pc, n_wm_pc=n_wm_pc)
        #print(f"codeword: := {d}")
        """ Sub-Block Interleaver"""
        #y = PC_Sub_Block_Interleaver.main_sub_block_interleaver(d=d, N=N)
        #print(f"d interleaved := {y}")

        """ Rate Matching by Circular Buffer    """
        e = PC_Rate_Matching.circular_buffer(y=d, matching_scheme=matching_scheme, matching_set=MS)
        #e, matching_scheme = PC_Rate_Matching.main_circular_buffer(y, E, N, N - E, K, QNF, QNI)
        print(f"e:=\n {e}")


        """ Channel Interleaver """
        f, I_BIL = PC_Channel_Interleaver.main_channel_interleaver(E=E, e=e, channel=physical_channel)
        print(f"f :=\n {f}")
        g.append(e)

        ################################################################################
        """                           Decoding                                       """
        ################################################################################
        """Channel de-Interleaver"""
        #ee = PC_Channel_Interleaver.inv_channel_interleaver(f=f, E=E, I_BIL=I_BIL)
        ee = d

        """ Rate de-Matching Circular Buffer    """
        yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=ee, matching_scheme=matching_scheme, MS=MS, p_cross=p_cross)
        #yy = PC_Rate_Matching.inv_circular_buffer(e=f, matching_scheme=matching_scheme, U=N-E, channel=soft_channel, p_cross=p_cross)
        print(f"yy :=\n {yy}")
        print(f"MS :=\n {MS}")
        """ Sub-Block de-Interleaver    """
        #dd = PC_Sub_Block_Interleaver.inv_sub_block_interleaver(yy, len(yy))
        dd = yy
        """ SC Decoder  """
        uu = BSC_SCL.BPSK_decoder(d=dd, N=N, frozen_set=QNF, p_cross=p_cross)
        #uu = LL_PC_BSC.SC_decoder(dd, N, list(QNF))
        print(f"uu := {uu}")
        breakpoint()

    g = g if C==1 else [elem for sublist in g for elem in sublist]
    if G%2 != 0:
        g.append(0)
