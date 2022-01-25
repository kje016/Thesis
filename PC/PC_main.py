# cd Desktop/Thesis/PySageMath/PC
from sage.all import *

import sys

import BSC_SCL
import HelperFunctions
import PC_Code_Block_Segmentation
import PC_CRC
import PC_CRC
import PC_Input_Bits_Interleaver
import PC_Subchannel_Allocation
import PC_Encoder
import PC_Sub_Block_Interleaver
import PC_Rate_Matching
import PC_Channel_Interleaver
import test_CRC

"""calling PC_main and its arguments 
* a = sys.argv[1] length of the original information bits before they are treated in the polar code process
* physical_channel = sys.argv[2] is the phsyical channel used for encryption
"""

# sage PC_main.py 10001000000110100110101111100100 PDCCH BSC
# sage PC_main.py 12 PDCCH BSC
if __name__ == "__main__":
    cntr = 0
    runs = 20
    for p in range(runs):
        G = 200  # TODO: have G to not be hard-coded
        physical_channel = sys.argv[2]  # PUCCH
        soft_channel = sys.argv[3]
        p_cross = 0.1  # TODO: for now, only BSC
        # a = [1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1]
        #a = [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0]
        # a = a + [0]*(12-len(a))
        #A = len(a)
        a, A = list(random_vector(GF(2), int(sys.argv[1]))), int(sys.argv[1])
        print(f"a:= {a, A}")

        " Mother polar code length and rate matching selection    "
        E = ceil(G/2)           # code length (after rate matching)
        K = A + 24         # '+ 24' since using the crc24 polnomial TODO: how to decide K?

        " computing n"
        n_min, n_max = 5, 9    # n_max = 10 for uplink, 9 for downlink.
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
            ctess, poltess = test_CRC.CRC(a, A)
            print(f"CRC are equal := {c==ctess}")
            # TODO: polynomial set to 24 for interleaver testing
            """ Interleaving the CRC bits """
            I_IL, c_ap = PC_Input_Bits_Interleaver.main_bit_interleaver(physical_channel, c, A)

            # TODO: I_IL set to True for interleaver testing
            """ Subchannel allocation """
            u, n_pc, n_wm_pc, QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.main(N=N, c_ap=c_ap, A=A, E=E, channel=physical_channel)

            """ Polar Code Encoding """
            d = PC_Encoder.main_encoder(u=u, N=N, n_pc=n_pc, n_wm_pc=n_wm_pc)

            """ Sub-Block Interleaver"""
            #y = PC_Sub_Block_Interleaver.main_sub_block_interleaver(d=d, N=N)
            #print(f"d interleaved := {y}")

            """ Rate Matching by Circular Buffer    """
            e = PC_Rate_Matching.circular_buffer(y=d, matching_set=MS)
            #print(f"e:=\n {e}")


            """ Channel Interleaver """
            f, I_BIL = PC_Channel_Interleaver.main_channel_interleaver(E=E, e=e, channel=physical_channel)
            #print(f"f :=\n {f}")
            g.append(e)

            ################################################################################
            """                           Decoding                                       """
            ################################################################################
            r = list(HelperFunctions.channel_noise(s=e, channel=soft_channel, p=p_cross))
            """Channel de-Interleaver"""
            # N, matching_scheme, d, p_cross, MS, QNF = 8, 'puncturing', [0, 0, 0, 1, 0, 0, 1, 0], 0.2, [0, 4], [0, 1, 2, 4]
            # ee = PC_Channel_Interleaver.inv_channel_interleaver(f=f, E=E, I_BIL=I_BIL)
            ee = r

            """ Rate de-Matching Circular Buffer    """
            yy = PC_Rate_Matching.inv_circular_buffer(N=N, ee=ee, matching_scheme=matching_scheme, MS=MS, p_cross=p_cross)
            #yy = PC_Rate_Matching.inv_circular_buffer(e=f, matching_scheme=matching_scheme, U=N-E, channel=soft_channel, p_cross=p_cross)

            """ Sub-Block de-Interleaver    """
            #dd = PC_Sub_Block_Interleaver.inv_sub_block_interleaver(yy, len(yy))
            dd = yy
            """ SC Decoder  """
            #uu = BSC_SCL.BPSK_decoder(d=dd, N=N, frozen_set=QNF, p_cross=p_cross)
            #uu = LL_PC_BSC.SC_decoder(dd, N, list(QNF))
            uu = BSC_SCL.decoder(d=dd, N=N, frozen_set=QNF, p_cross=p_cross)
            for dec in uu:
                #breakpoint()
                if vector(PC_CRC.bit_long_division(list(dec.inf_bits), polynomial)) == 0:
                    print(cntr)
                    cntr += 1
                    break
            #print(f"uu := {uu}")
    print(f"correct := {cntr}")

    g = g if C==1 else [elem for sublist in g for elem in sublist]
    if G%2 != 0:
        g.append(0)
