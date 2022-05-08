# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import sys
import csv
import datetime
import time
from random import sample
from sys import getsizeof
import gc
from scipy.stats import norm

import PC_Decoding
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Rate_Matching
import PC_Subchannel_Allocation
import HelperFunctions as HF
import PC_CRC
import test_CRC

# sage MC.py 12 0 BSC SCL
R = [1/2] # [1/2, 2/5, 1/3, 1/4,  1/5]   # Rate of the code
A_min = 12
runs = 1000
SNR = [0.6] #[0.01, 0.3, 0.2]   # really p_cross
#SNR = [1, 2, 3, 4, 5] #, 6]      # this is SNR
# SNR = [0.5, 1, 2, 3, 4]
A = int(sys.argv[1])
I_IL = int(sys.argv[2])
channel, decoder = sys.argv[3].upper(), sys.argv[4].upper()
pol = PC_CRC.get_pol(A, I_IL)
K = A + pol.degree()
n_min, n_max = 5, 10 - I_IL  # n_max = 10 for uplink, 9 for downlink.
N0 = None

for rate in R:
    E = ceil(K / rate)
    n = min(ceil(log(E, 2)), n_max)     # TODO: hvordan velges egt 'n'?
    N = 2 ** n
    QN0 = PC_Subchannel_Allocation.get_Q_N0(N)
    npc, n_wm_pc = PC_Subchannel_Allocation.get_n_pc_bits(K, E, I_IL)
    QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, E, npc, QN0)
    GN = PC_Encoder.gen_G(n)
    QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)
    del QN0; gc.collect()
    for i, snr in enumerate(SNR):
        if channel == 'AWGN':
            sigma = sqrt(1 / (2 * rate * 10 ** (snr/ 10)))
            p = 1 - norm.cdf(1 / sigma)  # error probability, from proposition 2.9
            N0 = 2*sigma**2
        BLER, BER = 0, 0
        false_negative = 0  # counts of how many times the codeword is in scout, but is not the most likely output
        for iteration in range(runs):
            if iteration % 100 == 0:
                print(iteration)
            a = random_vector(GF(2), A)
            c, G = test_CRC.CRC(a, A, pol)
            c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL=I_IL, c_seq=c)
            u = PC_Subchannel_Allocation.calc_u(N, QNI, c_ap, QNPC)
            d = vector(GF(2), u) * GN
            e = PC_Rate_Matching.circular_buffer(d, MS, matching_scheme)
            r = list(HF.channel_noise(s=e, channel=channel, p=sigma if channel == 'AWGN' else snr))
            scout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=matching_scheme, MS=MS,
                                             p_cross=snr, channel=channel + '_' + decoder, I_IL=I_IL, PI=PI, C=G)
            if decoder == 'SCL':
                if I_IL:
                    scout = scout[0].inf_bits
                    if scout == '':
                        scout = a + vector(GF(2), [1]*A)
                else:
                    crc_pass = False
                    for i, elem in enumerate(scout):
                        if test_CRC.CRC_check(elem.inf_bits, K, pol) == 0:
                            scout = vector(GF(2), elem.inf_bits)
                            crc_pass = True
                            break
                    if not crc_pass:
                        scout = vector(GF(2), scout[0].inf_bits)
            BER = BER + (a + scout[:A]).hamming_weight()
            BLER = BLER + sign((a + scout[:A]).hamming_weight())


        file_getter = channel + '_' + decoder
        with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\{file_getter}.csv', mode='a',
                  newline='') as file:
            result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            result_writer.writerow(
                [A, rate, K, N, runs, BER, BLER, false_negative, snr, I_IL, datetime.datetime.now()])
            gc.collect()
