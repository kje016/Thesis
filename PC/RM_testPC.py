# cd Desktop/Thesis/PySageMath/PC
import time

from sage.all import *

import csv
import datetime

import gc
from scipy.stats import norm

import PC_Decoding
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Rate_Matching
import PC_Subchannel_Allocation
import HelperFunctions as HF
#import test_CRC
import CRC

run_vals = [[1/2, 1, 22, 'AWGN'], [1/2, 2, 22, 'AWGN'], [1/2, 3, 22, 'AWGN'], [1/2, 4, 22, 'AWGN'],
            [1/2, 5, 22, 'AWGN']]
I_IL = 0
runs = 10000
decoder = 'SC'
rate = 1/2
snr = 2
channel = 'AWGN'

A = 22

pol = CRC.get_pol(A, I_IL)
K_min = A + pol.degree()
N0 = None
Emin = ceil(K_min / rate)
n = min(ceil(log(Emin, 2)), 10)  # TODO: hvordan velges egt 'n'?
N = 2 ** n
for elem in run_vals:
    snr = elem[1]
    K_max = floor(N*rate)
    for rep in range(N - Emin, 0, -2):
        A = 22 + (N-Emin - rep)
        K = A + pol.degree()
        N0 = None
        #E = ceil(K / rate)
        U = rep
        E = N-rep

        QN0 = PC_Subchannel_Allocation.get_Q_N0(N)
        npc, n_wm_pc = PC_Subchannel_Allocation.get_n_pc_bits(K, Emin, I_IL)
        QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, Emin, npc, QN0, U)
        GN = PC_Encoder.gen_G(n)
        QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)
        del QN0; gc.collect()
        if channel == 'AWGN':
            sigma = sqrt(1 / (2 * rate * 10 ** (snr / 10)))
            p = 1 - norm.cdf(1 / sigma)  # error probability, from proposition 2.9
            N0 = 2 * sigma ** 2

        BLER, BER = 0, 0
        start_time = time.time()
        for iteration in range(runs):
            if iteration % 100 == 0 or iteration == runs-1:
                print(time.time()-start_time)
                start_time = time.time()
                print(f"{BER}, {BLER}")
                print(iteration)
            a = random_vector(GF(2), A)
            c, G = CRC.CRC(a, A, pol)
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
                        scout = a + vector(GF(2), [1] * A)
                else:
                    crc_pass = False
                    for i, elem in enumerate(scout):
                        if CRC.CRC_check(elem.inf_bits, K, pol) == 0:
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
        #with open(f'Users\\kristian\\PycharmProjects\\Thesis\\PC\\Tests\\{file_getter}.csv', mode='a', newline='') as file:
            result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            result_writer.writerow(
                [A, rate, K, N, runs, BER, BLER, f'U={U}', snr, I_IL, datetime.datetime.now()])
            gc.collect()