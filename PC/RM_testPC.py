# cd Desktop/Thesis/PySageMath/PC
import time

from sage.all import *

import csv
import datetime

import gc
from scipy.stats import norm
import numpy as np

import PC_Decoding
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Rate_Matching
import PC_Subchannel_Allocation
import HelperFunctions as HF
#import test_CRC
import CRC

run_vals = [
[0.5, 'BEC'], [0.4, 'BEC'], [0.3, 'BEC'], [0.2, 'BEC'],[0.1, 'BEC'],
[0.02, 'BSC'],[0.04, 'BSC'], [0.06, 'BSC'], [0.08, 'BSC'], [0.1, 'BSC'],
[1, 'AWGN'], [2, 'AWGN'], [3, 'AWGN'], [4, 'AWGN'], [5, 'AWGN'],

[0.5, 'BEC'], [0.4, 'BEC'], [0.3, 'BEC'], [0.2, 'BEC'],[0.1, 'BEC'],
[0.02, 'BSC'],[0.04, 'BSC'], [0.06, 'BSC'], [0.08, 'BSC'], [0.1, 'BSC'],
[1, 'AWGN'], [2, 'AWGN'], [3, 'AWGN'], [4, 'AWGN'], [5, 'AWGN'],
]
I_IL = 0
PI = []
runs = 10000
decoder = 'SCL'
rate = 1/2

P = 11
A_min = 22
E_min, E_max = 64+1, 128
K_min, K_max = E_min*rate, E_max*rate
A_min, A_max = K_min-P, K_max-P

pol = CRC.get_pol(A_min, I_IL)

N0 = None

n = min(ceil(log(E_min, 2)), 10)  # TODO: hvordan velges egt 'n'?
N = 2 ** n
QN0 = PC_Subchannel_Allocation.get_Q_N0(N)
npc, n_wm_pc = PC_Subchannel_Allocation.get_n_pc_bits(K_min, E_min, I_IL)
GN = PC_Encoder.gen_G(n)

K_max = floor(N*rate)
E_loop = list(np.arange(8, E_max-E_min, 16))
E_loop.append(E_max-E_min+1)
E_loop = [E_loop[0]]
for elem in run_vals:
    snr = elem[0]
    channel = elem[1]
    print(elem)
    for U in E_loop:
        print(f'U:={U}')
        E = N-U
        K = floor(E*rate)
        A = K - P
        N0 = None
        QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, E, npc, QN0, U)
        QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)
        #del QN0; gc.collect()
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
            c, H = CRC.CRC(a, A, pol, I_IL, PI)
            c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL=I_IL, c_seq=c)
            u = PC_Subchannel_Allocation.calc_u(N, QNI, c_ap, QNPC)
            d = vector(GF(2), u) * GN
            e = PC_Rate_Matching.circular_buffer(d, MS, matching_scheme)
            r = list(HF.channel_noise(s=e, channel=channel, p=sigma if channel == 'AWGN' else snr))
            scout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=matching_scheme, MS=MS,
                                            p_cross=snr, channel=channel + '_' + decoder, I_IL=I_IL, PI=PI, H=H)
            if decoder == 'SCL':
                if I_IL:
                    scout = scout[0].inf_bits
                    if scout == '':
                        scout = a + vector(GF(2), [1] * A)

            BER = BER + (a + scout[:A]).hamming_weight()
            BLER = BLER + sign((a + scout[:A]).hamming_weight())

        file_getter = channel + '_' + decoder
        with open(f'Tests/{file_getter}.csv', mode='a',
                  newline='') as file:
        #with open(f'Users\\kristian\\PycharmProjects\\Thesis\\PC\\Tests\\{file_getter}.csv', mode='a', newline='') as file:
            result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            result_writer.writerow(
                [A, rate, K, N, runs, BER, BLER, f'U={U}, E={E}', snr, I_IL, datetime.datetime.now()])
            gc.collect()