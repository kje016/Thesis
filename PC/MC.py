# cd Desktop/Thesis/PySageMath/PC
import time
import os

from sage.all import *
import sys
import csv
import datetime
#import time
#from random import sample
#from sys import getsizeof
import gc
from scipy.stats import norm

import PC_Decoding
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Rate_Matching
import PC_Subchannel_Allocation
import HelperFunctions as HF
import CRC



run_vals = [
[1/2, 1, 40, 'AWGN', 'SCL'], [1/2, 2, 40, 'AWGN', 'SCL'], [1/2, 3, 40, 'AWGN', 'SCL'], [1/2, 4, 40, 'AWGN', 'SCL'],
[1/2, 5, 40, 'AWGN', 'SCL'],

[1/2, 0.1, 40, 'BSC', 'SCL'], [1/2, 0.08, 40, 'BSC', 'SCL'], [1/2, 0.06, 40, 'BSC', 'SCL'], [1/2, 0.04, 40, 'BSC', 'SCL'],
[1/2, 0.02, 40, 'BSC', 'SCL'],

[1/3, 1, 18, 'AWGN', 'SCL'], [1/3, 2, 18, 'AWGN', 'SCL'], [1/3, 3, 18, 'AWGN', 'SCL'], [1/3, 4, 18, 'AWGN', 'SCL'],
[1/3, 5, 18, 'AWGN', 'SCL'],

[1/3, 0.14, 18, 'BSC', 'SCL'], [1/3, 0.12, 18, 'BSC', 'SCL'], [1/3, 0.1, 18, 'BSC', 'SCL'],
[1/3, 0.08, 18, 'BSC', 'SCL'], [1/3, 0.06, 18, 'BSC', 'SCL']
]
"""
run_vals = [ [2/5, 0.14, 19, 'BSC'], [2/5, 0.12, 19, 'BSC'], [2/5, 0.10, 19, 'BSC'],
[2/5, 0.08, 19, 'BSC'], ]
"""
I_IL = 1
runs = 20000
#decoder = 'SC'
for elem in run_vals:
    rate = elem[0]
    snr = elem[1]
    A = elem[2]
    channel = elem[3]
    decoder = elem[4]

    pol = CRC.get_pol(A, I_IL)
    K = A + pol.degree()
    N0 = None
    false_negative = 'NaN'

    E = ceil(K / rate)
    n = min(ceil(log(E, 2)), 10)     # TODO: hvordan velges egt 'n'?
    N = 2 ** n
    QN0 = PC_Subchannel_Allocation.get_Q_N0(N)
    npc, n_wm_pc = PC_Subchannel_Allocation.get_n_pc_bits(K, E, I_IL)
    QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, E, npc, QN0, 0)
    GN = PC_Encoder.gen_G(n)
    QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)
    del QN0; gc.collect()

    if channel == 'AWGN':
        sigma = sqrt(1 / (2 * rate * 10 ** (snr / 10)))
        p = 1 - norm.cdf(1 / sigma)  # error probability, from proposition 2.9
        N0 = 2 * sigma ** 2

    BLER, BER = 0, 0
    start_time = time.time()
    print(elem, N-E)
    for iteration in range(runs):
        if iteration % 100 == 0 or iteration == runs-1:
            print(time.time()-start_time)
            start_time = time.time()
            print(f"{BER}, {BLER}")
            print(iteration)
            a = random_vector(GF(2), A)
            #a = vector(GF(2), [1])
            c, G = CRC.CRC(a, A, pol)
            c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL=I_IL, c_seq=c)
            #print(f'c:\n{c_ap}')
            u = PC_Subchannel_Allocation.calc_u(N, QNI, c_ap, QNPC)
            d = vector(GF(2), u) * GN
            #print(f'd:\n{d}')
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
    #with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\{file_getter}.csv', mode='a',
    #          newline='') as file:
    with open(f'Tests/{file_getter}.csv', mode='a', newline='') as file:
        result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        result_writer.writerow(
            [A, rate, K, N, runs, BER, BLER, f'RM testing: U={len(MS)}', snr, I_IL, datetime.datetime.now()])
        gc.collect()
