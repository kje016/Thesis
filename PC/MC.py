# cd Desktop/Thesis/PySageMath/PC
# cd Desktop/Thesis/PC
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

"""
[1/2, 1, 40, 'AWGN', 'SCL'], [1/2, 2, 40, 'AWGN', 'SCL'], [1/2, 3, 40, 'AWGN', 'SCL'], [1/2, 4, 40, 'AWGN', 'SCL'],
[1/2, 5, 40, 'AWGN', 'SCL'],

run_vals = [
[1/2, 0.5, 1, 'BEC', 'SCL'], [1/2, 0.4, 1, 'BEC', 'SCL'],[1/2, 0.3, 1, 'BEC', 'SCL'],[1/2, 0.2, 1, 'BEC', 'SCL'],
[1/2, 0.1, 1, 'BEC', 'SCL'],

[2/5, 0.6, 8, 'BEC', 'SCL'],[2/5, 0.5, 8, 'BEC', 'SCL'],[2/5, 0.4, 8, 'BEC', 'SCL'],[2/5, 0.3, 8, 'BEC', 'SCL'],
[2/5, 0.2, 8, 'BEC', 'SCL'],[2/5, 0.1, 8, 'BEC', 'SCL'],

[2/5, 1, 19, 'AWGN', 'SCL'],[2/5, 2, 19, 'AWGN', 'SCL'],[2/5, 3, 19, 'AWGN', 'SCL'],[2/5, 4, 19, 'AWGN', 'SCL'],
[2/5, 5, 19, 'AWGN', 'SCL'],
[2/5, 0.14, 19, 'BSC', 'SCL'],[2/5, 0.12, 19, 'BSC', 'SCL'],[2/5, 0.1, 19, 'BSC', 'SCL'],[2/5, 0.08, 19, 'BSC', 'SCL'],
[2/5, 0.06, 19, 'BSC', 'SCL'],[2/5, 0.04, 19, 'BSC', 'SCL'],
[2/5, 0.6, 19, 'BEC', 'SCL'],[2/5, 0.5, 19, 'BEC', 'SCL'],[2/5, 0.4, 19, 'BEC', 'SCL'],[2/5, 0.3, 19, 'BEC', 'SCL'],
[2/5, 0.2, 19, 'BEC', 'SCL'],
[1/2, 1, 21, 'AWGN', 'SCL'],[1/2, 2, 21, 'AWGN', 'SCL'],[1/2, 3, 21, 'AWGN', 'SCL'],[1/2, 4, 21, 'AWGN', 'SCL'],
[1/2, 5, 21, 'AWGN', 'SCL'],

run_vals = [
[2/5, 0.12, 27, 'BSC', 'SCL'], [2/5, 0.1, 27, 'BSC', 'SCL'],[2/5, 0.08, 27, 'BSC', 'SCL'],[2/5, 0.06, 27, 'BSC', 'SCL'],
[2/5, 0.04, 27, 'BSC', 'SCL'],

[2/5, 1, 27, 'AWGN', 'SCL'],[2/5, 2, 27, 'AWGN', 'SCL'],[2/5, 3, 27, 'AWGN', 'SCL'],[2/5, 4, 27, 'AWGN', 'SCL'],
[2/5, 5, 27, 'AWGN', 'SCL'],

[1/2, 0.1, 40, 'BSC', 'SCL'],[1/2, 0.08, 40, 'BSC', 'SCL'],[1/2, 0.06, 40, 'BSC', 'SCL'],[1/2, 0.04, 40, 'BSC', 'SCL'],
[1/2, 0.02, 40, 'BSC', 'SCL'],

[1/2, 1, 40, 'AWGN', 'SCL'],[1/2, 2, 40, 'AWGN', 'SCL'],[1/2, 3, 40, 'AWGN', 'SCL'],[1/2, 4, 40, 'AWGN', 'SCL'],
[1/2, 5, 40, 'AWGN', 'SCL'],
[1/2, 0.5, 8, 'BEC', 'SCL', 1],[1/2, 0.4, 8, 'BEC', 'SCL', 1],[1/2, 0.3, 8, 'BEC', 'SCL', 1],[1/2, 0.2, 8, 'BEC', 'SCL', 1],
[1/2, 0.1, 8, 'BEC', 'SCL', 1],
[1/2, 0.5, 21, 'BEC', 'SCL',0],[1/2, 0.4, 21, 'BEC', 'SCL',0],[1/2, 0.3, 21, 'BEC', 'SCL',0],[1/2, 0.2, 21, 'BEC', 'SCL',0],
[1/2, 0.1, 21, 'BEC', 'SCL',0],
]
"""
run_vals = [
['', 1, 21, 'AWGN', 'SC'],['', 2, 21, 'AWGN', 'SC'],['', 3, 21, 'AWGN', 'SC'],['', 4, 21, 'AWGN', 'SC'],
['', 5, 21, 'AWGN', 'SC']
]

runs = 10000
I_IL = 0
#decoder = 'SC' or 'SCL'
for elem in run_vals:
    rate = 1/2
    snr = elem[1]
    A = elem[2]
    channel = elem[3] # 'AWGN'/ 'BSC' / 'BEC'
    decoder = elem[4] # 'SC'/ 'SCL'

    pol = CRC.get_pol(A, I_IL)
    K = A + pol.degree()
    PI = PC_Input_Bits_Interleaver.get_pi(K)
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

    BLER, BER, ET = 0, 0, 0
    start_time = time.time()
    print(elem, len(MS))
    for iteration in range(runs):
        if iteration % 100 == 0 or iteration == runs-1:
            print(time.time()-start_time)
            start_time = time.time()
            print(f"{BER}, {BLER}, {ET}")
            print(iteration)
            a = random_vector(GF(2), A)
            #a = vector(GF(2), [1])
            c, H = CRC.CRC(a, A, pol, I_IL, PI)
            c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL=I_IL, c_seq=c)
            #print(f'c:\n{c_ap}')
            u = PC_Subchannel_Allocation.calc_u(N, QNI, c_ap, QNPC)
            d = vector(GF(2), u) * GN
            #print(f'd:\n{d}')
            e = PC_Rate_Matching.circular_buffer(d, MS, matching_scheme)
        r = list(HF.channel_noise(s=e, channel=channel, p=sigma if channel == 'AWGN' else snr))
        scout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=matching_scheme, MS=MS,
                                        p_cross=snr, channel=channel + '_' + decoder, I_IL=I_IL, PI=PI, H=H)
        if decoder == 'SCL':
            if I_IL:
                scout = scout[0].inf_bits
                if scout == '':
                    scout = a + vector(GF(2), [1] * A)
                    ET += 1

        BER = BER + (a + scout[:A]).hamming_weight()
        BLER = BLER + sign((a + scout[:A]).hamming_weight())
    file_getter = channel + '_' + decoder
    #with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\{file_getter}.csv', mode='a',
    #          newline='') as file:
    with open(f'Tests/{file_getter}.csv', mode='a', newline='') as file:
        result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        result_writer.writerow(
            [A, rate, K, N, runs, BER, BLER, f'U={len(MS)},E:{E},H-check', snr, I_IL, datetime.datetime.now()])
        gc.collect()
