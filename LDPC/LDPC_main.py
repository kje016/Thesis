# cd Desktop/Thesis/PySageMath/LDPC
# cd PycharmProjects/Thesis/LDPC
import time
import threading

from sage.all import *
import gc
import csv
import datetime
from scipy.stats import norm
import sys

import CRC
import Parameter_Functions as PF
import LDPC_Encoding
import LDPC_Rate_Matching
import LDPC_HelperFunctions as HF
import BG_main

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

lss = {0: [2, 4, 8, 16, 32, 64, 128, 256], 1: [3, 6, 12, 24, 48, 96, 192, 384],
       2: [5, 10, 20, 40, 80, 160, 320], 3: [7, 14, 28, 56, 112, 224], 4: [9, 18, 36, 72, 144, 288],
       5: [11, 22, 44, 88, 176, 352], 6: [13, 26, 52, 104, 208], 7: [15, 30, 60, 120, 240]}


#SNR = vector(RealField(10), [1, 1.5, 2, 2.5, 3, 3.5, 5, 4.5, 5, 5.5, 6])
#SNP = vector(RealField(4), [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45])
"""
[1/2, 7, 21, 'AWGN', 0.95, 0.4, True],[1/2, 6, 21, 'AWGN', 0.95, 0.4, True],[1/2, 5, 21, 'AWGN', 0.95, 0.4, True],
[1/2, 4, 21, 'AWGN', 0.95, 0.4, True],[1/2, 3, 21, 'AWGN', 0.95, 0.4, True],[1/2, 2, 21, 'AWGN', 0.95, 0.4, True],
[1/2, 1, 21, 'AWGN', 0.95, 0.4, True],

[1/2, 7, 53, 'AWGN', 0.95, 0.4, True],[1/2, 6, 53, 'AWGN', 0.95, 0.4, True],[1/2, 5, 53, 'AWGN', 0.95, 0.4, True],
[1/2, 4, 53, 'AWGN', 0.95, 0.4, True],[1/2, 3, 53, 'AWGN', 0.95, 0.4, True],[1/2, 2, 53, 'AWGN', 0.95, 0.4, True],
[1/2, 1, 53, 'AWGN', 0.95, 0.4, True],

[1/3, 7, 21, 'AWGN', 0.95, 0.4, True],[1/3, 6, 21, 'AWGN', 0.95, 0.4, True],[1/3, 5, 21, 'AWGN', 0.95, 0.4, True],
[1/3, 4, 21, 'AWGN', 0.95, 0.4, True],[1/3, 3, 21, 'AWGN', 0.95, 0.4, True],[1/3, 2, 21, 'AWGN', 0.95, 0.4, True],
[1/3, 1, 21, 'AWGN', 0.95, 0.4, True],

[1/3, 7, 53, 'AWGN', 0.95, 0.4, True],[1/3, 6, 53, 'AWGN', 0.95, 0.4, True],[1/3, 5, 53, 'AWGN', 0.95, 0.4, True],
[1/3, 4, 53, 'AWGN', 0.95, 0.4, True],[1/3, 3, 53, 'AWGN', 0.95, 0.4, True],[1/3, 2, 53, 'AWGN', 0.95, 0.4, True],
[1/3, 1, 53, 'AWGN', 0.95, 0.4, True],


[1/3, 0.1, 21, 'BEC', 1, 0, False],[1/3, 0.2, 21, 'BEC', 1, 0, False],[1/3, 0.3, 21, 'BEC', 1, 0, False],
[1/3, 0.4, 21, 'BEC', 1, 0, False],[1/3, 0.5, 21, 'BEC', 1, 0, False],[1/3, 0.6, 21, 'BEC', 1, 0, False],

[1/3, 0.1, 53, 'BEC', 1, 0, False], [1/3, 0.2, 53, 'BEC', 1, 0, False], [1/3, 0.3, 53, 'BEC', 1, 0, False],[1/3, 0.4, 53, 'BEC', 1, 0, False],
[1/3, 0.5, 53, 'BEC', 1, 0, False], [1/3, 0.6, 53, 'BEC', 1, 0, False],


[1/3, 0.02, 21, 'BSC', 1, 0, False], [1/3, 0.04, 21, 'BSC', 1, 0, False], [1/3, 0.06, 21, 'BSC', 1, 0, False],
[1/3, 0.08, 21, 'BSC', 1, 0, False], [1/3, 0.1, 21, 'BSC', 1, 0, False],
[1/3, 0.02, 21, 'BSC', 0.95, 0.4, False], [1/3, 0.04, 21, 'BSC', 0.95, 0.4, False], [1/3, 0.06, 21, 'BSC', 0.95, 0.4, False],
[1/3, 0.08, 21, 'BSC', 0.95, 0.4, False], [1/3, 0.1, 21, 'BSC', 0.95, 0.4, False],

[1/3, 0.02, 53, 'BSC', 1, 0, False], [1/3, 0.04, 53, 'BSC', 1, 0, False], [1/3, 0.06, 53, 'BSC', 1, 0, False],
[1/3, 0.08, 53, 'BSC', 1, 0, False], [1/3, 0.1, 53, 'BSC', 1, 0, False],
[1/3, 0.02, 53, 'BSC', 0.95, 0.4, False], [1/3, 0.04, 53, 'BSC', 0.95, 0.4, False], [1/3, 0.06, 53, 'BSC', 0.95, 0.4, False],
[1/3, 0.08, 53, 'BSC', 0.95, 0.4, False], [1/3, 0.1, 53, 'BSC', 0.95, 0.4, False],

[1/2, 0.02, 21, 'BSC', 1, 0, False], [1/2, 0.04, 21, 'BSC', 1, 0, False], [1/2, 0.06, 21, 'BSC', 1, 0, False],
[1/2, 0.08, 21, 'BSC', 1, 0, False], [1/2, 0.1, 21, 'BSC', 1, 0, False],
[1/2, 0.02, 21, 'BSC', 0.95, 0.4, False], [1/2, 0.04, 21, 'BSC', 0.95, 0.4, False], [1/2, 0.06, 21, 'BSC', 0.95, 0.4, False],
[1/2, 0.08, 21, 'BSC', 0.95, 0.4, False], [1/2, 0.1, 21, 'BSC', 0.95, 0.4, False],

[1/2, 0.02, 53, 'BSC', 1, 0, False], [1/2, 0.04, 53, 'BSC', 1, 0, False], [1/2, 0.06, 53, 'BSC', 1, 0, False],
[1/2, 0.08, 53, 'BSC', 1, 0, False], [1/2, 0.1, 53, 'BSC', 1, 0, False],
[1/2, 0.02, 53, 'BSC', 0.95, 0.4, False], [1/2, 0.04, 53, 'BSC', 0.95, 0.4, False], [1/2, 0.06, 53, 'BSC', 0.95, 0.4, False],
[1/2, 0.08, 53, 'BSC', 0.95, 0.4, False], [1/2, 0.1, 53, 'BSC', 0.95, 0.4, False],


[1/2, 0.02, 21, 'BSC', 0.95, 0.4, True], [1/2, 0.04, 21, 'BSC', 0.95, 0.4, True], [1/2, 0.06, 21, 'BSC', 0.95, 0.4, True],
[1/2, 0.08, 21, 'BSC', 0.95, 0.4, True], [1/2, 0.1, 21, 'BSC', 0.95, 0.4, True],

[1/2, 0.02, 53, 'BSC', 0.95, 0.4, True], [1/2, 0.04, 53, 'BSC', 0.95, 0.4, True], [1/2, 0.06, 53, 'BSC', 0.95, 0.4, True],
[1/2, 0.08, 53, 'BSC', 0.95, 0.4, True], [1/2, 0.1, 53, 'BSC', 0.95, 0.4, True],

[1/3, 0.02, 21, 'BSC', 0.95, 0.4, True], [1/3, 0.04, 21, 'BSC', 0.95, 0.4, True], [1/3, 0.06, 21, 'BSC', 0.95, 0.4, True],
[1/3, 0.08, 21, 'BSC', 0.95, 0.4, True], [1/3, 0.1, 21, 'BSC', 0.95, 0.4, True],

[1/3, 0.02, 53, 'BSC', 0.95, 0.4, True], [1/3, 0.04, 53, 'BSC', 0.95, 0.4, True], [1/3, 0.06, 53, 'BSC', 0.95, 0.4, True],
[1/3, 0.08, 53, 'BSC', 0.95, 0.4, True], [1/3, 0.1, 53, 'BSC', 0.95, 0.4, True],
[1/2, 7, 21, 'AWGN', 0.95, 0.4, True],[1/2, 6, 21, 'AWGN', 0.95, 0.4, True],[1/2, 5, 21, 'AWGN', 0.95, 0.4, True],
[1/2, 4, 21, 'AWGN', 0.95, 0.4, True],[1/2, 3, 21, 'AWGN', 0.95, 0.4, True],[1/2, 2, 21, 'AWGN', 0.95, 0.4, True],
[1/2, 1, 21, 'AWGN', 0.95, 0.4, True],

[1/2, 7, 53, 'AWGN', 0.95, 0.4, True],[1/2, 6, 53, 'AWGN', 0.95, 0.4, True],[1/2, 5, 53, 'AWGN', 0.95, 0.4, True],
[1/2, 4, 53, 'AWGN', 0.95, 0.4, True],[1/2, 3, 53, 'AWGN', 0.95, 0.4, True],[1/2, 2, 53, 'AWGN', 0.95, 0.4, True],
[1/2, 1, 53, 'AWGN', 0.95, 0.4, True],

[1/3, 7, 21, 'AWGN', 0.95, 0.4, True],[1/3, 6, 21, 'AWGN', 0.95, 0.4, True],[1/3, 5, 21, 'AWGN', 0.95, 0.4, True],
[1/3, 4, 21, 'AWGN', 0.95, 0.4, True],[1/3, 3, 21, 'AWGN', 0.95, 0.4, True],[1/3, 2, 21, 'AWGN', 0.95, 0.4, True],
[1/3, 1, 21, 'AWGN', 0.95, 0.4, True],

[1/3, 7, 53, 'AWGN', 0.95, 0.4, True],[1/3, 6, 53, 'AWGN', 0.95, 0.4, True],[1/3, 5, 53, 'AWGN', 0.95, 0.4, True],
[1/3, 4, 53, 'AWGN', 0.95, 0.4, True],[1/3, 3, 53, 'AWGN', 0.95, 0.4, True],[1/3, 2, 53, 'AWGN', 0.95, 0.4, True],
[1/3, 1, 53, 'AWGN', 0.95, 0.4, True],


[1/3, 0.1, 21, 'BEC', 1, 0, False],[1/3, 0.2, 21, 'BEC', 1, 0, False],[1/3, 0.3, 21, 'BEC', 1, 0, False],
[1/3, 0.4, 21, 'BEC', 1, 0, False],[1/3, 0.5, 21, 'BEC', 1, 0, False],[1/3, 0.6, 21, 'BEC', 1, 0, False],

[1/3, 0.1, 53, 'BEC', 1, 0, False], [1/3, 0.2, 53, 'BEC', 1, 0, False], [1/3, 0.3, 53, 'BEC', 1, 0, False],[1/3, 0.4, 53, 'BEC', 1, 0, False],
[1/3, 0.5, 53, 'BEC', 1, 0, False], [1/3, 0.6, 53, 'BEC', 1, 0, False],


[1/3, 0.02, 21, 'BSC', 1, 0, False], [1/3, 0.04, 21, 'BSC', 1, 0, False], [1/3, 0.06, 21, 'BSC', 1, 0, False],
[1/3, 0.08, 21, 'BSC', 1, 0, False], [1/3, 0.1, 21, 'BSC', 1, 0, False],
[1/3, 0.02, 21, 'BSC', 0.95, 0.4, False], [1/3, 0.04, 21, 'BSC', 0.95, 0.4, False], [1/3, 0.06, 21, 'BSC', 0.95, 0.4, False],
[1/3, 0.08, 21, 'BSC', 0.95, 0.4, False], [1/3, 0.1, 21, 'BSC', 0.95, 0.4, False],

[1/3, 0.02, 53, 'BSC', 1, 0, False], [1/3, 0.04, 53, 'BSC', 1, 0, False], [1/3, 0.06, 53, 'BSC', 1, 0, False],
[1/3, 0.08, 53, 'BSC', 1, 0, False], [1/3, 0.1, 53, 'BSC', 1, 0, False],
[1/3, 0.02, 53, 'BSC', 0.95, 0.4, False], [1/3, 0.04, 53, 'BSC', 0.95, 0.4, False], [1/3, 0.06, 53, 'BSC', 0.95, 0.4, False],
[1/3, 0.08, 53, 'BSC', 0.95, 0.4, False], [1/3, 0.1, 53, 'BSC', 0.95, 0.4, False],


[1/2, 0.02, 21, 'BSC', 1, 0, False], [1/2, 0.04, 21, 'BSC', 1, 0, False], [1/2, 0.06, 21, 'BSC', 1, 0, False],
[1/2, 0.08, 21, 'BSC', 1, 0, False], [1/2, 0.1, 21, 'BSC', 1, 0, False],
[1/2, 0.02, 21, 'BSC', 0.95, 0.4, False], [1/2, 0.04, 21, 'BSC', 0.95, 0.4, False], [1/2, 0.06, 21, 'BSC', 0.95, 0.4, False],
[1/2, 0.08, 21, 'BSC', 0.95, 0.4, False], [1/2, 0.1, 21, 'BSC', 0.95, 0.4, False],

[1/2, 0.02, 53, 'BSC', 1, 0, False], [1/2, 0.04, 53, 'BSC', 1, 0, False], [1/2, 0.06, 53, 'BSC', 1, 0, False],
[1/2, 0.08, 53, 'BSC', 1, 0, False], [1/2, 0.1, 53, 'BSC', 1, 0, False],
[1/2, 0.02, 53, 'BSC', 0.95, 0.4, False], [1/2, 0.04, 53, 'BSC', 0.95, 0.4, False], [1/2, 0.06, 53, 'BSC', 0.95, 0.4, False],
[1/2, 0.08, 53, 'BSC', 0.95, 0.4, False], [1/2, 0.1, 53, 'BSC', 0.95, 0.4, False],


[1/2, 0.02, 21, 'BSC', 0.95, 0.4, True], [1/2, 0.04, 21, 'BSC', 0.95, 0.4, True], [1/2, 0.06, 21, 'BSC', 0.95, 0.4, True],
[1/2, 0.08, 21, 'BSC', 0.95, 0.4, True], [1/2, 0.1, 21, 'BSC', 0.95, 0.4, True],

[1/2, 0.02, 53, 'BSC', 0.95, 0.4, True], [1/2, 0.04, 53, 'BSC', 0.95, 0.4, True], [1/2, 0.06, 53, 'BSC', 0.95, 0.4, True],
[1/2, 0.08, 53, 'BSC', 0.95, 0.4, True], [1/2, 0.1, 53, 'BSC', 0.95, 0.4, True],

[1/3, 0.02, 21, 'BSC', 0.95, 0.4, True], [1/3, 0.04, 21, 'BSC', 0.95, 0.4, True], [1/3, 0.06, 21, 'BSC', 0.95, 0.4, True],
[1/3, 0.08, 21, 'BSC', 0.95, 0.4, True], [1/3, 0.1, 21, 'BSC', 0.95, 0.4, True],

[1/3, 0.02, 53, 'BSC', 0.95, 0.4, True], [1/3, 0.04, 53, 'BSC', 0.95, 0.4, True], [1/3, 0.06, 53, 'BSC', 0.95, 0.4, True],
[1/3, 0.08, 53, 'BSC', 0.95, 0.4, True], [1/3, 0.1, 53, 'BSC', 0.95, 0.4, True],
"""
runs = 10000
lim = 1000

runs_vals =[
[1/2, 2, 300, 'AWGN', 0.95, 0.4, False],
[1/2, 2, 400, 'AWGN', 0.95, 0.4, False],
[1/2, 2, 500, 'AWGN', 0.95, 0.4, False],

[1/2, 3, 300, 'AWGN', 0.95, 0.4, False],
[1/2, 3, 400, 'AWGN', 0.95, 0.4, False],
[1/2, 3, 500, 'AWGN', 0.95, 0.4, False],

[1/2, 4, 300, 'AWGN', 0.95, 0.4, False],
[1/2, 4, 400, 'AWGN', 0.95, 0.4, False],
[1/2, 4, 500, 'AWGN', 0.95, 0.4, False],
]

def non_zero_matrix(input_matrix):
    output = []
    for i in range(input_matrix.nrows()):
        temp = {}
        for j in input_matrix.row(i).nonzero_positions():
            temp.update({j: 1})
        output.append(temp)
    return output


for elem in runs_vals:
    print(elem)
    rate = elem[0]
    snr = elem[1]
    A = elem[2]
    channel = elem[3]
    gamma, lam = elem[4], elem[5]
    use_core = elem[6]

    N0, sig = None, None
    pol = CRC.get_pol(A)
    B = A + pol.degree()
    bg = PF.det_BG(A, rate)
    L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
    K_ap = B_ap // C
    Kb = PF.determine_kb(B=B, bg=bg)
    Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss, K_ap=K_ap)
    BG, Bi = BG_main.create_BG(Zc, iLS, bg)
    H = HF.Protograph(BG, Zc)

    del BG;
    gc.collect()


    if channel == 'AWGN':
        sig = sqrt(1 / (2 * rate * 10 ** (snr / 10)))
        p = 1 - norm.cdf(1 / sig)  # error probability, from proposition 2.9
        N0 = 2*sig** 2
        print(f"sigma:{sig}, N0:{N0}, pross:{p}")
    BLER, BER, FAR, AVGit = 0, 0, 0, 0
    start_time = time.time()
    for iterations in range(runs):
        if BLER >= lim:
            break
        if iterations % 100 == 0:
            print(time.time() - start_time)
            start_time = time.time()
            print(f"BLER:{BLER}, BER:{BER}")
            print(iterations)
            a = random_vector(GF(2), A)
            #print(f'a:{a}')
            c, G = CRC.CRC(a, A, pol)
            crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=c)   # crk := padding the codeword
            D = vector(GF(2), crk)
            u = LDPC_Encoding.Encoding(H=H, Bi=Bi, Zc=Zc, D=D, K=K, kb=Kb, BG=bg)
            e, HRM = LDPC_Rate_Matching.RM_main(u=u, Zc=Zc, H=H, K=K, K_ap=K_ap, rate=rate, B=B, channel=channel)
            HNZ = non_zero_matrix(HRM)
        r = HF.channel_noise(s=e, channel=channel, p=sig if channel == 'AWGN' else snr)
        llr_r = LDPC_Rate_Matching.fill_w_llr(r=r, Zc=Zc, K=K, K_ap=K_ap, p=snr, N0=N0, channel=channel, HRM=HRM)
        if channel == 'BEC':
            import minsum_BEC
            aa, suces, iter = minsum_BEC.minsum_BEC(HRM, llr_r)
            #print(suces, iter)
        else:
            import LDPC_MinSum
            aa, suces, iter = LDPC_MinSum.minsum_SPA(H=HRM, HNZ=HNZ, llr=llr_r, rcore=4 * Zc, lam=lam,gamma=gamma, Zc=Zc,K=K, N0=N0, use_core=use_core)
        crc_check = CRC.CRC_check(aa[:B], len(aa[:B]), pol)
        BER += (aa[:A]+a).hamming_weight()
        BLER += sign((aa[:A]+a).hamming_weight())
        FAR += sign((aa[:A]+a).hamming_weight()) and not sign(crc_check.hamming_weight())
        AVGit += iter

        #print(iter)

    with open(f'Tests/{channel}.csv', mode='a',
              newline='') as file:
        result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        result_writer.writerow(
            [A, rate, B, HRM.ncols(), iterations, BER, BLER, snr, AVGit, f'OMS,gamma:{gamma},lam{lam},AMS:{use_core}', datetime.datetime.now()])
        gc.collect()

