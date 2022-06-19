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
[3/10, 6, 53, 'AWGN', 1, 0],[3/10, 5, 53, 'AWGN', 1, 0],[3/10, 4, 53, 'AWGN', 1, 0],
[3/10, 3, 53, 'AWGN', 1, 0],[3/10, 2, 53, 'AWGN', 1, 0],[3/10, 1, 53, 'AWGN', 1, 0],

[3/10, 6, 53, 'AWGN', 1, 0.2],[3/10, 5, 53, 'AWGN', 1, 0.2],[3/10, 4, 53, 'AWGN', 1, 0.2],
[3/10, 3, 53, 'AWGN', 1, 0.2],[3/10, 2, 53, 'AWGN', 1, 0.2],[3/10, 1, 53, 'AWGN', 1, 0.2],

[3/10, 6, 53, 'AWGN', 0.95, 0.4],[3/10, 5, 53, 'AWGN', 0.95, 0.4],[3/10, 4, 53, 'AWGN', 0.95, 0.4],
[3/10, 3, 53, 'AWGN', 0.95, 0.4],[3/10, 2, 53, 'AWGN', 0.95, 0.4],[3/10, 1, 53, 'AWGN', 0.95, 0.4],

[3/10, 6, 53, 'AWGN', 0.95, 0.2],[3/10, 5, 53, 'AWGN', 0.95, 0.2],[3/10, 4, 53, 'AWGN', 0.95, 0.2],
[3/10, 3, 53, 'AWGN', 0.95, 0.2],[3/10, 2, 53, 'AWGN', 0.95, 0.2],[3/10, 1, 53, 'AWGN', 0.95, 0.2],

[3/10, 6, 53, 'AWGN', 1, 0],[3/10, 5, 53, 'AWGN', 1, 0],[3/10, 4, 53, 'AWGN', 1, 0],
[3/10, 3, 53, 'AWGN', 1, 0],[3/10, 2, 53, 'AWGN', 1, 0],[3/10, 1, 53, 'AWGN', 1, 0],

[3/10, 6, 53, 'AWGN', 1, 0.2],[3/10, 5, 53, 'AWGN', 1, 0.2],[3/10, 4, 53, 'AWGN', 1, 0.2],
[3/10, 3, 53, 'AWGN', 1, 0.2],[3/10, 2, 53, 'AWGN', 1, 0.2],[3/10, 1, 53, 'AWGN', 1, 0.2],

[3/10, 6, 53, 'AWGN', 0.95, 0.4],[3/10, 5, 53, 'AWGN', 0.95, 0.4],[3/10, 4, 53, 'AWGN', 0.95, 0.4],
[3/10, 3, 53, 'AWGN', 0.95, 0.4],[3/10, 2, 53, 'AWGN', 0.95, 0.4],[3/10, 1, 53, 'AWGN', 0.95, 0.4],

[3/10, 6, 53, 'AWGN', 0.95, 0.2],[3/10, 5, 53, 'AWGN', 0.95, 0.2],[3/10, 4, 53, 'AWGN', 0.95, 0.2],
[3/10, 3, 53, 'AWGN', 0.95, 0.2],[3/10, 2, 53, 'AWGN', 0.95, 0.2],[3/10, 1, 53, 'AWGN', 0.95, 0.2],

[1/2, 6, 53, 'AWGN', 0.95, 0.4],[1/2, 5, 53, 'AWGN', 0.95, 0.4],[1/2, 4, 53, 'AWGN', 0.95, 0.4],
[1/2, 3, 53, 'AWGN', 0.95, 0.4],[1/2, 2, 53, 'AWGN', 0.95, 0.4],[1/2, 1, 53, 'AWGN', 0.95, 0.4],

[1/2, 6, 53, 'AWGN', 0.95, 0.2],[1/2, 5, 53, 'AWGN', 0.95, 0.2],[1/2, 4, 53, 'AWGN', 0.95, 0.2],
[1/2, 3, 53, 'AWGN', 0.95, 0.2],[1/2, 2, 53, 'AWGN', 0.95, 0.2],[1/2, 1, 53, 'AWGN', 0.95, 0.2],
"""
runs = 10000
lim = 1000

runs_vals =[
[1/3, 7, 21, 'AWGN', 1, 0],[1/3, 6, 21, 'AWGN', 1, 0],[1/3, 5, 21, 'AWGN', 1, 0],[1/3, 4, 21, 'AWGN', 1, 0],
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
    #breakpoint()
    for iterations in range(runs):
        if BLER >= lim:
            break
        if iterations % 50 == 0:
            print(f"BLER:{BLER}, BER:{BER}, FAR:{FAR}")
        if iterations % 100 == 0:
            print(iterations)
            print(time.time() - start_time)
            start_time = time.time()
            #print(f"BLER:{BLER}, BER:{BER}")
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
            aa, suces, iter = LDPC_MinSum.minsum_SPA(H=HRM, HNZ=HNZ, llr=llr_r, rcore=4 * Zc, lam=lam,gamma=gamma, Zc=Zc, K=K, r=r, N0=N0)
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
            [A, rate, B, HRM.ncols(), iterations, BER, BLER, snr, AVGit, f'OMS, gamma:{gamma}, :lam{lam}', datetime.datetime.now()])
        gc.collect()

