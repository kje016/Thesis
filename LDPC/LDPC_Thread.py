# cd Desktop/Thesis/PySageMath/LDPC
import time
import threading

from sage.all import *

import gc
import csv
import datetime
from scipy.stats import norm

import CRC
import Parameter_Functions as PF
import LDPC_Encoding
import LDPC_Rate_Matching
import LDPC_HelperFunctions as HF


lss = {0: [2, 4, 8, 16, 32, 64, 128, 256], 1: [3, 6, 12, 24, 48, 96, 192, 384],
       2: [5, 10, 20, 40, 80, 160, 320], 3: [7, 14, 28, 56, 112, 224], 4: [9, 18, 36, 72, 144, 288],
       5: [11, 22, 44, 88, 176, 352], 6: [13, 26, 52, 104, 208], 7: [15, 30, 60, 120, 240]}

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

runs = 1000

""""[20, 0.1, 1/2, 'BSC'],  [20, 0.09, 1/2, 'BSC'], [20, 0.08, 1/2, 'BSC'], [20, 0.07, 1/2, 'BSC'],
              [20, 0.06, 1/2, 'BSC'], [20, 0.05, 1/2, 'BSC'], [20, 0.04, 1/2, 'BSC'], [20, 0.03, 1/2, 'BSC'],
              [20, 0.02, 1/2, 'BSC'], 
"""

def thread_decoder(e, channel, sig, snr, Zc, K, K_ap, HRM, B, a):
    r = HF.channel_noise(s=e, channel=channel, p=sig if channel == 'AWGN' else snr)
    llr_r = LDPC_Rate_Matching.fill_w_llr(r=r, Zc=Zc, K=K, K_ap=K_ap, p=N0 if channel == 'AWGN' else snr,
                                          channel=channel)
    if channel == 'BEC':
        import minsum_BEC
        aa, suces, iters = minsum_BEC.minsum_BEC(HRM, llr_r)
    else:
        import LDPC_MinSum
        aa, suces, iters = LDPC_MinSum.minsum_SPA(HRM, llr_r, channel, sig, 4 * Zc)
    crc_check = CRC.CRC_check(aa[:B], len(aa[:B]), pol)

    BER[0] = BER[0] + (aa[:A] + a).hamming_weight()
    BLER[0] = BLER[0] + sign(crc_check.hamming_weight())
    FAR[0] = FAR[0] + sign((aa[:A] + a).hamming_weight()) and not sign(crc_check.hamming_weight())
    AVGit[0] = AVGit[0] + iters


runs_vals =[['', 1], ['', 2], ['', 3], ['', 4], ['', 5]]
A = 21
rate = 2/5
channel = 'AWGN'
pol = CRC.get_pol(A)
B = A + pol.degree()

bg = PF.det_BG(A, rate)
L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
K_ap = B_ap // C
Kb = PF.determine_kb(B=B, bg=bg)
Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss, K_ap=K_ap)
BG = HF.get_base_matrix(bg, iLS, Zc)
#BGB = BG.matrix_from_rows_and_columns(list(range(4)), list(range(10, 10+4)))
H = HF.Protograph(BG, Zc)


# sage LDPC_main.py 20 bsc
for elem in runs_vals:
    snr = elem[1]

    if channel == 'AWGN':
        sig = sqrt(1 / (2 * rate * 10 ** (snr / 10)))
        p = 1 - norm.cdf(1 / sig)  # error probability, from proposition 2.9
        N0 = 2 * sig ** 2
        print(f"sigma:{sig}, N0:{N0}, pross:{p}")
    BLER, BER, FAR, AVGit = [0], [0], [0], [0]
    start_time = time.time()
    a = random_vector(GF(2), A)
    c, G = CRC.CRC(a, A, pol)

    crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=c)  # crk := padding the codeword
    D = vector(GF(2), crk)
    u = LDPC_Encoding.Encoding(H=H, Zc=Zc, D=D, K=K, kb=Kb, BG=bg)
    for iterations in range(runs):
        if runs % 100 == 0:
            print(time.time() - start_time)
            start_time = time.time()
            print(f"BLER:{BLER}, BER:{BER}")
            print(iterations)
            a = random_vector(GF(2), A)
            c, G = CRC.CRC(a, A, pol)

            crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=c)  # crk := padding the codeword
            D = vector(GF(2), crk)
            u = LDPC_Encoding.Encoding(H=H, Zc=Zc, D=D, K=K, kb=Kb, BG=bg)
            e, HRM = LDPC_Rate_Matching.RM_main(u=u, Zc=Zc, H=H, K=K, K_ap=K_ap, rate=rate, B=B, channel=channel)

        thread1 = threading.Thread(target=thread_decoder, args=(e, channel, sig, snr, Zc, K, K_ap, HRM, B, a))
        thread2 = threading.Thread(target=thread_decoder, args=(e, channel, sig, snr, Zc, K, K_ap, HRM, B, a))
        thread3 = threading.Thread(target=thread_decoder, args=(e, channel, sig, snr, Zc, K, K_ap, HRM, B, a))

        thread1.start()
        thread2.start()
        thread3.start()

        thread1.join()
        thread2.join()
        thread3.join()

    with open(f'PycharmProjects\\Thesis\\PySageMath\\LDPC\\Tests\\{channel}.csv', mode='a',
              newline='') as file:
        result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        result_writer.writerow(
            [A, rate, K, len(e), runs, BER, BLER, snr, AVGit, 'threads=3', datetime.datetime.now()])
        gc.collect()
    """
    with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\LDPC\\Tests\\{channel}.csv', mode='a',
              newline='') as file:
        result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        result_writer.writerow(
            [A, rate, K, H.ncols(), runs, BER, BLER, snr, iter, suces, datetime.datetime.now()])
        gc.collect()
    """

