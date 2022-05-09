# cd Desktop/Thesis/PySageMath/LDPC
import gc
import sys
import csv
import datetime

import CRC
import LDPC_MinSum
import Parameter_Functions as PF
import LDPC_Encoding
import LDPC_Rate_Matching
import LDPC_HelperFunctions as HF
import minsum_BEC
import OMS


from sage.all import *

lss = {0: [2, 4, 8, 16, 32, 64, 128, 256], 1: [3, 6, 12, 24, 48, 96, 192, 384],
       2: [5, 10, 20, 40, 80, 160, 320], 3: [7, 14, 28, 56, 112, 224], 4: [9, 18, 36, 72, 144, 288],
       5: [11, 22, 44, 88, 176, 352], 6: [13, 26, 52, 104, 208], 7: [15, 30, 60, 120, 240]}

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()


SNR = vector(RealField(10), [1, 1.5, 2, 2.5, 3, 3.5, 5, 4.5, 5, 5.5, 6])
SNP = vector(RealField(4), [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45])
R = [1/2] # [1/2, 2/5, 1/3, 1/4,  1/5]   # Rate of the code
runs = 100
# sage LDPC_main.py 20 1/2 bsc
if __name__ == "__main__":
    A = int(sys.argv[1])
    channel = sys.argv[2].upper()
    pol = CRC.get_pol(A)
    B = A + pol.degree()

    #sigma = vector(RealField(10), map(lambda z: sqrt(1 / (2 * R * 10 ** (z / 10))), SNR))
    #N0 = 2 * sigma[0] ** 2

    bg = PF.det_BG(A, R)
    B = A + crc24a.degree()
    L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
    K_ap = B_ap // C
    Kb = PF.determine_kb(B=B, bg=bg)
    Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss, K_ap=K_ap)
    BG = HF.get_base_matrix(bg, iLS, Zc)
    H = HF.Protograph(BG, Zc)
    del BG; gc.collect()

    for snr in SNR:
        if channel == 'AWNG':
            sigma = sqrt(1 / (2 * rate * 10 ** (snr / 10)))
            p = 1 - norm.cdf(1 / sigma)  # error probability, from proposition 2.9
            N0 = 2 * sigma ** 2
        BLER, BER = 0, 0

        for iteration in range(runs):
            if iteration % 100 == 0:
                print(iteration)

            a = list(random_vector(GF(2), int(sys.argv[1])))
            c, G = CRC.CRC(a, A, pol)
            # crk := padding the codeword
            crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=b)   # TODO: testing for C > 1 & need to split crk
            D = vector(GF(2), crk)
            u = LDPC_Encoding.Encoding(H=H, Zc=Zc, D=D, K=K, kb=Kb, BG=bg)
            e, HRM = LDPC_Rate_Matching.RM_main(u=u, Zc=Zc, H=H, K=K, K_ap=K_ap, R=R)

            r = HF.channel_noise(e, channel, snr)
            # if 'AWGN' -> channel_noise(e, 'AWGN', sigma)
            # if 'BSC' || 'BSC' -> channel_noise(e, 'BSC'/'BSC', cross_p)
            llr_r = LDPC_Rate_Matching.fill_w_llr(r, Zc, K, K_ap, snr, H.ncols() - H.nrows(), channel)
            # tess = OMS.OMS(Zc=Zc, H=HRM, r=llr_r)
            if channel == 'BEC':
                aa, is_codeword = minsum_BEC.minsum_BEC(HRM, llr_r)
            else:
                aa, is_codeword = LDPC_MinSum.minsum_SPA(HRM, llr_r, N0, channel, snr, 4*Zc)
            if is_codeword:
                crc_check = CRC.CRC_check(aa[:B], len(aa[:B]), pol)

        file_getter = channel + '_' + decoder
        with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\{file_getter}.csv', mode='a',
                  newline='') as file:
            result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            result_writer.writerow(
                [A, R, K, H.ncols(), runs, BER, BLER, snr, datetime.datetime.now()])
            gc.collect()