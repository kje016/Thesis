# cd Desktop/Thesis/PySageMath/LDPC
import sys
import CRC
import LDPC_Decoding
import Parameter_Functions as PF
import LDPC_Encoding
import LDPC_Rate_Matching
import LDPC_HelperFunctions as HF
import minsum_BEC
from scipy.stats import norm

from sage.all import *

lss = {0: [2, 4, 8, 16, 32, 64, 128, 256], 1: [3, 6, 12, 24, 48, 96, 192, 384],
       2: [5, 10, 20, 40, 80, 160, 320], 3: [7, 14, 28, 56, 112, 224], 4: [9, 18, 36, 72, 144, 288],
       5: [11, 22, 44, 88, 176, 352], 6: [13, 26, 52, 104, 208], 7: [15, 30, 60, 120, 240]}

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24a = x**24 + x**23 + x**18 + x**17 + x**14 + x**11 + x**10 + x**7 + x**6 + x**5 + x**4 + x**3 + x + x**0
crc24b = x**24 + x**23 + x**6 + x**5 + x + x**0
crc24c = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0

SNR = vector(RealField(10), [1, 1.5, 2, 2.5, 3, 3.5, 5, 4.5, 5, 5.5, 6])
# sage LDPC_main.py 20 1/2 AWGN
if __name__ == "__main__":
    runs = 40
    for i in range(runs):
        a, A = list(random_vector(GF(2), int(sys.argv[1]))), int(sys.argv[1])
        # a = vector(GF(2), [1]*A)
        print(f"a := {a}")
        R = [int(x) for x in sys.argv[2].split('/')]
        R = R[0] / R[1]
        channel = sys.argv[3].upper()

        sigma = vector(RealField(10), map(lambda z: sqrt(1/(2*R*10**(z/10))), SNR))
        N0 = 2*sigma[0]**2

        bg = PF.det_BG(A, R)


        b = CRC.main_CRC(a, crc24a)
        B = len(b)
        L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
        K_ap = B_ap // C
        Kb = PF.determine_kb(B=B, bg=bg)
        Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss, K_ap=K_ap)
        # print(f"Zc := {Zc}")
        crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=b)
        _, D = PF.get_d_c(Zc=Zc, K=K, C=crk)
        D = vector(GF(2), D)

        X, H, BG = LDPC_Encoding.Encoding(bg=bg, iLS=iLS, Zc=Zc, D=D, K=K, kb=Kb)
        e, HRM = LDPC_Rate_Matching.RM_main(D=X, Zc=Zc, H=H, K=K, K_ap=K_ap, R=R)
        r = HF.channel_noise(e, channel, 0.1)
        # if 'AWGN' -> channel_noise(e, 'AWGN', sigma)
        # if 'BSC' || 'BSC' -> channel_noise(e, 'BSC'/'BSC', cross_p)
        llr_r = LDPC_Rate_Matching.fill_e(r, Zc, K, K_ap, 0.1, H.ncols()-H.nrows(), channel)
        if channel == 'BEC':
            aa, is_codeword = minsum_BEC.minsum_BEC(HRM, llr_r)
        else:
            aa, is_codeword = LDPC_Decoding.spa_main(HRM, llr_r, N0, channel, 0.1)
        #print(f"H*v_hat := {is_codeword}")
        if is_codeword:
            crc_check = CRC.CRC_check(aa[:B], crc24a)
            print(f"crc_check := {vector(GF(2), crc_check) == 0}")
            print(f" aa = {aa[:A]}")
        print(f"!!         RUNS:= {i}          !!")