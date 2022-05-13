from sage.all import *
import csv
import CRC
import LDPC_HelperFunctions as HF
import LDPC_Encoding
import LDPC_Rate_Matching
import Parameter_Functions as PF


lss = {0: [2, 4, 8, 16, 32, 64, 128, 256], 1: [3, 6, 12, 24, 48, 96, 192, 384],
       2: [5, 10, 20, 40, 80, 160, 320], 3: [7, 14, 28, 56, 112, 224], 4: [9, 18, 36, 72, 144, 288],
       5: [11, 22, 44, 88, 176, 352], 6: [13, 26, 52, 104, 208], 7: [15, 30, 60, 120, 240]}

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

rate = 1 / 3
for iteration in range(20, 1000, 10):
    A = iteration
    pol = CRC.get_pol(A)
    B = A + pol.degree()

    # sigma = vector(RealField(10), map(lambda z: sqrt(1 / (2 * R * 10 ** (z / 10))), SNR))
    # N0 = 2 * sigma[0] ** 2
    bg = PF.det_BG(A, rate)
    L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
    K_ap = B_ap // C
    Kb = PF.determine_kb(B=B, bg=bg)
    Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss, K_ap=K_ap)
    BG = HF.get_base_matrix(bg, iLS, Zc)
    BGB = BG.matrix_from_rows_and_columns(list(range(4)), list(range(10, 10 + 4)))
    print(bg, iLS, Zc)
    H = HF.Protograph(BG, Zc)
    a = random_vector(GF(2), A)
    c, G = CRC.CRC(a, A, pol)
    # crk := padding the codeword
    crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=c)  # TODO: testing for C > 1 & need to split crk
    D = vector(GF(2), crk)
    u = LDPC_Encoding.Encoding(H=H, Zc=Zc, D=D, K=K, kb=Kb, BG=bg)
    e, HRM = LDPC_Rate_Matching.RM_main(u=u, Zc=Zc, H=H, K=K, K_ap=K_ap, rate=rate, B=B)

