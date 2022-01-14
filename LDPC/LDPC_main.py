import sys
import CRC
import Parameter_Functions as PF
import LDPC_Encoding
import LDPC_Rate_Matching
from sage.all import *
# cd Desktop/Thesis/PySageMath/LDPC

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

# sage LDPC_main.py 20 1/2
if __name__ == "__main__":
    a, A = list(random_vector(GF(2), int(sys.argv[1]))), int(sys.argv[1])
    a = vector(GF(2), [1]*A)
    print(f"a := {a}")
    R = [int(x) for x in sys.argv[2].split('/')]
    R = R[0] / R[1]

    bg = PF.det_BG(A, R)
    b = CRC.main_CRC(a, crc24a)
    B = len(b)
    L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
    K_ap = B_ap // C
    Kb = PF.determine_kb(B=B, bg=bg)
    Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss,K_ap=K_ap)
    crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=b)
    _, D = PF.get_d_c(Zc=Zc, K=K, C=crk)
    D = vector(GF(2), D)
    X, H, BG = LDPC_Encoding.Encoding(bg=bg, iLS=iLS, Zc=Zc, D=D, K=K, kb=Kb)
    E = LDPC_Rate_Matching.ty(D=X, Zc=Zc, BG=BG, H=H, K=K, K_ap=K_ap, B=B, R=R)
    E = LDPC_Rate_Matching.RM_main(D=X, Zc=Zc, BG=BG, H=H, K=K, K_ap=K_ap, B=B, R=R)
    breakpoint()
