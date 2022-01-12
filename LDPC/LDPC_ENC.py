# cd Desktop/Thesis/PySageMath/LDPC
# personal access token: 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
from sage.all import *
import LDPC_HelperFunctions as HF
import Parameter_Functions as PF


lss = {0: [2, 4, 8, 16, 32, 64, 128, 256], 1: [3, 6, 12, 24, 48, 96, 192, 384],
       2: [5, 10, 20, 40, 80, 160, 320], 3: [7, 14, 28, 56, 112, 224], 4: [9, 18, 36, 72, 144, 288],
       5: [11, 22, 44, 88, 176, 352], 6: [13, 26, 52, 104, 208], 7: [15, 30, 60, 120, 240]}


def mul_sh(a, vec):
    if a == -1:
        return vector(GF(2), [0]*len(vec))
    res = vec[a:]+ vec[:a]
    return vector(GF(2), res)


def calc_lambdas(kb, H, Z, D, K):
    result, A = [], H.matrix_from_rows_and_columns(list(range(4 * Z)), list(range(K)))
    for j in range(4):
        lam_j = vector(GF(2), Z)
        for l in range(kb):
            a = H.matrix_from_rows_and_columns(list(range((j*4), (j+1)*4)), list(range(l*Z, (l+1)*Z)))
            lam_j = lam_j + (a*D[l*Z:(l+1)*Z])
        result.append(lam_j)
    return result


B, bg = 20, 2

b_bits = [1,1,1,0, 0,1,0,1, 1,0,0,1, 0,1,1,0, 1,0,1,0]
L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
K_ap = B_ap //C

kb = PF.determine_kb(B=B_ap, bg=bg)
Zc, iLS, K = PF.det_Z(bg, kb, lss, K_ap)
crk = HF.calc_crk(C=C, K_ap=K_ap, K=K, L=L, b_bits=b_bits)
_, D = HF.get_d_c(Zc=Zc, K=K, C=crk)
D = vector(GF(2), D)

BG = HF.get_base_matrix(bg, iLS, Zc)
H = HF.Protograph(BG, Zc)
mb, nb = BG.nrows(), BG.ncols()
Am, Bm = H.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(10 * Zc))), H.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range((10 * Zc), (10 * Zc) + 4 * Zc)))
Cm, Dm = H.matrix_from_rows_and_columns(list(range((mb-4) * Zc)), list(range(kb * Zc))), H.matrix_from_rows_and_columns(list(range((mb - 4) * Zc)), list(range(kb * Zc, (kb * Zc) + 4 * Zc)))
BGA, BGB = BG.matrix_from_rows_and_columns(list(range(Zc)), list(range(kb+Zc))), BG.matrix_from_rows_and_columns(list(range(Zc)), list(range(kb+Zc, kb+Zc+Zc)))

# Pc (core parity): can be calculated from submatrices A & B
lambdas = calc_lambdas(kb, H, Zc, D, K)
pc1_shift = sum(lambdas)
pc1 = vector(GF(2), [pc1_shift[-1]] + list(pc1_shift)[:len(pc1_shift)-1])   # TODO: find correct shift
pc2 = lambdas[0] + pc1
pc4 = lambdas[3] + pc1
pc3 = lambdas[1] + pc2
"""
pc4 = lambdas[3] + pc1_shift
pc3 = lambdas[1] + pc2
"""
"""
pc3 = lambdas[1] + pc2
pc4 = lambdas[2] + pc1_shift
"""
"""
pc4 = lambdas[3] + pc1_shift
pc3 = lambdas[2] + pc4
"""
"""
pc3 = lambdas[1] + pc2
pc4 = lambdas[3] + pc1_shift
"""
"""
pc3 = lambdas[1] + pc2
pc4 = lambdas[3] + pc3
"""
# Pa (additional parity): calculated from information core and parit bits using C & D
def calc_pa(H, Pc, D):
    CD = H.matrix_from_rows_and_columns(list(range(4*Zc, mb*Zc)), list(range((10 * Zc) + (4 * Zc))))
    return CD * vector(GF(2), list(D)+list(Pc))


I = D  # len of kb*Z
Pc = vector(GF(2), list(pc1)+list(pc2)+list(pc3)+list(pc4))    # len of 4*Z
Pa = calc_pa(H=H, Pc=Pc, D=D)     # len of (mb-4)*Z
# x = [i pc pa]
X = vector(GF(2), list(I)+list(Pc)+list(Pa))
print(f"H*X := {H*X}")
print(f"pc1 := {pc1}")
print(f"pc1_shift {pc1_shift}")
print()
print(f"Am*D := {Am*D}")
print(f"Bm*Pc:= {Bm*Pc}")
print(f"res  :={Am*D+Bm*Pc}")
print(f"H*X ==0 := {H*X==0}")

