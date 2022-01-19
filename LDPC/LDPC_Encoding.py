# cd Desktop/Thesis/PySageMath/LDPC
# personal access token: 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
from sage.all import *
import LDPC_HelperFunctions as HF


def mul_sh(a, vec):
    if a == -1:
        return vector(GF(2), [0]*len(vec))
    res = vec[a:] + vec[:a]
    return vector(GF(2), res)


def calc_lambdas(kb, H, Z, D, K):
    result, A = [], H.matrix_from_rows_and_columns(list(range(4 * Z)), list(range(K)))
    for j in range(4):
        lam_j = vector(GF(2), Z)
        for l in range(kb):
            a = H.matrix_from_rows_and_columns(list(range((j*Z), (j+1)*Z)), list(range(l*Z, (l+1)*Z)))
            lam_j = lam_j + (a*D[l*Z:(l+1)*Z])
        result.append(lam_j)
    return result


# Pa (additional parity): calculated from information core and parity bits using C & D
def calc_pa(H, Pc, D, Zc, K):
    CD = H.matrix_from_rows_and_columns(list(range(4*Zc, H.nrows())), list(range(K + (4 * Zc))))
    return CD * vector(GF(2), list(D)+list(Pc))


def Encoding(bg, iLS, Zc, D, K, kb):
    BG = HF.get_base_matrix(bg, iLS, Zc)
    H = HF.Protograph(BG, Zc)
    # Pc (core parity): can be calculated from submatrices A & B
    lambdas = calc_lambdas(kb, H, Zc, D, K)
    pc1_shift = sum(lambdas)
    pc1 = vector(GF(2), [pc1_shift[-1]] + list(pc1_shift)[:len(pc1_shift)-1])   # TODO: find correct shift

    pc2 = lambdas[0] + pc1
    pc3 = lambdas[1] + pc2
    pc4 = lambdas[3] + pc1

    I = D  # len of kb*Z
    Pc = vector(GF(2), list(pc1)+list(pc2)+list(pc3)+list(pc4))    # len of 4*Z
    Pa = calc_pa(H=H, Pc=Pc, D=D, Zc=Zc, K=K)     # len of (mb-4)*Z

    # x = [i pc pa]
    X = vector(GF(2), list(I)+list(Pc)+list(Pa))
    # print(f"H*X ==0 := {H*X==0}")
    # print(f"mb = {BG.nrows()}: nb = {BG.ncols()} : kb = {BG.ncols()-BG.nrows()}")S
    return X, H, BG


"""
    mb = 42
    nb = 52
    kb = 10
    Zc = 8
    Am = H.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(K)))
    Bm = H.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(K, K + 4*Zc)))
    Cm = H.matrix_from_rows_and_columns(list(range((4*Zc, H.nrows())), list(range(kb * Zc)))
    Dm = H.matrix_from_rows_and_columns(list(range((mb - 4) * Zc)), list(range(kb * Zc, (kb * Zc) + 4 * Zc)))
    
    BG2A, BG2B = BG.matrix_from_rows_and_columns(list(range(Zc)), list(range(10))), BG.matrix_from_rows_and_columns(list(range(4)), list(range(10, 10+4)))
    BG1A , BG1B = BG.matrix_from_rows_and_columns(list(range(Zc)), list(range(22))), BG.matrix_from_rows_and_columns(list(range(4)), list(range(22, 22+4)))
    
    print(f"Am*D:= {Am*D+Bm*Pc} \n Bm*Pc := {Bm*Pc} \n Am*D+Bm*Pc := {Am*D+Bm*Pc}")
"""


