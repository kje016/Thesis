# cd Desktop/Thesis/PySageMath/LDPC
from sage.all import *
import LDPC_HelperFunctions as HF


def calc_lambdas(kb, H, Z, D, K):
    result, A = [], H.matrix_from_rows_and_columns(list(range(4 * Z)), list(range(kb*Z)))
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


def Encoding(H, Bi, Zc, D, K, kb, BG):
    # Pc (core parity): can be calculated from submatrices A & B
    lambdas = calc_lambdas(kb, H, Zc, D, K)
    breakpoint()
    PX = HF.P(Bi[1], Zc)
    if BG == 1:
        if Bi == (0,1):
            pc1 = sum(lambdas)
            pc1_shift = vector(GF(2), list(pc1)[1:] + [pc1[0]])
            pc2 = lambdas[0] + pc1_shift
            pc4 = lambdas[3] + pc1_shift
            pc3 = lambdas[2] + pc4
        else:
            pc1 = sum(lambdas)
            pc1_shift = vector(GF(2), list(pc1[Bi[1]:])+ list(pc1[:Bi[1]]))
            pc2 = lambdas[0] + pc1_shift
            pc4 = lambdas[3] + pc1_shift
            pc3 = lambdas[2] + pc4
        """
        if Bi == 0:
            
        """
    else:
        if Bi == 1:
            pc1 = sum(lambdas)
            pc1_shift = vector(GF(2), list(pc1)[1:] + [pc1[0]])
            pc2 = lambdas[0] + pc1_shift
            pc3 = lambdas[1] + pc2
            pc4 = lambdas[3] + pc1_shift
        else:
            pc1_shift = sum(lambdas)
            pc1 = vector(GF(2), [pc1_shift[-1]] + list(pc1_shift)[:len(pc1_shift) - 1])
            pc2 = lambdas[0] + pc1
            pc3 = lambdas[1] + pc2
            pc4 = lambdas[3] + pc1


    Pc = vector(GF(2), list(pc1)+list(pc2)+list(pc3)+list(pc4))    # len of 4*Z
    Pa = calc_pa(H=H, Pc=Pc, D=D, Zc=Zc, K=K)     # len of (mb-4)*Z

    # x = [i pc pa]
    X = vector(GF(2), list(D[:])+list(Pc)+list(Pa))
    # print(f"H*X ==0 := {H*X==0}")
    # print(f"mb = {BG.nrows()}: nb = {BG.ncols()} : kb = {BG.ncols()-BG.nrows()}")S
    if (H*X).hamming_weight() != 0:
        breakpoint()
    return X
