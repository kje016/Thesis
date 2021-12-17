# cd Desktop/Thesis/PySageMath/LDPC
# personal access token: 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
from sage.all import *
import LDPC_HelperFunctions as HF
import CRC
# { a0 : [j0,...,j7], a1 : [j0, ..., j7],...}
lifting_size_set = {2: [0, 1, 2, 3, 4, 5, 6, 7], 3: [0, 1, 2, 3, 4, 5, 6, 7], 5: [0, 1, 2, 3, 4, 5, 6],
                    7: [0, 1, 2, 3, 4, 5], 9: [0, 1, 2, 3, 4, 5], 11: [0, 1, 2, 3, 4, 5],
                    13: [0, 1, 2, 3, 4], 15: [0, 1, 2, 3, 4]}

for key, value in lifting_size_set.items():
    crnt_lift = [key * 2 ** j for j in value]
    lifting_size_set.update({key : crnt_lift})


def E(perm, z):
    if perm == -1:
        return zero_matrix(GF(2), z, z)
    if perm == 0:
        return identity_matrix(GF(2), z)
    else:
        res_matrix = []
        for i in range(z):
            row = [0]*z
            row[(perm+i)%z]= 1
            res_matrix.append(row)
        return matrix(GF(2), res_matrix)


def Protograph(base_matrix, z):
    protograph = Matrix(GF(2), len(base_matrix.rows())*z, len(base_matrix.columns())*z)
    for i_row, row in enumerate(base_matrix):
        for i_col, col in enumerate(row):
            protograph.set_block(i_row*z, i_col*z, E(row[i_col], z))
    return protograph


B = 20
b_bits = [1] * B
bg = 2
L, C, B_ap = HF.get_param(bg=bg, B=B)
K_ap, kb = HF.determine_kb(B=B, B_ap=B_ap, C=C, bg=bg)

# TODO: calculate crk
kb = HF.determine_kb_BG2(B)
Z, iLS, K = HF.determine_Z(bg, kb, lifting_size_set, K_ap)
#crk = HF.calc_crk(C=C, K_ap=K_ap, K=K, L=L, b_bits=b_bits)
BG = HF.get_base_matrix(bg, iLS, Z)
H = Protograph(BG, Z)
A = H.matrix_from_rows_and_columns(list(range(4 * Z)), list(range(kb * Z)))
B = H.matrix_from_rows_and_columns(list(range(4*Z)), list(range((kb*Z), (kb*Z)+4*Z)))
# C = H.matrix_from_rows_and_columns(list(range((mb-4)*Z)), list(range(kb*Z)))
# D = H.matrix_from_rows_and_columns(list(range((mb-4)*Z)), list(range(kb*Z, (kb*Z)+4*Z)))
b_bits.extend([0] * ((kb * Z) - len(b_bits)))     #padding inf_bits
b_bits = vector(GF(2), b_bits)
breakpoint()
def lambda_it(kb, H, Z, b_bits):
    result, A = [], H.matrix_from_rows_and_columns(list(range(4 * Z)), list(range(kb * Z)))
    for j in range(4):
        lam_ji = vector(GF(2), Z)
        for l in range(kb):
            a = A.matrix_from_rows_and_columns(list(range((j*4), (j+1)*4)), list(range(l*Z, (l+1)*Z)))
            lam_ji += (a*b_bits[l*Z:(l*Z)+Z])
        result.append(lam_ji)
    return result


# Initialisation of x = [i pc pa]
mb = BG.ncols()-kb
nb = kb + mb
len_pc, len_pa = 4*Z, (mb-4)*Z

lambdas = lambda_it(kb, H, Z, b_bits)
# Pc (core parity): can be calculated from submatrices A & B
pc1 = sum(lambdas)
if B[2, B.ncols()-B.rows()] == -1:

pc1_shift = vector(GF(2), [pc1[-1]] + list(pc1)[:len(pc1)-1])
pc2 = lambdas[0] + pc1_shift
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
pc3 = lambdas[1] + pc2
pc4 = lambdas[2] + pc1_shift
# Pa (additional parity): calculated from information core and parit bits using C & D
def calc_pa(H, Pc, b_bits):
    C = H.matrix_from_rows_and_columns(list(range((mb-4)*Z)), list(range(kb*Z)))
    D = H.matrix_from_rows_and_columns(list(range((mb-4)*Z)), list(range(kb*Z, (kb*Z)+4*Z)))
    return C*b_bits + D*Pc


I = b_bits  # len of kb*Z
Pc = vector(GF(2), flatten(list(pc1)+list(pc2)+list(pc3)+list(pc4)))    # len of 4*Z
Pa = calc_pa(H=H, Pc=Pc, b_bits=b_bits)     # len of (mb-4)*Z
# x = [i pc pa]
X = vector(GF(2), list(I)+list(Pc)+list(Pa))
print(X)
print()
print(A*b_bits + B*Pc)
breakpoint()
