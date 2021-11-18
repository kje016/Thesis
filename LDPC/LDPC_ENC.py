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

n_inf_bits = 20
inf_bits = [1] * n_inf_bits
bg = 2


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


kb = HF.determine_kb_BG2(n_inf_bits)
Z, iLS = HF.determine_Z(kb, lifting_size_set, n_inf_bits)
BG = HF.get_base_matrix(bg, iLS, Z)
H = Protograph(BG, Z)

inf_bits.extend([0] * ((kb*Z)-len(inf_bits)) )     #padding inf_bits
inf_bits = vector(GF(2), inf_bits)


def calc_lambdas(kb, Z, H):
    A = H.matrix_from_rows_and_columns(list(range(4 * Z)), list(range(kb * Z)))
    B = H.matrix_from_rows_and_columns(list(range(4*Z)), list(range(4*Z, (4*Z)*2)))
    res_lambda = []
    for j in range(0, A.nrows(), Z):
        temp_comp = vector(GF(2), 4)
        for i in range(0, len(inf_bits), Z):
            a = A.matrix_from_rows_and_columns(list(range(j, j + Z)), list(range(i, i + Z)))
            #b = B.matrix_from_rows_and_columns(list(range(j, j+Z)), list(range(j,j+Z)))
            temp_comp = temp_comp+(a*inf_bits[i:i + Z])
        res_lambda.append(temp_comp)
    return res_lambda


# Initialisation of x = [i pc pa]
mb = H.ncols()//Z-kb
nb = kb + mb
len_pc, len_pa = 4*Z, (mb-4)*Z

lambdas = calc_lambdas(kb, Z, H)
breakpoint()
# Pc (core parity): can be calculated from submatrices A & B
pc1 = vector(GF(2),4)
for x in lambdas:
    pc1 = pc1 + x

breakpoint()
pc1 = list(pc1)
pc1.append(pc1.pop(0))
pc1 = vector(GF(2), pc1)
pc2 = lambdas[0] + pc1
pc3 = lambdas[1] + pc2
pc4 = lambdas[3] + pc1


# Pa (additional parity): calculated from information core and parit bits using C & D
def calc_pa():
    nrows, ncols = [i for i in range(4 * Z, H.nrows())], [i for i in range((kb * Z) + 4*Z)]
    CD = H.matrix_from_rows_and_columns(nrows, ncols)
    vec = vector(GF(2), list(inf_bits)+list(Pc))
    return CD*vec


# x = [i pc pa]
I = inf_bits
Pc = vector(GF(2), flatten(list(pc1)+list(pc2)+list(pc3)+list(pc4)))
Pa = calc_pa()

X = vector(GF(2), list(I)+list(Pc)+list(Pa))

print(X)
breakpoint()

