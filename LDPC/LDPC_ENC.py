# cd Desktop/Thesis/PySageMath/LDPC
from sage.all import *
# { a0 : [j0,...,j7], a1 : [j0, ..., j7],...}
lifting_size_set = {2: [0, 1, 2, 3, 4, 5, 6, 7], 3: [0, 1, 2, 3, 4, 5, 6, 7], 5: [0, 1, 2, 3, 4, 5, 6],
                    7: [0, 1, 2, 3, 4, 5], 9: [0, 1, 2, 3, 4, 5], 11: [0, 1, 2, 3, 4, 5],
                    13: [0, 1, 2, 3, 4], 15: [0, 1, 2, 3, 4]}

for key, value in lifting_size_set.items():
    crnt_lift = [key * 2 ** j for j in value]
    lifting_size_set.update({key : crnt_lift})

information_bits = 20  # evt len(codeword)
bg = 2


def get_base_matrix(bg, ils, zc):
    matrix = []
    f = open("base_matrices\\" + f"NR_{bg}_{ils}_{zc}.txt", "r")
    txt_m = f.read().split('\n')
    for row in txt_m:
        row = list(map(int, [c for c in row.split(' ') if c!= '']))
        matrix.append(row)
    matrix.pop(-1)          # matrix contains an empty list at the end
    return Matrix(matrix)


def determine_kb_BG2(k_bits):
    if k_bits > 640:
        return 10
    if (k_bits > 560 and k_bits <= 640):
        return 9
    if (k_bits > 192 and k_bits <= 560):
        return 8
    else:
        return 6


# Z = a * 2**j
def determine_Z(input_kb):
    lifted_set = {}
    #breakpoint()
    for key,value in lifting_size_set.items():
        crnt_lift = [x for x in value if input_kb*x >= information_bits]
        crnt_lift = min(crnt_lift)
        lifted_set.update({key : crnt_lift})

    lifting_size = min(lifted_set.values())
    set_index = list(lifted_set.values())
    set_index = set_index.index(lifting_size)
    return lifting_size,set_index


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
            sub_matrix = E(row[i_col], z)
            protograph.set_block(i_row*z, i_col*z, sub_matrix)

    return protograph


kb = determine_kb_BG2(information_bits)
Z, iLS = determine_Z(kb)

BG = get_base_matrix(bg, iLS, Z)
H = Protograph(BG, Z)


def calc_lamb():
    get_rows = [i for i in range(4 * Z)]
    get_cols = [i for i in range(kb * Z)]
    A = H.matrix_from_rows_and_columns(get_rows, get_cols)
    res_lambda = []
    for j in range(0, A.nrows(), Z):
        temp_comp = vector(GF(2),4)
        for i in range(0, len(inf_bits), Z):
            a = A.matrix_from_rows_and_columns([x for x in range(j, j + Z)], [x for x in range(i, i + Z)])
            inf_trans = inf_bits[i:i + Z]
            lam = a*inf_trans
            temp_comp = temp_comp+lam
        res_lambda.append(temp_comp)
    return res_lambda


def vectors_to_vector(input_list):
    output_vec = []
    for vec in input_list:
        output_vec.extend(vec.list())
    return vector(GF(2), output_vec)

# Initialisation of x = [i pc pa]
inf_bits = [1] * information_bits
inf_bits.extend([0] * ((kb*Z)-len(inf_bits)) )     #padding inf_bits
inf_bits = vector(GF(2), inf_bits)

mb = H.ncols()//Z-kb
nb = kb + mb
len_pc, len_pa = 4*Z, (mb-4)*Z


lambdas = calc_lamb()
pc1 = vector(GF(2),4)
for x in lambdas:
    pc1 = pc1 + x

pc1 = list(pc1)
pc1.append(pc1.pop(0))
pc1 = vector(GF(2), pc1)
pc2 = lambdas[0] + pc1
pc3 = lambdas[1] + pc2
pc4 = lambdas[3] + pc1


def calc_pa():
    nrows = [i for i in range(4 * Z, H.nrows())]
    ncols = [i for i in range((kb * Z) + 4*Z)]
    CD = H.matrix_from_rows_and_columns(nrows, ncols)
    vec = vector(GF(2), list(inf_bits)+list(Pc))
    return CD*vec


# x = [i pc pa]
I = inf_bits
Pc = vectors_to_vector([pc1,pc2,pc3,pc4])
Pa = calc_pa()

X = vector(GF(2), list(I)+list(Pc)+list(Pa))

print(X)
breakpoint()

