from sage.all import *


def P(perm, z):
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
            protograph.set_block(i_row * z, i_col * z, P(row[i_col], z))
    return protograph


def calc_crk(C, K_ap, K, L, b_bits):
    s, crk = 0, []
    for r in range(C):
        for k in range(K_ap-L):
            crk.append(b_bits[s])
            s += 1
        if C > 1:
            print("calc_crk() need implementing when C>1")
    crk.extend([None]*(K-K_ap))
    return crk


def get_base_matrix(bg, ils, zc):
    matrix = []
    f = open("base_matrices\\" + f"NR_{bg}_{ils}_{zc}.txt", "r")
    txt_m = f.read().split('\n')
    for row in txt_m:
        row = list(map(int, [c for c in row.split(' ') if c!= '']))
        matrix.append(row)
    matrix.pop(-1)          # matrix contains an empty list at the end
    return Matrix(matrix)


# C := Input to channel coding
# D := Bits after encoding
def get_d_c(Zc, K, C):
    D = []
    for k in range(2*Zc, K):
        if C[k] != None:
            D.append(C[k])
        else:
            C[k] = 0
            D.append(None)
    return D, C



def print_H_square(H, i_row, i_col, zc):
    print(H.matrix_from_rows_and_columns(list(range(i_row*zc, (i_row+1)*zc)), list(range(i_col*zc, (i_col+1)*zc))))
