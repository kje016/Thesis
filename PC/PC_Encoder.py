# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


def Kronecker_Product(base_matrix, current_matrix):
    G1N = matrix([base_matrix[0, 0] * g for g in current_matrix])
    G2N = matrix([base_matrix[0, 1] * g for g in current_matrix])
    G3N = matrix([base_matrix[1, 0] * g for g in current_matrix])
    G4N = matrix([base_matrix[1, 1] * g for g in current_matrix])

    result_matrix = block_matrix([[G1N, G2N], [G3N, G4N]])
    return result_matrix


def gen_G(n):
    MS = MatrixSpace(GF(2), 2)
    G2 = MS.matrix([1, 0, 1, 1])
    # Channel transformation Matrix: GN
    GN = Kronecker_Product(G2, G2)
    for i in range(n - 2):
        GN = Kronecker_Product(G2, GN)
    if n == 1:
        return G2
    return GN
