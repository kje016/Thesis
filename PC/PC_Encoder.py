# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


def Kronecker_Product(base_matrix, current_matrix):
    G1N = matrix([base_matrix[0, 0] * g for g in current_matrix])
    G2N = matrix([base_matrix[0, 1] * g for g in current_matrix])
    G3N = matrix([base_matrix[1, 0] * g for g in current_matrix])
    G4N = matrix([base_matrix[1, 1] * g for g in current_matrix])

    result_matrix = block_matrix([[G1N, G2N], [G3N, G4N]])
    return result_matrix


def main_encoder(u, N, n_pc, n_wm_pc):
    # Basic polarization kernel: G2
    MS = MatrixSpace(GF(2), 2)
    G2 = MS.matrix([1, 0, 1, 1])

    # Channel transformation Matrix: GN
    GN = Kronecker_Product(G2, G2)
    for i in range(1, log(N, 2)-1):
        GN = Kronecker_Product(G2, GN)
    breakpoint()
    dv = vector(GF(2), u) *GN
    return list(dv)

c = main_encoder([0, 0, 0, 1, 0, 0, 1, 1], 8, 0, 0)