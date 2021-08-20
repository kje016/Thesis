from sage.all import *
import sys

def main():
    # the realiability sequence is 1-indexed
    with open('Reliability_Sequence.txt') as f:
        reliability_sequence = f.read()
    reliability_sequence = reliability_sequence.replace(' ', '')
    reliability_sequence = reliability_sequence.split(",")
    reliability_sequence = [int(x) for x in reliability_sequence]

    K = int(len(sys.argv[1]))
    n = (log(K,2) +1)
    print("K: ", K)
    c = list([1]); c = [int(x) for x in c]
    N = 2 ** n
    # Basic polarization kernel: G2
    field = GF(2);
    type(field)
    MS = MatrixSpace(field, 2)
    G = MS.matrix([1, 0, 1, 1])

    def Kronecker_Product(base_matrix, current_matrix):
        dim = current_matrix.nrows() * 2
        G1N = MatrixSpace(field, dim)
        G2N = MatrixSpace(field, dim)
        G3N = MatrixSpace(field, dim)
        G4N = MatrixSpace(field, dim)

        G1N = matrix([base_matrix[0, 0] * g for g in current_matrix])
        G2N = matrix([base_matrix[0, 1] * g for g in current_matrix])
        G3N = matrix([base_matrix[1, 0] * g for g in current_matrix])
        G4N = matrix([base_matrix[1, 1] * g for g in current_matrix])

        result_matrix = block_matrix([[G1N, G2N], [G3N, G4N]])
        return result_matrix

    # Channel transformation Matrix: GN
    GN = Kronecker_Product(G, G)
    for i in range(1, n - 1):
        GN = Kronecker_Product(G, GN)

    # Frozen Set
    # 'elem-1' since reliability_sequence is 1-indexed
    F = []
    for index, elem in enumerate(reliability_sequence):
        if elem <= N:
            F.append(elem - 1)
        if len(F) >= K:
            break

    # Ui = 0 if i elem F
    u = [1] * N
    c_geti = 0

    for index, elem in enumerate(u):
        if index in F:
            u[index] = 0
        else:
            u[index] = c[c_geti]
            c_geti += 1

    uv = vector(field, u)
    dv = uv*GN
    d = list(dv)
    print("Frozen set:", F)

    return d

if __name__ == "__main__":
    d = main()
    print("codeword:")
    print(d)