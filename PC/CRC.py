from sage.all import *

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24 = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0


def get_pol(A, I_IL):
    if I_IL == 1:
        return crc24
    else:
        if A < 12:
            return x**0
        elif 12 <= A <= 19:
            return crc6
        elif 20 <= A <= 1706:
            return crc11
        else:
            return crc24


def gen_CRC_mat(A, pol):
    P = pol.degree()
    C = zero_matrix(GF(2), A, P)
    C[-1] = vector(GF(2), pol.list()[::-1][1:])
    k = A-2
    while k >= 0:
        for i in range(P-1):
            C[k, i] = C[k+1, i+1] + C[k+1, 0] * pol[P - (i+1)]
        C[k, P-1] = C[k+1, 0] * pol[0]
        k = k-1
    return C



def CRC(a, A, pol, I_IL, PI):
    if len(a) < 12 and I_IL == 0:
        return a, None
    P = pol.degree()
    C = gen_CRC_mat(A, pol)
    res = (vector(GF(2), a)*C).list()
    CRC = vector(GF(2), list(a) + res)
    H = block_matrix(GF(2), [[C.transpose(), matrix.identity(C.ncols())]])
    if I_IL:
        H = column_matrix(GF(2), [H.column(a) for a in PI])

    return CRC, H


def CRC_check(a, A, pol):
    P = pol.degree()
    C = zero_matrix(GF(2), A, P)
    C[-1] = vector(GF(2), pol.list()[::-1][1:])
    k = A-2
    while k >= 0:
        for i in range(P-1):
            C[k, i] = C[k+1, i+1] + C[k+1, 0] * pol[P - (i+1)]
            C[k, P-1] = C[k+1, 0] * pol[0]
        k = k-1
    res = vector(Matrix(GF(2), a)*C)
    return res


def test_I_CRC(c, cap, a, PI, pol):
    cap = vector(GF(2), cap)
    A = len(cap)-pol.degree()
    C = gen_CRC_mat(A, pol)
    #G = block_matrix(GF(2), [[matrix.identity(C.nrows()), C]])
    H = block_matrix(GF(2), [[C.transpose(), matrix.identity(C.ncols())]])
    Hperm = column_matrix(GF(2), [H.column(a) for a in PI])
    """
    iPI = [PI.index(a) for a in PI if a >= A]
    h1 = column_matrix(GF(2), [H.column(a) for a in PI[:iPI[0]+1]])
    for i in range(len(iPI)):
        hi = column_matrix(GF(2), [H.column(a) for a in PI[:iPI[i]+1]])
        ci = cap[:iPI[i]+1]
        print(f'h1*c={hi*ci}')
    """
    #Hperm = Matrix(GF(2), [Hperm[a] for a in PI])
    return Hperm
