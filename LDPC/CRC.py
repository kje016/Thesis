# cd Desktop/Thesis/PySageMath/LDPC
from sage.all import *

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24a = x**24 + x**23 + x**18 + x**17 + x**14 + x**11 + x**10 + x**7 + x**6 + x**5 + x**4 + x**3 + x + x**0
crc24b = x**24 + x**23 + x**6 + x**5 + x + x**0
crc24c = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0


def get_pol(A):
    if A <= 3824:
        return crc16
    else:
        return crc24a


def CRC(a, A, pol):
    if len(a) < 12:
        return a, None
    P = pol.degree()
    C = zero_matrix(GF(2), A, P)
    C[-1] = vector(GF(2), pol.list()[::-1][1:])
    k = A-2
    while k >= 0:
        for i in range(P-1):
            C[k, i] = C[k+1, i+1] + C[k+1, 0] * pol[P - (i+1)]
        C[k, P-1] = C[k+1, 0] * pol[0]
        k = k-1
    res = (vector(GF(2), a)*C).list()
    CRC = vector(GF(2), list(a) + res)
    return CRC, C

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
