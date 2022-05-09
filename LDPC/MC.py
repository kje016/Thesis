from sage.all import *

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()


def CRC(a, A, pol):
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