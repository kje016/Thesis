from sage.all import *

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()


def CRC(a, A, pol):
    P = pol.degree()+1
    C = zero_matrix(GF(2), A, P)
    C[-1] = vector(GF(2), pol.list()[::-1])
    t = A-1
    while t >= 0:
        for j in range(1, P):
            C[t-1, j-1] = C[t, j] + C[t, 0] * pol[P - j]
            C[t - 1, P - 1] = C[t, 0] * pol[0]
        t = t-1

    res = vector(Matrix(GF(2), a)*C).list()
    CRC = vector(GF(2), list(a) + res)
    return CRC


def CRC_check(a, A, pol):
    P = pol.degree()
    C = zero_matrix(GF(2), A, P)
    C[-1] = vector(GF(2), pol.list()[::-1][1:])
    t = A-1
    while t >= 0:
        for j in range(1, P):
            C[t-1, j-1] = C[t, j] + C[t, 0] * pol[P - j]
            C[t - 1, P - 1] = C[t, 0] * pol[0]
        t = t-1

    res = vector(Matrix(GF(2), a)*C).list()
    return res
