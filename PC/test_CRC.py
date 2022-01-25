from sage.all import *

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()
g = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0
P = g.degree()


def CRC(a, A):
    C = zero_matrix(GF(2),A, P)
    C[-1] = vector(GF(2), g.list()[::-1][1:])
    t = A-1
    while t >= 0:
        for j in range(1, P):
            C[t-1, j-1] = C[t, j] + C[t, 0] * g[P - j]
            C[t - 1, P - 1] = C[t, 0] * g[0]
        t = t-1

    res = vector(Matrix(GF(2), a)*C).list()
    CRC = vector(GF(2), list(a) + res)
    return CRC, g


def CRC_check(a, A):
    C = zero_matrix(GF(2),A, P)
    C[-1] = vector(GF(2), g.list()[::-1][1:])
    t = A-1
    while t >= 0:
        for j in range(1, P):
            C[t-1, j-1] = C[t, j] + C[t, 0] * g[P - j]
            C[t - 1, P - 1] = C[t, 0] * g[0]
        t = t-1

    res = vector(Matrix(GF(2), a)*C).list()
    return res
