from sage.all import *

import PC_CRC

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
    res = vector(Matrix(GF(2), a)*C).list()
    CRC = vector(GF(2), list(a) + res)
    return CRC


def Chun_mat(a, A, pol):
    g = pol
    n = g.degree()
    m = sum([a*x**i for i,a in enumerate(a)])
    #m = (x**2 + x + x**0);
    if m.degree() < 11:
        m = m * x **(11 - m.degree());
    k = m.degree() + 1;

    # More efficient version (from the paper)
    D = zero_matrix(k, n);
    D[k - 1] = [g[n - i] for i in list(range(1, n+1))];

    t = k - 1;
    while t > 0:
        for j in range(1, n):
            D[t - 1, j - 1] = Mod(D[t, j] + D[t, 0] * g[n - j], 2);
            D[t - 1, n - 1] = D[t, 0] * g[0];
        t = t - 1;
    mm = [m[k - 1 - i] for i in range(k)]   # TODO: this reverses the list?

    return Matrix(GF(2), mm)*D, D


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
    res = vector(Matrix(GF(2), a)*C).list()
    return res


def ICRC_check(a, A, pol):
    P = pol.degree()
    C = zero_matrix(GF(2), A, P)
    C[-1] = vector(GF(2), pol.list()[::-1][1:])
    k = A-2
    while k >= 0:
        for i in range(P-1):
            C[k, i] = C[k+1, i+1] + C[k+1, 0] * pol[P - (i+1)]
            C[k, P-1] = C[k+1, 0] * pol[0]
        k = k-1
    res = vector(Matrix(GF(2), a)*C).list()
    return res


"""
    From "Interleaved CRC for Polar Codes"
"""

