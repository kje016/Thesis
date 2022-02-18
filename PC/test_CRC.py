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
    CRC = vector(GF(2), a + res)
    return CRC, C


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
    testing for G matrix
"""
def as_pol(list):
    return sum([c*x**(e) for e,c in enumerate(list[::-1])])


def G_crc(a, pol):
    cc = identity_matrix(GF(2), pol.degree()+len(a))
    I, C = identity_matrix(pol.degree()), Matrix(pol.degree(), pol.degree())
    C[-1] = pol.list()[::-1][1:]

    j = pol.degree()-2
    while j >= 0:
        pol_list = (x * as_pol(C[j+1])).mod(pol).list()[::-1]
        C[j] = [0]*(pol.degree()-len(pol_list)) + pol_list
        j = j-1
    genmat = block_matrix([[I, C]])
    HI = identity_matrix(GF(2), genmat.ncols()-genmat.nrows())
    H = block_matrix([[C.transpose(), HI]])
    breakpoint()
    return vector(GF(2), Matrix(GF(2), a)* genmat), genmat
