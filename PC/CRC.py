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


def ICRC_check(a, A, iPI, K, PI):
    crcbits= []
    for elem in iPI[::-1]:
        crcbits.insert(0, a[elem])

    infbits = [i for i in a]
    for a in iPI[::-1]:
        print(a)
        infbits.pop(a)

    cword = vector(GF(2), list(infbits) + list(crcbits))
    breakpoint()
    #cword = list(map(lambda x: a.append(a.pop(x)), iPI))
    #cword = vector(GF(2), list(a:[]))
    pol = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0
    P = pol.degree()
    C = zero_matrix(GF(2), K, P)
    C[-1] = vector(GF(2), pol.list()[::-1][1:])
    k = K-2
    while k >= 0:
        for i in range(P-1):
            C[k, i] = C[k+1, i+1] + C[k+1, 0] * pol[P - (i+1)]
            C[k, P-1] = C[k+1, 0] * pol[0]
        k = k-1
    cc = Matrix(GF(2), [C[a] for a in PI])
    breakpoint()
    testi = cword*cc
    Cap = cc[:A].rows()
    #Cap.append(vector(GF(2), [a for a in cc[iPI[0]]]))
    Cap = Matrix(GF(2), Cap)

    res = vector(GF(2), a) *Cap
    print(f'res:{res}')
    breakpoint()
    return res

