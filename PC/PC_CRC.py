# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24 = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0


# TODO: dbl check if this is correct
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


def bit_long_division(a, pol):
    if a.count(1) == 0:
        return [0]*pol.degree()

    pad = [0] * (pol.degree())
    remainder, divisor = a + pad, pol.list()[::-1]
    pos = remainder.index(1)
    while pos < len(a):
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = len(a)
    return remainder[len(a):]


def CRC_checksum(c, pol):
    atemp = c[:-pol.degree()]
    remainder, divisor = c[:], pol.list()[::-1]
    try:
        pos = remainder.index(1)
    except:
        return c
    while pos < len(c[:-pol.degree()]):
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = len(atemp)
    return remainder


def Itess(c, pol):
    if c.count(1) == 0:
        breakpoint()
        return [0]*pol.degree()

    remainder, divisor = c[:], pol.list()[::-1]
    pos = remainder.index(1)
    while pos < len(c)-1:
        calc = [b+c for b, c in zip(remainder[pos:], pol.list()[::-1])]
        remainder = remainder[:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = oo
    return remainder


def ICRC_check(c, PI, pol, A):
    if c.count(1) == 0:
        return [0]*pol.degree()
    infs, checks = [], []
    for i in range(len(c)):
        if PI[i] >= A:
            checks.append(c[i])
        else:
            infs.append(c[i])
    remainder, divisor = infs + [0]*(A-len(infs)) + checks + [0]*(pol.degree()-len(checks)), pol.list()[::-1]
    pos = remainder.index(1)
    while pos < len(infs):
        calc = [b+c for b, c in zip(remainder[pos:], pol.list()[::-1])]
        remainder = remainder[:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = oo
    return remainder


# bit_long_division returns the remainder, so that in CRC_calc() the message.extend
# acts as the hardware interleave r
def CRC_calc(message, polynomial):
    if polynomial is None:
        return message
    output = message + bit_long_division(message, polynomial)
    tess = CRC_checksum(output, polynomial)
    return vector(GF(2), output)
