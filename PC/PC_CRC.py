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
    pad = [0] * (pol.degree()+1)
    atemp = a + pad
    A = len(atemp)
    remainder, divisor = atemp[:], pol.list()[::-1]
    try:
        pos = remainder.index(1)    # TODO: crashes in 0 codeword
    except:
        [0]*pol.degree()

    while pos < len(a):
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = len(a)
    return remainder[len(a):-1]


def CRC_checksum(c, pol):
    atemp = c[:-pol.degree()+1]
    remainder, divisor = c[:], pol.list()[::-1]
    try:
        pos = remainder.index(1)    # TODO: crashes in 0 codeword
    except:
        [0]*len(remainder)
    while pos < len(c[:-pol.degree()+1]):
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = len(atemp)
    return remainder



def interleave_check(a, pol):
    pad = [0] * (pol.degree()+1)
    atemp = a + pad
    A = len(atemp)
    remainder, divisor = atemp[:], pol.list()[::-1]
    pos = remainder.index(1)    # TODO: crashes in 0 codeword
    while pos < len(a):
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        try:
            pos = remainder.index(1)
        except:
            pos = len(a)
    return remainder[len(a)+1:-1]

# bit_long_division returns the remainder, so that in CRC_calc() the message.extend
# acts as the hardware interleave r
def CRC_calc(message, polynomial):
    if polynomial is None:
        return message
    output = message + bit_long_division(message, polynomial)
    return vector(GF(2), output)


