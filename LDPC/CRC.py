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


def bit_long_division(a, pol):
    A = len(a)
    remainder, divisor = a[:], pol.list()[::-1]
    while 1 in remainder[:A-pol.degree()]:
        pos = remainder.index(1)
        calc = [b+c for b, c in zip(remainder[pos:pos+(pol.degree()+1)], divisor)]

        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
    return remainder[A-pol.degree():]


# bit_long_division returns the remainder, so that in CRC_calc() the message.extend
# acts as the hardware interleave r
def CRC_calc(message, polynomial):
    if polynomial is None:
        return message
    pad = [0]*polynomial.degree()
    message.extend(bit_long_division(message+pad, polynomial))
    return message


def main_CRC(a, polynomial):
    a = list(a)
    output = CRC_calc(a, polynomial)
    return output



