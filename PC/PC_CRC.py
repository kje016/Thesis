# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
#crc24a = x**24 + x**23 + x**18 + x**17 + x**14 + x**11 + x**10 + x**7 + x**6 + x**5 + x**4 + x**3 + x + x**0
#crc24b = x**24 + x**23 + x**6 + x**5 + x + x**0
crc24 = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0


def get_pol(A_seq, channel):
    if channel == "PUCCH" or channel == "PUSCH":
        if len(A_seq) < 12:
            return None
        if 12 <= len(A_seq) <= 12:
            return crc6
        elif 20 <= len(A_seq) <= 1706:
            return crc11
    else:
        return crc24
    # TODO: check if in specification there are more 24-bit crc's for PC
    #elif channel == "PBCH":
    #    return crc24a
    #elif channel == "PBCH":
    #    return crc24a
    #elif channel == "PDCCH":
    #    return crc24a


def bit_long_division(A_seq, pol):
    A = len(A_seq)
    remainder = A_seq[:]
    pos = remainder.index(1)
    divisor = pol.list()[::-1]
    while 1 in remainder[:A-pol.degree()]:
        calc = remainder[pos:pos+(pol.degree()+1)]
        calc = [a+b for a, b in zip(calc, divisor)]

        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        pos = remainder.index(1)
    return remainder[A-pol.degree():]


# bit_long_division returns the remainder, so that in CRC_calc() the message.extend
# acts as the hardware interleave r
def CRC_calc(message, polynomial):
    if polynomial is None:
        return message
    pad = [0]*polynomial.degree()
    message.extend(bit_long_division(message+pad, polynomial))
    return message


def main_CRC(A_seq, channel):
    A_seq = [int(x) for x in A_seq]
    polynomial = get_pol(A_seq, channel)
    output = CRC_calc(A_seq, polynomial)
    return output, polynomial
