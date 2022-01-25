# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24 = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0


def get_pol(A_seq, channel):
    #if channel == "PUCCH" or channel == "PUSCH":
    if len(A_seq) < 12:
        return None
    if 12 <= len(A_seq) <= 19:
        return crc6
    elif 20 <= len(A_seq) <= 1706:
        return crc11
    else:
        return crc24


def bit_long_division(a, pol):
    pad = [0] * pol.degree()
    atemp = a + pad
    A = len(atemp)
    remainder, divisor = atemp[:], pol.list()[::-1]
    pos = remainder.index(1)
    while pos < A-pol.degree():
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        pos = remainder.index(1)
    # print(f"crc_calc := {remainder[A-pol.degree():]}")
    return remainder[A-pol.degree():]


def crc_shift_register(a, pol):
    register, divisor = a[:pol.degree()], pol.list()[::-1]
    pos = len(register)
    while pos < len(a)-len(divisor):
        output = register[0]
        pos += 1
        register = a[pos:pos+len(register)]
        if output == 1:
            register = [a+b for a, b in zip(register, divisor)]
    return register

# bit_long_division returns the remainder, so that in CRC_calc() the message.extend
# acts as the hardware interleave r
def CRC_calc(message, polynomial):
    if polynomial is None:
        return message
    message.extend(bit_long_division(message, polynomial))
    return message


def main_CRC(a, channel):
    polynomial = get_pol(a, channel)
    polynomial = crc24      # TODO: not only use crc24?
    output = CRC_calc(a, polynomial)
    return vector(GF(2), output), polynomial


