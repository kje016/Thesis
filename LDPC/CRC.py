# cd Desktop/Thesis/PySageMath/LDPC
import sys
from sage.all import *
"""
From 5G standard.

"""
var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24a = x**24 + x**23 + x**18 + x**17 + x**14 + x**11 + x**10 + x**7 + x**6 + x**5 + x**4 + x**3 + x + x**0
crc24b = x**24 + x**23 + x**6 + x**5 + x + x**0
crc24c = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0

#CRC_s = {crc6: 6, crc11: 11, crc16: 16, crc24a: 24, crc24b: 24, crc24c: 24}


def pad_remainder(k_seq, input_pol):
    return vector(GF(2), k_seq[:] + [0]*input_pol.degree())


def bit_long_division(A_seq, message, pol):
    A = len(A_seq)
    remainder = list(message[:])
    pos = remainder.index(1)
    divisor = pol.list()[::-1]
    while 1 in remainder[0:A-1] or not pos > A-1:
        #breakpoint()
        calc = remainder[pos:pos+(pol.degree()+1)]
        calc = [a+b for a, b in zip(calc, divisor)]

        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        pos = remainder.index(1)
        print(f"{remainder[:A]}: {pos}")
        print()
    return remainder[A:]

def CRC_calc(message, polynomial):
    pad = [0]*polynomial.degree()
    return  message + bit_long_division(message+pad, polynomial)


def paramater_decision(bg, B, Zc):
    Kcb = 0
    L = 0
    C = 0
    K = 0
    B_ap = 0
    """ getting maximum code block size Kb"""
    if bg == 1:
        Kcb = 8448
    else:
        Kcb = 3840


    """ Total number of code blocks 'C' is determined by:    """
    if B <= Kcb:
        L = 0
        C = 1
        B_ap = B
    else:
        L = 24
        C = ceil(B/(Kcb-L))
        B_ap = B + C*L


    """ number of bits 'K' in each code block: """
    K_ap = B_ap//C

    """
    Find the minimum value of 'Z' denoted as Zc s.t Kb * Zc >= K_ap
    set K = 22*Zc for BG1
    set K = 10*Zc for BG2
    """
    if bg == 1:
        K = 22 * Zc
    else:
        K = 10 * Zc

    return Kcb, L, C, K, K_ap, B_ap


def calc_crk(K_ap, input_pol, B_seq, C):
    L = input_pol.degree()
    B = len(B_seq)
    code_block = [B_seq[i:K_ap] for i in range(0, K_ap-L, B)]
    if C > 1:
        for k in code_block:
            code_block.append(calc_crk(CRC_calc(k, crc24b)))
    return code_block


def main(A_seq, bg, Kb, Zc):
    breakpoint()
    B_seq = CRC_calc(A_seq, crc24a)
    B = len(B_seq)
    Kcb, L, C, K, K_ap, B_ap = paramater_decision(bg, B, Zc)
    C_seq = calc_crk(K_ap=K_ap, input_pol=crc24b, B_seq=B_seq, C=C)
    print(f"bit_sequence after CRC := {B}")
    return C_seq



if __name__ == "__main__":
    main(A_seq=sys.argv[1], bg=sys.argv[2], Kb=sys.argv[3], Zc=sys.argv[4])




