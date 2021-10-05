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

CRC_s = {crc6 : 6, crc11: 11, crc16: 16, crc24a:24, crc24b: 24, crc24c: 24}
#PC_crc = {(12, 19): crc6}

test_pol = x**3+x+x**0
m = [1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]
B = len(m)


""" getting maximum code block size Kb"""
if bg ==1:
    Kcb = 8448
else:
    Kcb = 3840


""" Total number of code blocks 'C' is determined by:    """
if B <= kcb:
    L = 0
    C = 1
    B_ap = b
else:
    L = 24
    C = ceil(b/(kcb-L))
    B_ap = b + C*L


""" number of bits 'K' in each code block: """
K_ap = B_ap/C
if bg == 1:
    Kb = 22
else:       # bg2
    if B > 640:
        Kb = 10
    elif B > 560:
        Kb = 9
    elif B > 192:
        Kb = 8
    else:
        Kb = 6


"""
Find the minimum value of 'Z' denoted as Zc s.t Kb * Zc >= K_ap
set K = 22*Zc for BG1
set K = 10*Zc for BG2
"""


def pad_remainder(k_seq, input_pol):
    return vector(GF(2), k_seq[:] + [0]*input_pol.degree())

def bit_long_division(message, pol):
    remainder = list(message[:])
    pos = remainder.index(1)
    divisor = pol.list()[::-1]
    while 1 in remainder[0:len(m)-1] or not pos > len(m)-1:
        #breakpoint()
        calc = remainder[pos:pos+(pol.degree()+1)]
        calc = [a+b for a, b in zip(calc, divisor)]

        remainder = remainder[0:pos] + calc + remainder[pos+pol.degree()+1:]
        pos = remainder.index(1)
        print(f"{remainder}: {pos}")
        print()
    breakpoint()
    return remainder[-pol.degree():]



cword = pad_remainder(m, test_pol)
check = bit_long_division(cword, test_pol)
m.extend(check)
breakpoint()

L, C, B_ap = get_tot_c(Bm )

K_ap, Kb = bits_per_block()


