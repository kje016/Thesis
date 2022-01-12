# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions
######################################################################################
var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

crc6 = x**6 + x**5 + x**0
crc11 = x**11 + x**10 + x**9 + x**5 + x**0
crc16 = x**16 + x**12 + x**5 + x**0
crc24 = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0

c1 = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
c2 = [1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# here: A is the index+1 of the last information bit
def div(a, A, pol):
    #A = len(a)
    remainder, divisor = a[:], pol.list()[::-1]
    breakpoint()
    while 1 in remainder[:A+1]:
        pos = remainder.index(1)
        print(pos, remainder[:A+1])
        calc = [b+c for b, c in zip(remainder[pos:], divisor)]
        remainder = remainder[:pos] + calc + remainder[pos+len(calc):]

    return remainder
######################################################################################


K_max_il = 164
pattern = "0 2 4 7 9 14 19 20 24 25 26 28 31 34 42 45 49 50 51 53 54 56 58 59 61 62 65 66 " \
          "67 69 70 71 72 76 77 81 82 83 87 88 89 91 93 95 98 101 104 106 108 110 111 113 115 118 119 120 " \
          "122 123 126 127 129 132 134 138 139 140 1 3 5 8 10 15 21 27 29 32 35 43 46 52 55 57 60 63 " \
          "68 73 78 84 90 92 94 96 99 102 105 107 109 112 114 116 121 124 128 130 133 135 141 6 11 16 22 30 " \
          "33 36 44 47 64 74 79 85 97 100 103 117 125 131 136 142 12 17 23 37 48 75 80 86 137 143 13 18 " \
          "38 144 39 145 40 146 41 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163"
interleaver_pattern = [int(x) for x in list(pattern.split(' '))]
#int_test = HelperFunctions.standard_to_list('Interleaving_pattern.txt')


def get_pi(K):
    PI, h = [], K_max_il-K
    for m in range(K_max_il):
        if interleaver_pattern[m] >= h:
            PI.append(interleaver_pattern[m]-h)
    return PI


def interleaver(flag_param, c_seq):
    if not flag_param:
        return c_seq
    PI = get_pi(len(c_seq))
    output = [c_seq[value] for value in PI]
    print(f"PI:={PI}")
    tess = div(output, max(PI.index(0), PI.index(1), PI.index(2)), crc24)
    breakpoint()
    return output


def main_bit_interleaver(channel, c):
    flag = True     # channel in ["PBCCH", "PDCCH"]    # interleaver is used in PBCH, PDCCH, bypassed for PUCCH & PUSCH
    c_ap = interleaver(flag, c)
    return flag, c_ap
