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

