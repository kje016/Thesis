from sage.all import *


def RM_main(u, Zc, H, K, K_ap, R, B):
    Kb, colpunct, punct = H.ncols() - H.nrows(), 2 * Zc, floor((K - K_ap) // Zc) * Zc
    A_ap = K - colpunct - punct   # A is the amount of crc bits after removing 2 cols and padding
    E = ceil((B/R)/Zc) * Zc
    pbits = u[K: K + (E-A_ap)]   # getting the parity bits and
    e = list(u[colpunct: B+(K-K_ap-punct)]) + list(pbits)

    te = vector(GF(2), list(u[:colpunct]) + list(e[:A_ap]) + [0]*(punct) + list(pbits))
    Hm = H.matrix_from_rows_and_columns(list(range(E-A_ap)), list(range(K + (E-A_ap))))
    HH = H.matrix_from_rows_and_columns(list(range(E+colpunct+punct - K)), list(range(K + (E-A_ap))))
    breakpoint()
    return vector(ZZ, e), Hm


# punctured bits are unknown
# shortened filler bits are known to be 0
def fill_w_llr(r, Zc, K, K_ap, p, Kb, channel):
    # TODO: hard-coded to be 2??
    punct, short = 2 * Zc, floor((K - K_ap) // Zc) * Zc
    A = Kb - short - punct
    llr = 2/p**2 if channel== 'AWGN' else log((1-p)/p)
    punct_inf = [0] * (2 * Zc)
    if channel == 'AWGN':
        inf_bits = r[:A]*llr
        short_bits = [-oo]*floor((K-K_ap)//Zc) * Zc
        p_bits = r[A:]*llr
    elif channel == 'BSC':
        inf_bits = r[:A]*llr
        short_bits = [-llr]*floor((K-K_ap)//Zc) * Zc
        p_bits = r[A:]*llr
    else:   # channel = BEC
        inf_bits = map(lambda y: 0 if y == 2 else y*oo, r[:A])
        short_bits = [-oo] * (floor((K - K_ap) // Zc) * Zc)
        p_bits = map(lambda y: 0 if y == 2 else y*oo, r[A:])
    te = vector(RealField(10), punct_inf + list(inf_bits) + list(short_bits) + list(p_bits))
    breakpoint()
    return vector(RealField(10), punct_inf + list(inf_bits) + list(short_bits) + list(p_bits))
