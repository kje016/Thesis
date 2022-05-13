from sage.all import *


def RM_main(u, Zc, H, K, K_ap, rate, B):
    colpunct, punct = 2 * Zc, floor((K - K_ap) // Zc) * Zc
    A_ap = K - colpunct - punct   # A is the amount of crc bits after removing 2 cols and padding
    E = ceil((B / rate) / Zc) * Zc
    pbits = u[K: K + (E-A_ap)]   # getting the parity bits and
    e = list(u[colpunct: B+(K-K_ap-punct)]) + list(pbits)
    Hm = H.matrix_from_rows_and_columns(list(range(E-A_ap)), list(range(K + (E-A_ap))))
    #te = vector(GF(2), list(u[:colpunct]) + list(e[:A_ap]) + [0]*(punct) + list(pbits))
    #print(f'len reconstructed e= {len(te)}')
    #HH = H.matrix_from_rows_and_columns(list(range(E+colpunct+punct - K)), list(range(K + (E-A_ap))))
    #if (Hm*te).hamming_weight() != 0:
        #breakpoint()
    return vector(RealField(10), e), Hm


# punctured bits are unknown
# shortened filler bits are known to be 0
def fill_w_llr(r, Zc, K, K_ap, p, channel):
    # TODO: hard-coded to be 2??
    colpunct, punct = 2 * Zc, floor((K - K_ap) // Zc) * Zc
    A = K - colpunct - punct
    llr = 2/p**2 if channel== 'AWGN' else log((1-p)/p)
    col_inf = [0] * (2 * Zc)
    if channel == 'AWGN':
        inf_bits = list(r[:A]*llr)
        punct_inf = [-oo]*punct
        p_bits = list(r[A:]*llr)
    elif channel == 'BSC':
        inf_bits = list(r[:A]*llr)
        punct_inf = [-oo]*punct
        p_bits = list(r[A:]*llr)
    else:   # channel = BEC
        inf_bits = list(map(lambda y: 0 if y == 2 else y*oo, r[:A]))
        punct_inf = [-oo] * punct
        p_bits = list(map(lambda y: 0 if y == 2 else y*oo, r[A:]))
    #te = vector(RealField(19),punct_inf + list(inf_bits) + list(short_bits) + list(p_bits))
    return vector(RealField(10), col_inf + inf_bits + punct_inf + p_bits)
