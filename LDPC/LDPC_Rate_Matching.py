from sage.all import *


def RM_main(u, Zc, H, K, K_ap, rate, B, channel):
    colpunct, punct = 0 if channel == 'BEC' else 2 * Zc, floor((K - K_ap) // Zc) * Zc
    colpunct= 0# TODO
    A_ap = K - colpunct - punct   # A is the amount of crc bits after removing 2 cols and padding
    E = ceil((B / rate) / Zc) * Zc
    pbits = u[K: K + (E-A_ap)]   # getting the parity bits and
    e = list(u[colpunct: B+(K-K_ap-punct)]) + list(pbits)
    Hm = H.matrix_from_rows_and_columns(list(range(E-A_ap)), list(range(K + (E-A_ap))))
    return vector(RealField(10), e), Hm


# punctured bits are unknown
# shortened filler bits are known to be 0
def fill_w_llr(r, Zc, K, K_ap, p, channel):
    colpunct, punct = 0 if channel == 'BEC' else 2 * Zc, floor((K - K_ap) // Zc) * Zc
    colpunct = 0 # TODO
    A = K - colpunct - punct
    llr = -((4*1)/p) if channel== 'AWGN' else log((p/(1-p)))   #where p = N0
    col_inf = [0] * colpunct
    if channel == 'AWGN':
        inf_bits = list(r[:A]*llr)
        punct_inf = [oo]*punct
        p_bits = list(r[A:]*llr)
    elif channel == 'BSC':
        inf_bits = list(r[:A]*llr)
        punct_inf = [oo]*punct
        p_bits = list(r[A:]*llr)
    else:   # channel = BEC
        inf_bits = list(map(lambda y: 0 if y == 2 else y*-oo, r[:A]))
        punct_inf = [oo] * punct
        p_bits = list(map(lambda y: 0 if y == 2 else y*-oo, r[A:]))
    return vector(RealField(10), col_inf + inf_bits + punct_inf + p_bits)#, vector(RealField(10), col_inf + list(r[:A]) + [-oo]*punct + list(r[A:]))
