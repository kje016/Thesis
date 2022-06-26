from sage.all import *


def RM_main(u, Zc, H, K, K_ap, rate, B, channel):
    colpunct, punct = 0 if channel == 'BEC' else 2 * Zc, floor((K - K_ap) // Zc) * Zc
    #colpunct = 0# TODO
    A_ap = K - colpunct - punct   # A is the amount of crc bits after removing 2 cols and padding
    E = ceil((B / rate) / Zc) * Zc
    pbits = u[K: K + (E-A_ap)]   # getting the parity bits and
    e = list(u[colpunct: B+(K-K_ap-punct)]) + list(pbits)

    ev = vector(GF(2), list(u[:K+ 4*Zc]))
    Hm = H.matrix_from_rows_and_columns(list(range(E-A_ap)), list(range(K + (E-A_ap))))
    Ht = H.matrix_from_rows_and_columns(list(range(4*Zc)), list(range(K + 4*Zc)))
    if Ht*ev != 0:
        breakpoint()
    return vector(RealField(10), e), Hm


# punctured bits are unknown
# shortened filler bits are known to be 0
def fill_w_llr(r, Zc, K, K_ap, p, N0, channel, HRM):
    colpunct, punct = 0 if channel == 'BEC' else 2 * Zc, floor((K - K_ap) // Zc) * Zc
    #colpunct = 0# TODO
    A = K - colpunct - punct
    llr = (4/N0) if channel== 'AWGN' else log(((1-p) / p), 2)   #where p = N0
    col_inf = [0] * colpunct
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
    return vector(RealField(10), col_inf + inf_bits + punct_inf + p_bits)#, vector(RealField(10), col_inf + list(r[:A]) + [-oo]*punct + list(r[A:]))
