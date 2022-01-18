from sage.all import *


def RM_main(D, Zc, H, K, K_ap, R):
    E = ceil(len(D)/Zc*R) * Zc
    Kb, punct, short = H.ncols()-H.nrows(), 2*Zc, floor((K-K_ap)//Zc) * Zc
    # floor so it doesnt take more filler bits than there actullay is filler bits

    A = Kb - short - punct  # A is the amount of information bits after removing 2 cols and padding
    pbits = D[Kb: Kb + (E-A)]   # getting the parity bits and
    e = list(D[2*Zc: A + 2*Zc]) + list(pbits)

    Hm = H.matrix_from_rows_and_columns(list(range(E-A)), list(range(Kb + (E-A))))
    #ee = fill_e(e, Zc, K, K_ap, A)
    return vector(ZZ, e), Hm


# punctured bits are unknown
# shortened filler bits are known to be 0
def fill_e(r, Zc, K, K_ap, p, Kb, channel):
    punct, short = 2 * Zc, floor((K - K_ap) // Zc) * Zc
    A = Kb - short - punct
    llr = log((1-p)/p)
    if channel == 'AWGN':
        punct_inf = [0]*(2*Zc)
        inf_bits = r[:A]
        short_bits = [-1]*floor((K-K_ap)//Zc) * Zc
        p_bits = r[A:]
    elif channel == 'BSC':
        punct_inf = [0]*(2*Zc)
        inf_bits = map(lambda y: y*llr, r[:A])
        short_bits = [-llr]*floor((K-K_ap)//Zc) * Zc
        p_bits = map(lambda y: y*llr, r[A:]) #LLR_fun(r[A:], 'BSC', 0.1)
    else:   # channel = BEC
        punct_inf = [0] * (2 * Zc)
        inf_bits = map(lambda y: 0 if y == 2 else y*oo, r[:A])
        short_bits = [-oo] * (floor((K - K_ap) // Zc) * Zc)
        p_bits = map(lambda y: 0 if y == 2 else y*oo, r[A:])
    return vector(RealField(10), punct_inf + list(inf_bits) + list(short_bits) + list(p_bits))
