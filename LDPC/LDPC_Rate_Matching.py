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
    return e, Hm


# punctured bits are unknown
# shortened filler bits are known to be 0
def fill_e(r, Zc, K, K_ap, p, Kb, channel):
    punct, short = 2 * Zc, floor((K - K_ap) // Zc) * Zc
    A = Kb - short - punct
    llr = log((1-p)/p)

    if channel == 'AWGN':
        punct_inf = LLR_fun([0.5]*(2*Zc), 'BSC', 0.1)
        inf_bits = r[:A]
        short_bits = vector(RealField(10), [-1]*floor((K-K_ap)//Zc) * Zc)
        p_bits = r[A:]
    else:
        punct_inf = LLR_fun([0.5]*(2*Zc), 'BSC', 0.1)
        inf_bits = LLR_fun(r[:A], 'BSC', 0.1)
        short_bits = [llr]*floor((K-K_ap)//Zc) * Zc
        p_bits = LLR_fun(r[A:], 'BSC', 0.1)
    return vector(RealField(10), list(punct_inf) + list(inf_bits) + list(short_bits) + list(p_bits))


# p_cross gjelder egt bare for BSC. evnt bytte til channel_metric/channel_characteristic
def LLR_fun(e, channel, p_cross):
    F = RealField(10)
    if channel == 'BSC':
        llr1 = log(p_cross/(1-p_cross))
        return vector(F, list(map(lambda x: x - 1, 2 * vector(F, e)))) * llr1