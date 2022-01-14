from sage.all import *



def RM_main(D, Zc, BG, H, K, K_ap, B, R):
    E = ceil(len(D)/Zc*R) * Zc
    Kb = H.ncols()-H.nrows()
    punct = 2*Zc
    short = floor((K-K_ap)//Zc) * Zc  # floor so it doesnt take more fille bits than there actullay is filler bits
    A = Kb - short - punct

    pbits = D[Kb: Kb + (E-A)]
    e = list(D[2*Zc: A + 2*Zc]) + list(pbits)

    Hm = H.matrix_from_rows_and_columns(list(range(E-A)), list(range(Kb + (E-A))))
    ee = fill_e(e, Zc, K, K_ap, A)
    #pbits = D[B+ (K-K_ap):B + (K-K_ap) + (E-A+2*Zc)]

    breakpoint()


def fill_e(e, Zc, K, K_ap, A):
    res = [1]*(2*Zc) + e[:A] + [0]*(floor((K-K_ap)//Zc) * Zc) + e[A:]
    return vector(GF(2), res)