from sage.all import *



# B = len of information bits (A) + Crc-bits
def RM_main(D, Zc, BG, H, K, K_ap, B, R):
    E = ceil((len(D)*R))
    fill_num = K-K_ap - (K-K_ap)%Zc
    e = list(D[2*Zc:])
    #e = list(D[2*Zc:B]) + list(D[B+fill_num:])  # dropping first 2 inf-columns and filler bits of kb
    # e = list(D[2*Zc:E+(2*Zc)]) when only the first 2 cols are punctured
    e = e[:E]


    nbRM = (H.ncols() - E)//Zc + 2
    Kb = (H.ncols() - H.nrows()) // Zc
    mbRM = ceil(Kb / R) + 2  # + fill_num//Zc   # number of rows in H

    Hm = H.matrix_from_rows_and_columns(list(range(mbRM*Zc)), list(range(nbRM*Zc)))

    te = [1]*(2*Zc) + e
    ee = fill_E(e, B, Zc, K, K_ap)
    tess = list(D[:len(ee)])
    breakpoint()
    return e, Hm


def fill_E(e, B, Zc, K, K_ap):
    res = [1]*(2*Zc) + e[:B-(2*Zc)] + [0]*((K-K_ap - (K-K_ap)%Zc)) + e[B-(2*Zc):]
    return res


def ty(D, Zc, BG, H, K, K_ap, B, R):
    mb, nb = BG.nrows(), BG.ncols()
    kb = nb-mb
    k = kb*Zc
    nbRM = ceil(kb/R) + 2
    n = nbRM*Zc
    mbRM = nbRM - kb

    BGG = BG.matrix_from_rows_and_columns(list(range(mbRM)), list(range(nbRM)))
    Hm = H.matrix_from_rows_and_columns(list(range(mbRM*Zc)), list(range(nbRM*Zc)))
    breakpoint()


"""
BGG = BG.matrix_from_rows_and_columns(list(range(mbRM)), list(range(nbRM)))



[mb,bn] = size(B)
kb = nb-mb
R = 1/2
nbRM = ceil(kb/R) +2 + (K-Kap)
n = nbRM * Zc
mbRM = nbRM - kb
"""