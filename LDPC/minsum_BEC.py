# cd Desktop/INF244/Exercises/MA3
from sage.all import *


def extrensic_information(Lv, H):
    output = []
    for i, row in enumerate(Lv):
        sign = (-1)**(len(H[i].nonzero_positions()))
        signs = [-1 if elem < 0 else 1 for elem in row]
        P = product(signs)
        Lvi_signs = list([a*b for a, b in zip(signs, [P*sign]*len(row))])
        Lvi = [0 if H.row(i)[j] == 0 else Lvi_signs[j]*oo for j in range(len(row))]
        output.append(Lvi)
    return Matrix(output)


def comp_l_tot(lc):
    res = [sum(map(lambda a: sgn(a), lc.column(i))) for i in range(lc.ncols())]
    return vector(RealField(10), list(map(lambda a: 0 if a == 0 else a*oo, res)))


def min_fun(Lc, n_mins):
    mins = []
    for i, row in enumerate(Lc):
        mini = [abs(elem) for elem in row if abs(elem) > 0]
        mins.append(sorted(mini)[:n_mins])
    return mins


# Lj = [(4*sqrt(Ec)/N0)*r[j] for j in range(len(r))] # (4*sqrt(Ec)/N0)*r[j] = 1*r[j] = r[:] in this case
def minsum_BEC(H, r):
    Lj = r
    Lv = [0 if elem == 0 else RealNumber(elem)*Lj[j] for i in range(H.nrows()) for j, elem in enumerate(H.row(i))]
    Lv = Matrix(H.nrows(), H.ncols(), Lv)
    codeword, runs = False, 0
    while not codeword and runs < 20:
        Lc = extrensic_information(Lv=Lv, H=H)
        l_tot = comp_l_tot(Lc)
        v_hat = vector(GF(2), [0 if elem <= 0 else 1 for elem in l_tot])
        runs += 1

        # check if v_hat is a valid codeword
        if H * v_hat == 0:
            print(f"MinSum runs := {runs}")
            return v_hat, True

        col_res = [list(map(lambda a: sgn(a), Lc.column(i))) for i in range(Lc.ncols())]
        for i, elem in enumerate(r):
            if r[i] == 0:
                for j in range(Lc.nrows()):
                    if H[j, i] == 1:
                        belief = sum(col_res[i][:j]+col_res[i][j+1:])
                        Lv[j, i] = belief

    print(f"    MinSum runs := {runs}")
    return v_hat, codeword

