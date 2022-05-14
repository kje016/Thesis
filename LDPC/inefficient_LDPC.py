# cd Desktop/INF244/Exercises/MA3
from sage.all import *


def min_update_row(row, min_vals, signs):
    # Should not happen that we only have one min_value
    output = []
    for i, elem in enumerate(row):
        if elem != 0:
            if min_vals[0] == abs(elem):
                if len(min_vals) == 1:
                    output.append(0)
                else:
                    output.append(min_vals[1]*signs[i])
            else:
                output.append(min_vals[-len(min_vals)]*signs[i])
        else:
            output.append(0)
    return output


def extrensic_information(Lv, min_vals):
    output = []
    for i, row in enumerate(Lv):
        sign = (-1)**(len(row.nonzero_positions()))
        signs = [-1 if elem < 0 else 1 for elem in row]
        P = product(signs)
        Lvi_signs = list([a*b for a, b in zip(signs, [P*sign]*len(row))])
        Lvi = min_update_row(row, min_vals[i], Lvi_signs)
        output.append(Lvi)
    return Matrix(output)


def comp_l_tot(lc, r):
    return vector([sum(lc.column(i))+r[i] for i in range(len(r))])


def min_fun(Lc, n_mins):
    mins = []
    for i, row in enumerate(Lc):
        mini = [abs(elem) for elem in row if abs(elem) > 0]
        mins.append(sorted(mini)[:n_mins])
    return mins


def minsum_SPA(r, H, N0):
    Lj = [((4*1)/N0)*r[j]for j in range(len(r))] # (4*sqrt(Ec)/N0)*r[j] = 1*r[j] = r[:] in this case
    Lv = [RealNumber(elem)*Lj[j] for i in range(H.nrows()) for j, elem in enumerate(H.row(i))]
    Lv = Matrix(RealField(10), H.nrows(), H.ncols(), Lv)

    codeword, runs = False, 0
    while not codeword and runs < 20:
        min_vals = min_fun(Lv, 2)
        Lc = extrensic_information(Lv=Lv, min_vals=min_vals)
        l_tot = comp_l_tot(Lc, r)
        v_hat = vector(GF(2), [0 if elem <= 0 else 1 for elem in l_tot])
        runs += 1

        # check if v_hat is a valid codeword
        if H * v_hat == 0:
            return v_hat, True

        # update Lv
        for i, row in enumerate(Lc.rows()):
            for j, col in enumerate(Lc.columns()):
                if H[i, j] != 0:
                    col_res = list(col)[0:i] + list(col)[i + 1:len(col)]
                    Lv[i, j] = Lj[j] + sum(col_res)  # Lj[j]
    return v_hat, codeword



