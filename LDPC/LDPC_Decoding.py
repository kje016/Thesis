# cd Desktop/INF244/Exercises/MA3
from sage.all import *


def min_update_row(row, min_vals, signs, H_row):
    # Should not happen that we only have one min_value
    output = []
    for i, elem in enumerate(row):
        if H_row[i] != 0:
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


def extrensic_information(Lv, min_vals, H):
    output = []
    for i, row in enumerate(Lv):
        sign = (-1)**(len(H[i].nonzero_positions()))
        signs = [-1 if elem < 0 else 1 for elem in row]
        P = product(signs)
        Lvi_signs = list([a*b for a, b in zip(signs, [P*sign]*len(row))])
        Lvi = min_update_row(row, min_vals[i], Lvi_signs, H[i])
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


# Lj = [(4*sqrt(Ec)/N0)*r[j] for j in range(len(r))] # (4*sqrt(Ec)/N0)*r[j] = 1*r[j] = r[:] in this case
def minsum_SPA(H, r, N0, channel, sigma):
    if channel == 'AWGN':
        Lj = [(2/sigma)*rj for rj in r]
    else:
        Lj = list(r)

    Lv = [RealNumber(elem)*Lj[j] for i in range(H.nrows()) for j, elem in enumerate(H.row(i))]
    Lv = Matrix(RealField(10), H.nrows(), H.ncols(), Lv)
    codeword, runs = False, 0
    while not codeword and runs < 30:
        min_vals = min_fun(Lv, 2)
        Lc = extrensic_information(Lv=Lv, min_vals=min_vals, H=H)
        l_tot = comp_l_tot(Lc, r)
        v_hat = vector(GF(2), [0 if elem <= 0 else 1 for elem in l_tot])
        runs += 1

        # check if v_hat is a valid codeword
        if H * v_hat == 0:
            print(f"MinSum runs := {runs}")
            return v_hat, True

        # update Lv
        for i, row in enumerate(Lc.rows()):
            for j, col in enumerate(Lc.columns()):
                if H[i, j] != 0:
                    col_res = list(col)[0:i] + list(col)[i + 1:len(col)]
                    Lv[i, j] = Lj[j] + sum(col_res)  # Lj[j]
    print(f"    MinSum runs := {runs}")
    return v_hat, codeword


def spa_main(H, r, N0, channel, sigma):
    # RDF = RealDoubleField(), the elements are of double precision floating numbers
    v_hat, is_codeword = minsum_SPA(H, r, N0, channel, sigma)
    #m_hat = vector(GF(2), list(v_hat)[:G.nrows()])
    return v_hat, is_codeword

"""
K = 80
mb = 42
nb = 52
kb = 10
Zc = 8
Am = H.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(K//2)))
LVA = Lv.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(K//2)))
LCA = Lc.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(K//2)))


Bm = H.matrix_from_rows_and_columns(list(range(4 * Zc)), list(range(K, K + 4*Zc)))
Cm = H.matrix_from_rows_and_columns(list(range((4*Zc, H.nrows())), list(range(kb * Zc)))
Dm = H.matrix_from_rows_and_columns(list(range((mb - 4) * Zc)), list(range(kb * Zc, (kb * Zc) + 4 * Zc)))
"""