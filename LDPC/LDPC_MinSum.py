# cd Desktop/INF244/Exercises/MA3
from sage.all import *

def min_update_row(row, min_vals, signs, H_row):
    output, min_get = [0]*len(row), len(min_vals)
    scale, offset = 1, 0.4
    for i, elem in enumerate(H_row.nonzero_positions()):
        if min_vals[-min_get] == abs(row[elem]):
            output[elem] = signs[i] * max(scale*min_vals[min_get-1]+offset, 0)
        else:
            output[elem] = signs[i] * max(scale*min_vals[-min_get]+offset, 0)
    return output


def sum_approx(Lv, min_vals, H):
    output = []
    for i, row in enumerate(Lv):
        sign = (-1)**(len(H[i].nonzero_positions()))
        signs = [-1 if row[pos] < 0 else 1 for pos in H[i].nonzero_positions()]
        P = product(signs)
        Lvi_signs = list([a*b for a, b in zip(signs, [P*sign]*len(signs))])
        Lvi = min_update_row(row, min_vals[i], Lvi_signs, H.row(i))
        output.append(Lvi)
    return Matrix(output)


def min_fun(Lc, n_mins):
    mins = []
    for i, row in enumerate(Lc):
        non_zeros = list(map(lambda a: abs(row[a]), row.nonzero_positions()))
        minis = sorted(non_zeros[:n_mins])
        for j in range(n_mins, len(non_zeros)):
            if non_zeros[j] < minis[-1]:
                minis[-1] = non_zeros[j]
                minis.sort()
        mins.append(minis)
    return mins


def nz_sum_approx(Lv, min_vals, rcore):
    output = []
    for i, row in enumerate(Lv):
        vec = vector(RealField(10), (row.values()))
        sign = (-1)**(len(vec))
        signs = [-1 if vec[j] < 0 else 1 for j in range(len(vec))]
        P = product(signs)
        Lvi_signs = list([a * b for a, b in zip(signs, [P * sign] * len(signs))])
        slvi = [product(signs[:j])*product(signs[j+1:]) for j in range(len(signs))]
        Lvi = nz_update_row(vec, min_vals[i], Lvi_signs, i < rcore)
        output.append({k:v for k,v in zip(row.keys(), Lvi)})
    return output


def nz_update_row(row, min_vals, signs, offset):
    output, min_get = [0] * len(row), len(min_vals)
    if len(min_vals) == 1:
        min_get.vals(min_get[0])
    for i, elem in enumerate(row):
        if min_vals[-min_get] == abs(elem):
            output[i] = signs[i] * max(min_vals[min_get - 1] - (1 * offset), 0)
        else:
            output[i] = signs[i] * max(min_vals[-min_get] - (1 * offset), 0)
    return output


def nz_min_fun(Lc, n_mins):
    mins = []
    for i, row in enumerate(Lc):
        vec = vector(RealField(10), (row.values()))
        non_zeros = [abs(vec[a]) for a in vec.nonzero_positions()]
        minis = sorted(non_zeros[:n_mins])
        for j in range(n_mins, len(non_zeros)):
            if non_zeros[j] < minis[-1]:
                minis[-1] = non_zeros[j]
                minis.sort()
        mins.append(vector(RealField(10), minis))
    return mins


def nz_col_sum(NZMatrix, r):
    output = [0]*len(r)
    for index, row in enumerate(NZMatrix):
        for j, elem in row.items():
            output[j] += elem
    return vector(RealField(10), output)


def get_column_vectors(nzmatrix, length):
    output = [{} for i in range(length)]
    for index, row in enumerate(nzmatrix):
        for j, elem in row.items():
            output[j].update({index: elem})
    return output


# Lj = [(4*sqrt(Ec)/N0)*r[j] for j in range(len(r))] # (4*sqrt(Ec)/N0)*r[j] = 1*r[j] = r[:] in this case
def minsum_SPA(H, r, channel, sigma, rcore):
    if channel == 'AWGN':
        Lj = [(2/sigma)*rj for rj in r]
    else:
        Lj = list(r)
    lv = []
    for i in range(H.nrows()):
        temp = {}
        for j in H.row(i).nonzero_positions():
            temp.update({j: Lj[j]})
        lv.append(temp)
    codeword, runs = False, 0
    # breakpoint()
    while not codeword and runs < 30:
        min_vals = nz_min_fun(lv, 2)
        Lc = nz_sum_approx(lv, min_vals, rcore)
        ltot = nz_col_sum(Lc, Lj) + r
        vhat = vector(GF(2), [0 if elem <= 0 else 1 for elem in ltot])
        runs += 1
        # check if v_hat is a valid codeword
        if H * vhat == 0:
            #print(f"MinSum runs := {runs}")
            return vhat, True, runs

        # update Lv
        colvecs = get_column_vectors(Lc, len(r))
        for j, col in enumerate(colvecs):
            col_sum = sum(col.values())

            for i, elem in col.items():
                lv[i].update({j: col_sum - col.get(i) + Lj[j]})
    return vhat, False, runs

