# cd Desktop/Thesis/PySageMath/LDPC
from sage.all import *


def nz_sum_approx(Lv, min_vals, rcore):
    output = []
    for i, row in enumerate(Lv):
        vec = vector(RealField(10), (row.values()))
        sign = (-1)**(len(vec))
        #signs = [-1 if vec[j] < 0 else 1 for j in range(len(vec))]
        signs = [-1 if elem < 0 else 1 for elem in row]
        P = product(signs)
        Lvi_signs = list([a * b for a, b in zip(signs, [P*sign] * len(signs))])
        Lvi = nz_update_row(vec, min_vals[i], Lvi_signs, rcore) #TODO: 'i < rcore'
        output.append({k:v for k,v in zip(row.keys(), Lvi)})
    return output


def it_nz_sum_approx(input_vec, min_vals, rcore, lam, gamma):
    vec = vector(RealField(10), (input_vec.values()))
    sign = (-1) ** (len(vec))
    #breakpoint()
    #signs = [sgn(a) for a in vec]
    signs = [-1 if a < 0 else 1 for a in vec]
    P = product(signs)
    Lvi_signs = list([a * b for a, b in zip(signs, [P * sign] * len(signs))])
    Lvi = nz_update_row(vec, min_vals, Lvi_signs, rcore, lam, gamma)
    return {k:v for k,v in zip(input_vec.keys(), Lvi)}


def nz_update_row(row, min_vals, signs, offset, lam, gamma):
    output, min_get = [0] * len(row), len(min_vals)
    if len(min_vals) == 1:
        min_vals = vector(RealField(10), list(min_vals)*2)
    if len(min_vals) == 0:
        min_vals = [0.0001, 0.0001]
    #breakpoint()
    for i, elem in enumerate(row):
        if abs(elem) == oo:
            output[i] = elem
        elif min_vals[0] != abs(elem):
            output[i] = signs[i] * max(gamma * min_vals[0] - (lam * offset), 0)
        else:
            output[i] = signs[i] * max(gamma * min_vals[-1] - (lam * offset), 0)
    return output


def nz_min_fun(Lc, n_mins):
    mins = []
    for i, row in enumerate(Lc):
        vec = vector(RealField(10), (row.values()))
        non_zeros = [abs(vec[a]) for a in vec.nonzero_positions()]
        minis = sorted(non_zeros)[:n_mins]
        mins.append(vector(RealField(10), minis))
    return mins


def vec_mins(input_vec, n_mins):
    vec = vector(RealField(10), (input_vec.values()))
    non_zeros = [abs(vec[a]) for a in vec.nonzero_positions()]
    return vector(RealField(10), sorted(non_zeros)[:n_mins])


def nz_col_sum(NZMatrix, N):
    output = [0]*N
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
def minsum_SPA(H, HNZ, llr, r, rcore, lam, gamma, Zc, K, N0, use_core):
    lv = [{} for i in range(len(HNZ))]
    Lc = [{} for i in range(len(HNZ))]
    for i, row in enumerate(HNZ):
        for j, elem in row.items():
            lv[i].update({j: llr[j]})

    #r = vector(RealField(10), [a/(4/N0) for a in llr])
    codeword, runs = False, 0
    #breakpoint()
    while not codeword and runs < 20:
        for l in range(len(HNZ)):
            min_vals = vec_mins(lv[l], 2)
            Lc[l] = it_nz_sum_approx(lv[l], min_vals, True, lam, gamma) # l < rcore
        ltot = nz_col_sum(Lc, len(llr)) + llr
        vhat = vector(GF(2), [0 if elem <= 0 else 1 for elem in ltot])
        runs += 1
        # check if v_hat is a valid codeword
        if H.matrix_from_rows_and_columns(list(range(4*Zc)), list(range(K + 4*Zc))) * vhat[:K+4*Zc] == 0:
            #print(f"MinSum runs := {runs}")
            return vhat, True, runs

        # update Lv
        colvecs = get_column_vectors(Lc, len(llr))
        for j, col in enumerate(colvecs):
            col_sum = sum(col.values())

            for i, elem in col.items():
                lv[i].update({j: sign(col_sum)*(abs(col_sum) - col.get(i)) + llr[j]})
    return vhat, False, runs

