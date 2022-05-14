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


def nz_sum_approx(Lv):
    output = []
    breakpoint()
    for i, row in enumerate(Lv):
        vec = vector(RealField(10), (row.values()))
        sign = (-1)**(len(vec))
        signs = [-1 if vec[j] < 0 else 1 for j in range(len(vec))]
        P = product(signs)
        Lvi_signs = list([a * b for a, b in zip(signs, [P * sign] * len(signs))])
        Lvi = [Lvi_signs[j]*oo for j in range(len(vec))]
        output.append({k: v for k, v in zip(row.keys(), Lvi)})
    return output


def nz_tanh(Lv):
    output = []
    for i, row in enumerate(Lv):
        veci = vector(RealField(10), list(map(lambda a: sign(a), row.values())))
        Lvi = []
        if len(veci.nonzero_positions()) - len(veci) >= 2:
            Lvi = [0.0 for a in range(len(veci))]
        else:
            for j, elem in enumerate(veci):
                prodtemp = product(veci[:j])*product(veci[j+1:])
                Lvi.append(0 if prodtemp == 0 else prodtemp*oo)
        output.append({k: v for k, v in zip(row.keys(), Lvi)})
    return output


def nz_col_sum(NZMatrix, r):
    output = r[:]
    for index, row in enumerate(NZMatrix):
        col_pos = list(row.keys())
        nzelems = [col_pos[j] for j, elem in enumerate(col_pos) if r[elem] == 0]
        for elem in nzelems:
            output[elem] += sgn(row.get(elem))
    return vector(RealField(10), output)


def comp_l_tot(lc):
    res = [sum(map(lambda a: sgn(a), lc.column(i))) for i in range(lc.ncols())]
    return vector(RealField(10), list(map(lambda a: 0 if a == 0 else a*oo, res)))


def get_column_vectors(nzmatrix, length):
    output = [{} for i in range(length)]
    for index, row in enumerate(nzmatrix):
        for j, elem in row.items():
            output[j].update({index: sign(elem)})
    return output


def min_fun(Lc, n_mins):
    mins = []
    for i, row in enumerate(Lc):
        mini = [abs(elem) for elem in row if abs(elem) > 0]
        mins.append(sorted(mini)[:n_mins])
    return mins


# Lj = [(4*sqrt(Ec)/N0)*r[j] for j in range(len(r))] # (4*sqrt(Ec)/N0)*r[j] = 1*r[j] = r[:] in this case
def minsum_BEC(H, r):
    Lj = r
    lv = []
    for i in range(H.nrows()):
        temp = {}
        for j in H.row(i).nonzero_positions():
            temp.update({j: Lj[j]})
        lv.append(temp)
    codeword, runs = False, 0
    while not codeword and runs < 30:
        #Lc = nz_sum_approx(lv)
        Lc = nz_tanh(lv)
        ltot = nz_col_sum(Lc, Lj) + vector(map(lambda a: sgn(a), r))
        vhat = vector(GF(2), [0 if elem <= 0 else 1 for elem in ltot])

        runs += 1
        ltot = vector(RealField(10), [a * oo if a != 0 else 0 for a in vector(RealField(10), nz_col_sum(Lc, Lj))])
        testvec = [a if a != 0 else 2 for a in ltot]
        # check if v_hat is a valid codeword
        if H * vhat == 0:
            #print(f"MinSum runs := {runs, H*vhat == 0}")
            return vhat, True, runs

        # update Lv
        colvecs = get_column_vectors(Lc, len(r))
        for j, col in enumerate(colvecs):
            if r[j] == 0:
                col_sum = sum(col.values())
                if abs(col_sum) != 0:
                    for i, elem in col.items():
                        lv[i].update({j:  oo*(col_sum)})

    #print(f"MinSum runs := {runs, H*vhat == 0}")
    return vhat, codeword, runs

