# cd Desktop/INF244/Exercises/MA3
from sage.all import *


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
    output = vector(RealField(10), [0 if a == 0 else a*oo for a in output])
    return vector(RealField(10), output)


def get_column_vectors(nzmatrix, length):
    output = [{} for i in range(length)]
    for index, row in enumerate(nzmatrix):
        for j, elem in row.items():
            output[j].update({index: sign(elem)})
    return output


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
    LTOTS= []
    while not codeword and runs < 20:
        Lc = nz_tanh(lv)
        ltot = nz_col_sum(Lc, Lj) + vector(map(lambda a: sgn(a), r))
        vhat = vector(GF(2), [0 if elem <= 0 else 1 for elem in ltot])
        runs += 1
        #ltot = vector(RealField(10), [a * oo if a != 0 else 0 for a in vector(RealField(10), nz_col_sum(Lc, Lj))])
        if H * vhat == 0:    # check if v_hat is a valid codeword
            #print(f"MinSum runs := {runs, H*vhat == 0}")
            return vhat, True, runs
        if ltot in LTOTS:
            return vhat, codeword, runs

        # update Lv
        colvecs = get_column_vectors(Lc, len(r))
        for j, col in enumerate(colvecs):
            col_sum = sum(col.values())
            if r[j] == 0:
                if abs(col_sum) != 0:
                    for i, elem in col.items():
                        lv[i].update({j:  oo*(col_sum)})
        LTOTS.append(ltot)

    #print(f"MinSum runs := {runs, H*vhat == 0}")
    return vhat, codeword, runs

def PP(inputt):
    for elem in inputt:
        print(elem)