# cd Desktop/Thesis/PySageMath/LDPC
from sage.all import *

import LDPC_MinSum


def sign(x):
    return 1 if x >= 0 else -1


def closest_to(elem, list):
    closest = 0
    for i, a in enumerate(list[1:]):
        if abs(abs(elem)-a) < abs(abs(elem)-list[closest]):
            closest = i+1
        else:
            break
    return sign(elem)*list[closest]


def quantization_map(y, mu, A):
    gamma = []
    for i, elem in enumerate(y):
        gamma.append(closest_to(mu*elem, A))
    return gamma


def fun1(gamma, beta, alphas, pos):
    temp = alphas.copy()
    for i, row in enumerate(alphas):
        for k, v in row.items():
            temp[i].update({k: gamma[k] - beta[k].get(i + pos)})
    return temp


def fun2(alpha, l, R, min_vals):
    output = []
    for i, row in enumerate(alpha):
        vec = vector(RealField(10), row.values())
        sign = (-1)**(len(vec))
        signs = [-1 if vec[j] < 0 else 1 for j in range(len(vec))]
        P = product(signs)
        Lvi_signs = list([a * b for a, b in zip(signs, [P * sign] * len(signs))])
        Lvi = nz_update_row(vec, min_vals[i], Lvi_signs, l < 4)
        output.append({k:v for k,v in zip(row.keys(), Lvi)})
    return output


def nz_update_row(row, min_vals, signs, offset):
    output, min_get = [0] * len(row), len(min_vals)
    for i, elem in enumerate(row):
        if min_vals[-min_get] == abs(elem):
            output[i] = signs[i] * max(min_vals[min_get-1]-(1*offset), 0)
        else:
            output[i] = signs[i] * max(min_vals[-min_get]-(1*offset), 0)
    return output

# ensures messages are sent of q bits
def saturation_fun(alpha, Q):
    temp = alpha.copy()
    for i, elem in enumerate(alpha):
        for k,v in elem.items():
            temp[i].update({k: sign(v)*min(abs(v), Q)})
    return temp


def OMS(BG, Zc, H, r, channel, sigma):
    M, N = H.nrows(), H.ncols()
    R, C = M//Zc, N//Zc
    L = R
    gamma = list(r)
    """ Non-zero matrix"""
    alphas = []
    for i in range(H.nrows()):
        temp = {}
        for j in H.row(i).nonzero_positions():
            temp.update({j: gamma[j]})
        alphas.append(temp)
    beta = []
    for i in range(H.ncols()):
        temp = {}
        for j in H.column(i).nonzero_positions():
            temp.update({j : 0})
        beta.append(temp)
    next_beta = beta.copy()

    """Initialization"""
    q, qtilde = 4, 6
    Q, Qtilde = 2**(q-1)-1, 2**(qtilde-1)-1
    A, Atilde = list(range(-Q, Q+1)), list(range(-Qtilde, Qtilde+1))
    mu = 3.8
    gamma = list(r) #quantization_map(r, mu, list(range(0, Q+1)))

    runs = 0
    breakpoint()
    while runs < 30:
        for l in range(L):
            # Ml = list(range(l*(M/L)+1, (l+1)*(M/L)))
            min_vals = LDPC_MinSum.nz_min_fun(alphas[l*Zc:(l+1)*Zc], 2)
            alpha = fun1(gamma, beta, alphas[l*Zc: (l+1)*Zc], l*Zc)
            alphas[l*Zc: (l+1)*Zc] = alpha

            # tess1 = saturation_fun(tess, Q)
            next_beta[l*Zc: (l+1)*Zc] = fun2(alpha, l, R, min_vals)
        next_beta = LDPC_MinSum.get_column_vectors(next_beta, len(r))
        for l in range(L):
            colvecs = [list(a.values()) for a in beta[l*Zc: (l+1)*Zc]]
            for j, col in enumerate(colvecs):
                gamma[j] = sum(colvecs[j]) + gamma[j+(l*Zc)]
                    #alphas[i].update({j: col_sum - col.get(i) + gamma[j]})
        vhat = vector(GF(2), [0 if elem<= 0 else 1 for elem in gamma])
        runs += 1
        print(H*vhat)
        if H*vhat == 0:
            print(f"MinSum runs := {runs}")
            return vhat, True
        beta = next_beta.copy()
    print(f"failed")
    breakpoint()
    return vhat, False
