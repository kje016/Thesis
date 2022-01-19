from sage.all import *
import CRC

"From standard 6.2.2"
def det_BG(A, R):
    if A <= 292 or (A <= 3824 and R <= 0.67) or R <= 0.25:
        return 2
    else:
        return 1


def get_code_block_param(bg, B):
    Kcb = 8448 if bg==1 else 3840
    if B <= Kcb:
        L, C, B_ap = 0, 1, B
    else:
        L = 24
        C = ceil(B/(Kcb-L))
        B_ap = B + C * L
    return L, C, B_ap


def determine_kb(B, bg):
    if bg == 1:
        kb = 22
    else:
        if B > 640:
            kb = 10
        elif B > 560:
            kb = 9
        elif B > 192:
            kb = 8
        else:
            kb = 6
    return kb


# C := Input to channel coding
# D := Bits after encoding
def get_d_c(Zc, K, C):
    D = []
    for k in range(2*Zc, K):
        if C[k] != None:
            D.append(C[k])
        else:
            C[k] = 0
            D.append(None)
    return D, C


def det_Z(bg, kb, lifting_set, K_ap):
    lifted_set = {}
    for key, value in lifting_set.items():
        lifted_set.update({key: min([x for x in value if kb*x >= K_ap], default=oo)})
    lifting_size = min(lifted_set.values())
    set_index = list(lifted_set.values()).index(lifting_size)
    K = 22 * lifting_size if bg == 1 else 10 * lifting_size
    return lifting_size, set_index, K


def calc_crk(C, K_ap, K, L, b_bits):
    s, output = 0, []
    for r in range(C):
        tess = b_bits[r*(K_ap-L): (r+1)*(K_ap-L)]
        if C > 1:
            tess.extend(CRC.main_CRC(tess, CRC.crc24b))
        tess.extend([0]*(K-K_ap))
        output.append(tess)
    if len(output) == 1:
        return output[0]
    return output
