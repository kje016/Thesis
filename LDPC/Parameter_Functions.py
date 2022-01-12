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


def det_Z(bg, kb, lifting_set, K_ap):
    lifted_set = {}
    for key, value in lifting_set.items():
        lifted_set.update({key: min([x for x in value if kb*x >= K_ap], default=oo)})
    lifting_size = min(lifted_set.values())
    set_index = list(lifted_set.values()).index(lifting_size)
    K = 22 * lifting_size if bg == 1 else 10 * lifting_size
    return lifting_size, set_index, K


def calc_crk(C, K_ap, K, L, b_bits):
    s, crk = 0, []
    for r in range(C):
        for k in range(K_ap-L):
            crk.append(b_bits[s])
            s += 1
        if C > 1:
            print("calc_crk() not finished for C > 1")
            crk.append(CRC.main_CRC(b_bits[r*K_ap-L: (r+1)*K_ap-L], CRC.crc24b))
            for elem, a in enumerate(crk):
                crk[a] = crk[a].extend([0]*len(crk[a])-K_ap)
            print("calc_crk() need implementing when C>1")
    crk.extend([0]*(K-K_ap))
    return crk
