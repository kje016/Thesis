from sage.all import *


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
