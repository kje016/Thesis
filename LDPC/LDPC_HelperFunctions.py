from sage.all import *


def determine_kb_BG2(k_bits):
    if k_bits > 640:
        return 10
    elif (k_bits > 560 and k_bits <= 640):
        return 9
    elif (k_bits > 192 and k_bits <= 560):
        return 8
    else:
        return 6

def get_param(bg, B):
    Kcb = 8448 if bg==1 else 3840
    if B <= Kcb:
        L, C, B_ap = 0, 1, B
    else:
        L = 24
        C = ceil(B/(Kcb-L))
        B_ap = B + C * L
    return L, C, B_ap

def determine_kb(B, B_ap, C, bg):
    K_ap = B_ap//C    # number of bits in each code block
    kb = 0
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
    return K_ap, kb


def paramater_decision(bg, B, Zc):
    Kcb, L, C, K, B_ap = 0, 0, 0, 0, 0
    """ getting maximum code block size Kb"""
    if bg == 1:
        Kcb = 8448
    else:
        Kcb = 3840


    """ Total number of code blocks 'C' is determined by:    """
    if B <= Kcb:
        L, C, B_ap = 0, 1, B

    else:
        L = 24
        C = ceil(B/(Kcb-L))
        B_ap = B + C*L

    """ number of bits 'K' in each code block: """
    K_ap = B_ap//C

    """
    Find the minimum value of 'Z' denoted as Zc s.t Kb * Zc >= K_ap
    set K = 22*Zc for BG1
    set K = 10*Zc for BG2
    """
    if bg == 1:
        K = 22 * Zc
    else:
        K = 10 * Zc

    return Kcb, L, C, K, K_ap, B_ap


# Z = a * 2**j
def determine_Z(bg, input_kb, lifting_size_set, K_ap):
    lifted_set = {}
    for key, value in lifting_size_set.items():
        crnt_lift = min([x for x in value if input_kb*x >= K_ap])
        lifted_set.update({key: crnt_lift})
    lifting_size = min(lifted_set.values())
    set_index = list(lifted_set.values()).index(lifting_size)
    K = 22*lifting_size if bg == 1 else 10*lifting_size
    return lifting_size, set_index, K


def calc_crk(C, K_ap, K, L, b_bits):
    B= len(b_bits)
    s = 0
    crk = []
    for r in range(C):
        for k in range(K_ap-L):
            crk.append(b_bits[s])
            s += 1
        if C > 1:
            print("calc_crk() need implementing when C>1")
    crk.append([0]*(K_ap-(K-1)))
    return crk


def gen_bg(bg, Z):
    pass



def get_base_matrix(bg, ils, zc):
    matrix = []
    f = open("base_matrices\\" + f"NR_{bg}_{ils}_{zc}.txt", "r")
    txt_m = f.read().split('\n')
    for row in txt_m:
        row = list(map(int, [c for c in row.split(' ') if c!= '']))
        matrix.append(row)
    matrix.pop(-1)          # matrix contains an empty list at the end
    return Matrix(matrix)

