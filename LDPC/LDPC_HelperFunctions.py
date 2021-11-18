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


# Z = a * 2**j
def determine_Z(input_kb, lifting_size_set, inf_bits):
    lifted_set = {}
    for key, value in lifting_size_set.items():
        crnt_lift = min([x for x in value if input_kb*x >= inf_bits])
        lifted_set.update({key: crnt_lift})
    lifting_size = min(lifted_set.values())
    set_index = list(lifted_set.values()).index(lifting_size)
    return lifting_size, set_index


def get_base_matrix(bg, ils, zc):
    matrix = []
    f = open("base_matrices\\" + f"NR_{bg}_{ils}_{zc}.txt", "r")
    txt_m = f.read().split('\n')
    for row in txt_m:
        row = list(map(int, [c for c in row.split(' ') if c!= '']))
        matrix.append(row)
    matrix.pop(-1)          # matrix contains an empty list at the end
    return Matrix(matrix)

