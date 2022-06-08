# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions as HF

node_states = ['l', 'r', 'u']
F = RealField(10)


def decoder(d, N, frozen_set, p_cross, I_IL, PI, C):
    C_perm, CRCpos = [], []
    if I_IL:
        C_perm = Matrix(C[a] for a in [a for a in PI if a < C.nrows()])
        CRCpos = [PI.index(a) for a in PI if a >= C.nrows()]
    llr = -p_cross
    tree = HF.init_tree(N, d)
    """SCL initialization"""
    L = 8
    list_decoders = [HF.Decoder("", 0)]
    depth, done, node = 0, False, tree[0]
    node_i = tree.index(node)

    while not done:
        if depth == log(N, 2):
            node.state = node_states[2]
            is_frozen = node_i-(N-1) in frozen_set # alternatively var name,
            list_decoders = HF.bec_update_decoders(is_frozen, node.beliefs[0], llr,  list_decoders, L, C_perm, CRCpos)
            if not list_decoders:   # no surviving paths
                return [HF.Decoder('', +oo)]
            node.beliefs = HF.bec_uhat(node.beliefs, is_frozen)
            if tree.index(node) == len(tree)-1:
                done = True
            node.state = node_states[2]
            node_i = floor(abs((node_i - 1)) / 2)
            node, depth = tree[node_i], depth - 1
        elif node.state == "":  # step L
            tree[2 * node_i + 1].beliefs = HF.f_bec(node.beliefs)
            node.state = node_states[0]
            node_i = 2 * node_i + 1
            node, depth = tree[node_i], depth+1
        elif node.state == node_states[0]:  # step R
            tree[2 * node_i + 2].beliefs = HF.g_bec(node.beliefs, tree[2 * node_i + 1].beliefs)
            node.state = node_states[1]
            node_i = 2 * node_i + 2
            node, depth = tree[node_i], depth+1
        else:   # step U
            node.beliefs = vector(F, list(map(lambda x: HF.bec_xor(x[0], x[1]), zip(tree[2 * node_i + 1].beliefs,
                                                    tree[2 * node_i + 2].beliefs))) + list(tree[2 * node_i + 2].beliefs))
            node.state = node_states[2]
            node_i = floor((node_i - 1) / 2)
            node, depth = tree[node_i], depth-1

    if PI:
        for dec in list_decoders:
            deinterleave = [0]*len(dec.inf_bits)
            for i, elem in enumerate(PI):
                deinterleave[elem] = dec.inf_bits[i]
            dec.inf_bits = vector(GF(2), deinterleave[:C.nrows()])
    else:
        for dec in list_decoders:
            dec.inf_bits = vector(GF(2), dec.inf_bits)
    return list_decoders
