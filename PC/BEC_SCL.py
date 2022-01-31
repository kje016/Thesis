# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions as HF

node_states = ['l', 'r', 'u']
F = RealField(10)


def decoder(d, N, frozen_set, p_cross):
    llr = -p_cross
    tree = HF.init_tree(N, d)
    """SCL initialization"""
    L = 8
    list_decoders = [HF.Decoder("", 0)]
    depth, done, node = 0, False, tree[0]
    while not done:
        if depth == log(N, 2):
            node.state = node_states[2]
            is_frozen = tree.index(node)-(2**log(N, 2)-1) in frozen_set # alternatively var name,
            node.beliefs = HF.uhat(node.beliefs, is_frozen)
            list_decoders = HF.update_decoders(is_frozen, node.beliefs[0], llr,  list_decoders, L)
            if tree.index(node) == len(tree)-1:
                done = True
            node, depth = tree[floor(abs((tree.index(node)-1))/2)], depth-1
        elif node.state == "":  # step L
            tree[node.l_child].beliefs = HF.f_bec(node.beliefs)
            node.state = node_states[0]
            node, depth = tree[node.l_child], depth+1
        elif node.state == node_states[0]:  # step R
            tree[node.r_child].beliefs = HF.g_bec(node.beliefs, tree[node.l_child].beliefs)
            node.state = node_states[1]
            node, depth = tree[node.r_child], depth+1
        else:   # step U
            node.beliefs = vector(F, list(map(lambda x: HF.bec_xor(x[0], x[1]), zip(tree[node.l_child].beliefs, tree[node.r_child].beliefs))) + list(tree[node.r_child].beliefs))
            node.state = node_states[2]
            node, depth = tree[floor((tree.index(node)-1)/2)], depth-1

    for dec in list_decoders:
        bits = set(list(range(0, N))) - set(frozen_set)
        bits = [dec.inf_bits[a] for a in bits]
        dec.inf_bits = vector(GF(2), [str(a) for a in bits]) #''.join(a for a in bits)

    return list_decoders

