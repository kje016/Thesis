# cd Desktop/Thesis/PySageMath/PC
import gc

from sage.all import *
import HelperFunctions as HF

node_states = ['l', 'r', 'u']
F = RealField(7)


def decoder(d, N, frozen_set):
    tree = HF.init_tree(N, d)
    depth, done, node = 0, False, tree[0]
    vhat = []
    while not done:
        if depth == log(N, 2):
            node.state = node_states[2]
            is_frozen = tree.index(node)-(2**log(N, 2)-1) in frozen_set
            node.beliefs = HF.uhat(node.beliefs, is_frozen, F)
            if not is_frozen:
                vhat.append(node.beliefs[0])
            if tree.index(node) == len(tree)-1:
                done = True
            node, depth = tree[floor(abs((tree.index(node)-1))/2)], depth-1
        elif node.state == "":  # step L
            tree[node.l_child].beliefs = HF.ft(node.beliefs, F)
            node.state = node_states[0]
            node, depth = tree[node.l_child], depth+1
        elif node.state == node_states[0]:  # step R
            tree[node.r_child].beliefs = HF.gt(node.beliefs, tree[node.l_child].beliefs, F)
            node.state = node_states[1]
            node, depth = tree[node.r_child], depth+1
        else:   # step U
            node.beliefs = vector(F, list(
                map(lambda x: mod(x[0] + x[1], 2), zip(tree[node.l_child].beliefs, tree[node.r_child].beliefs))) + list(
                tree[node.r_child].beliefs))
            node.state = node_states[2]
            node, depth = tree[floor((tree.index(node)-1)/2)], depth-1
    breakpoint()
    del tree
    gc.collect()
    return vector(GF(2), vhat)
