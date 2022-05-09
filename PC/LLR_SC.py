# cd Desktop/Thesis/PySageMath/PC
import gc

from sage.all import *
import HelperFunctions as HF

node_states = ['l', 'r', 'u']
F = RealField(7)


def decoder(d, N, frozen_set):
    tree = HF.init_tree(N, d)
    depth, done, node = 0, False, tree[0]
    node_i = tree.index(node)
    vhat = []
    while not done:
        if depth == log(N, 2):
            is_frozen = node_i-(N-1) in frozen_set
            node.beliefs = HF.uhat(node.beliefs, is_frozen, F)
            if not is_frozen:
                vhat.append(node.beliefs[0])
            if node_i == len(tree)-1:
                done = True
            node.state = node_states[2]
            node_i = floor(abs((node_i-1))/2)
            node, depth = tree[node_i], depth-1
        elif node.state == "":  # step L
            tree[2*node_i+1].beliefs = HF.ft(node.beliefs, F)
            node.state = node_states[0]
            node_i = 2*node_i+1
            node, depth = tree[node_i], depth+1
        elif node.state == node_states[0]:  # step R
            tree[2*node_i+2].beliefs = HF.gt(node.beliefs, tree[2*node_i+1].beliefs, F)
            node.state = node_states[1]
            node_i = 2*node_i+2
            node, depth = tree[node_i], depth+1
        else:   # step U
            node.beliefs = vector(F, list(
                map(lambda x: mod(x[0] + x[1], 2), zip(tree[2*node_i+1].beliefs, tree[2*node_i+2].beliefs))) + list(
                tree[2*node_i+2].beliefs))
            node.state = node_states[2]
            node_i = floor((node_i-1)/2)
            node, depth = tree[node_i], depth-1
    del tree
    gc.collect()
    return vector(GF(2), vhat)
