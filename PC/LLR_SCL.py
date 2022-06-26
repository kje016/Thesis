# cd Desktop/Thesis/PySageMath/PC
import gc

from sage.all import *

import CRC
import HelperFunctions as HF
#import MC

node_states = ['l', 'r', 'u']
child = ['left', 'right']
F = RealField(10)


def decoder(d, N, frozen_set, I_IL, PI, H):
    pis = []
    if I_IL:
        pis = [PI.index(a) for a in PI if a >= len(PI)-24]
        C_perm = CRC.gen_CRC_mat(8, CRC.crc24)
        C_perm = [C_perm[a] for a in [a for a in PI if a < (H.ncols()-H.nrows())]]
        CRCpos = [PI.index(a) for a in PI if a >= 8]

    tree = HF.init_tree(N, d)
    F = RealField(10)
    """SCL initialization"""
    L = 8
    list_decoders = [HF.Decoder("", 0)]
    depth, done, node = 0, False, tree[0]
    node_i = tree.index(node)
    while not done:
        if depth == log(N, 2):
            is_frozen = node_i-(N-1) in frozen_set
            #list_decoders = HF.update_decoders(is_frozen, node.beliefs[0], list_decoders, L, H, PI, pis)
            list_decoders = HF.G_update(is_frozen, node.beliefs[0], list_decoders, L, C_perm, CRCpos)
            if not list_decoders:   # no surviving paths
                return [HF.Decoder('', +oo)]
            node.beliefs = vector(F, [HF.sign_rev(node.beliefs)*(not is_frozen)])

            if tree.index(node) == len(tree)-1:
                done = True
            node.state = node_states[2]
            node_i = floor(abs((node_i - 1)) / 2)
            node, depth = tree[node_i], depth - 1
        elif node.state == "":  # step L
            tree[2 * node_i + 1].beliefs = HF.ft(node.beliefs, F)
            node.state = node_states[0]
            node_i = 2 * node_i + 1
            node, depth = tree[node_i], depth + 1
            node.child = child[0]
        elif node.state == node_states[0]:  # step R
            tree[2 * node_i + 2].beliefs = HF.gt(node.beliefs, tree[2 * node_i + 1].beliefs, F)
            node.state = node_states[1]
            node_i = 2 * node_i + 2
            node, depth = tree[node_i], depth + 1
            node.child = child[1]
        else:   # step U
            node.beliefs = vector(F, list( map(lambda x: mod(x[0] + x[1], 2),
                zip(tree[2 * node_i + 1].beliefs, tree[2 * node_i + 2].beliefs))) + list(tree[2 * node_i + 2].beliefs))
            node.state = node_states[2]
            node_i = floor((node_i - 1) / 2)
            node, depth = tree[node_i], depth-1
    if PI:
        for dec in list_decoders:
            deinterleave = [0]*len(dec.inf_bits)
            for i, elem in enumerate(PI):
                deinterleave[elem] = dec.inf_bits[i]
            dec.inf_bits = vector(GF(2), deinterleave[:H.ncols()])
    else:
        for dec in list_decoders:
            dec.inf_bits = vector(GF(2), dec.inf_bits)
            if H*dec.inf_bits == 0:
                return dec.inf_bits
        return list_decoders[0].inf_bits
    del tree
    gc.collect()
    list_decoders.sort(key=lambda dec: dec.path_metric)
    return list_decoders
