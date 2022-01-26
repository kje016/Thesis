# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import HelperFunctions as HF


class Node:
    def __init__(self, l_child, r_child, state):
        self.l_child = l_child
        self.r_child = r_child
        self.beliefs = []
        self.state = state

    def __str__(self):
        return f'beliefs:{str(self.beliefs)}, state:{self.state}, left_child: {str(self.l_child)}, right_child:{str(self.r_child)}'


# @arg float rv: real value digit
def sign(rv):
    return -1 if rv<0 else 1


def ft(beliefs):
    result = []
    for a1, a2 in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)]):
        result.append((sign(a1)*sign(a2))*min(abs(a1), abs(a2)))
    return vector(F, result)



def gt(beliefs, beta):
    result = []
    for a1, a2, b in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)], beta):
        result.append(a2 + a1*(1-2*b))
    return vector(F, result)




def decoder(d, N, frozen_set, p_cross):
    llr = log(p_cross / (1 - p_cross))
    tree = init_tree(N, d)
    """SCL initialization"""
    L = 8   # TODO: list size of L = 8. From "The Development and Operations of...
    list_decoders = [HF.Decoder("", 0)]
    depth, done, node = 0, False, tree[0]
    while not done:
        if depth == log(N, 2):
            node.state = node_states[2]
            is_frozen = tree.index(node)-(2**log(N, 2)-1) in frozen_set # alternatively var name,
            list_decoders = HF.update_decoders(is_frozen, node.beliefs[0], llr,  list_decoders, L)
            node.beliefs = vector(F, [HF.sign_rev(node.beliefs)])

            if tree.index(node) == len(tree)-1:
                done = True
            node, depth = tree[floor(abs((tree.index(node)-1))/2)], depth-1
        elif node.state == "":  # step L
            tree[node.l_child].beliefs = ft(node.beliefs)
            node.state = node_states[0]
            node, depth = tree[node.l_child], depth+1
        elif node.state == node_states[0]:  # step R
            tree[node.r_child].beliefs = gt(node.beliefs, tree[node.l_child].beliefs)
            node.state = node_states[1]
            node, depth = tree[node.r_child], depth+1
        else:   # step U
            node.beliefs = vector(F, list(map(lambda x: mod(x[0]+ x[1],2), zip(tree[node.l_child].beliefs, tree[node.r_child].beliefs))) + list(tree[node.r_child].beliefs))
            node.state = node_states[2]
            node, depth = tree[floor((tree.index(node)-1)/2)], depth-1

    for dec in list_decoders:
        bits = set(list(range(0, N))) - set(frozen_set)
        bits = [dec.inf_bits[a] for a in bits]
        dec.inf_bits = vector(GF(2), [str(a) for a in bits]) #''.join(a for a in bits)

    # TODO: not always the case that the most likely codeword was the sent symbol
    return list_decoders



node_states = ['l', 'r', 'u']
F = RealField(10)
#message = list(BPSK_decoder(vector(F, [1, 0, 1, 0, 0, 1, 0, 1]), 8, [0, 1, 2, 4], 0.2))
#print(f"llr_r:={message}")

