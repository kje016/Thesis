# cd Desktop/Thesis/PySageMath/PC
from sage.all import *

node_states = ['l', 'r', 'u']
F = RealField(10)


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
    return -1 if rv < 0 else 1


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


def init_tree(N, r):
    tree, d, n = [], 0, 1
    while n < N:
        tree.extend([Node(i, i+1, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
        d, n = d+1, n*2
    tree.extend([Node(None, None, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
    tree[0].beliefs = r
    return tree


def SC_decoder(d, N, frozen_set, p_cross):
    tree = init_tree(N, d)
    depth, done, node = 0, False, tree[0]
    while not done:
        if depth == log(N, 2):
            if tree.index(node)-(2**log(N, 2)-1) in frozen_set:
                node.beliefs = vector(F, [0])
            else:
                if node.beliefs >= 0:
                    node.beliefs = vector(F, [0])
                else:
                    node.beliefs = vector(F, [1])
            if tree.index(node) == len(tree)-1:
                done = True
            node.state = node_states[2]
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

    message = []
    for i in range(len(tree)-1, len(tree)-N-1, -1):
        message.extend(tree[i].beliefs)
    print(message)
    for x in reversed(frozen_set):
        message.pop(x)
    return vector(GF(2), message)
