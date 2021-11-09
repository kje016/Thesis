# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


# TODO: This is for BEC
class Node:
    def __init__(self, l_child, r_child, state):
        self.l_child = l_child
        self.r_child = r_child
        self.beliefs = []
        self.state = state

    def __str__(self):
        return f'beliefs:{str(self.beliefs)}, state:{self.state}, left_child: {str(self.l_child)}, right_child:{str(self.r_child)}'


def xor(a1, a2):
    if error_symbol in [a1, a2]:
        return error_symbol
    else:
        return Integer(a1).__xor__(Integer(a2))


def f(beliefs):
    result = []
    for x, y in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)]):
        result.append(xor(x, y))
    return vector(F, result)


def g(beliefs, beta):
    result = []
    for x0, x1, x2 in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)], beta):
        result.append( min(x1,xor(x0,x2)) )
    return vector(F, result)


def init_tree(N, r):
    tree, d, n = [], 0, 1
    while n < N:
        tree.extend([Node(i, i+1, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
        d, n = d+1, n*2
    tree.extend([Node(None, None, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
    tree[0].beliefs = r
    return tree


def BPSK_decoder(d, N, frozen_set):
    tree = init_tree(N, d)
    depth, node = 0, tree[0]
    while tree[0].state != "u":
        if depth == log(N, 2):
            if tree.index(node)-(2**log(N, 2)-1) in frozen_set:
                node.beliefs = vector(F, [F(0)])
            node.state = node_states[2]
            node, depth = tree[floor(abs((tree.index(node)-1))/2)], depth-1
        elif node.state == "":
            tree[node.l_child].beliefs = f(node.beliefs)
            node.state = node_states[0]
            node, depth = tree[node.l_child], depth+1 # step l
        elif node.state == node_states[0]:
            tree[node.r_child].beliefs = g(node.beliefs, tree[node.l_child].beliefs)
            node.state = node_states[1]
            node, depth = tree[node.r_child], depth+1# step r
        else:
            node.beliefs = vector(F, list(map(lambda x: xor(x[0], x[1]), zip(tree[node.l_child].beliefs, tree[node.r_child].beliefs))) + list(tree[node.r_child].beliefs))
            node.state = node_states[2]
            node, depth = tree[floor((tree.index(node)-1)/2)], depth-1# step u

    message = list(tree[0].beliefs)
    for x in reversed(frozen_set):
        message.pop(x)
    print(f"sc output:= {message}")
    return message

node_states = ['l', 'r', 'u']
error_symbol = 2
F = GF(3)
message = list(BPSK_decoder(vector(F, [2, 0, 2, 0, 2, 1, 0, 1]), 8, [0, 1, 2, 4]))


