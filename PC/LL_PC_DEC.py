# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial


# TODO: This is for BSC
class Node:
    def __init__(self, parent, l_child, r_child, state):
        self.parent = parent
        self.l_child = l_child
        self.r_child = r_child
        self.beliefs = []
        self.state = state

    def __str__(self):
        return f'parent :{str(self.parent)}, beliefs:{str(self.beliefs)}, state:{self.state}, left_child: {str(self.l_child)}, right_child:{str(self.r_child)}'


def sign(rv):
    return -1 if rv < 0 else 1


def f(beliefs):
    result = []
    for x, y in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)]):
        result.append((sign(x)*sign(y))*min(x, y))
    return vector(RealField(5), result)


def g(beliefs, beta):
    breakpoint()
    breakpoint()
    result = []
    for x0, x1, x2 in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)], beta):
        result.append( InfinitePolynomial(RealField(5), (x1 + (1-2*x2)*x0)) )
    return vector(P, result)


def init_tree(N, r):
    tree, d, n = [], 0, 1
    while n < N:
        tree.extend([Node(d, i, i+1, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
        d += 1
        n *= 2
    tree.extend([Node(d, None, None, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
    tree[0].beliefs = r
    return tree


node_states = ['l', 'r', 'u']


def BPSK_decoder(d, N, F):
    tree = init_tree(N, d)
    depth, node = 0, tree[0]
    while tree[0].state != "u":
        if depth == log(N, 2):
            if tree.index(node)-(2**log(N, 2)-1) in F:
                node.beliefs = vector(GF(2), [0])
            #breakpoint()
            node.state = node_states[2]
            node, depth = tree[floor(abs((tree.index(node)-1))/2)], depth-1
            #breakpoint()
        elif node.state == "":
            tree[node.l_child].beliefs = f(node.beliefs)
            node.state = node_states[0]
            node, depth = tree[node.l_child], depth+1 # step l
            #breakpoint()
        elif node.state == node_states[0]:
            #breakpoint()
            tree[node.r_child].beliefs = g(node.beliefs, tree[node.l_child].beliefs)
            node.state = node_states[1]
            node, depth = tree[node.r_child], depth+1# step r
            #breakpoint()
        else:
            node.beliefs = vector(GF(2), list(tree[node.l_child].beliefs + tree[node.r_child].beliefs) + list(tree[node.r_child].beliefs))
            node.state = node_states[2]
            node, depth = tree[floor((tree.index(node)-1)/2)], depth-1# step u
            #breakpoint()

    breakpoint()
    return tree[0].beliefs


