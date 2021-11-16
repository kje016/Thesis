# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


# TODO: This is for BSC
class Node:
    def __init__(self, l_child, r_child, state):
        self.l_child = l_child
        self.r_child = r_child
        self.beliefs = []
        self.state = state

    def __str__(self):
        return f'beliefs:{str(self.beliefs)}, state:{self.state}, left_child: {str(self.l_child)}, right_child:{str(self.r_child)}'


def sign(rv):
    if rv <= 0:
        return 0
    return 1


def f(beliefs):
    result = []
    for x, y in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)]):
        result.append((sign(x)*sign(y))*min(abs(x),abs(y)))
    return vector(F, result)


def g(beliefs, beta):
    result = []
    for a1, a2, b in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)], beta):
        #breakpoint()
        result.append( (a2 + (1-2*b)*a1) )
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
    depth, done, node = 0, False, tree[0]
    while not done:
        if depth == log(N, 2):
            #breakpoint()
            if tree.index(node)-(2**log(N, 2)-1) in frozen_set:
                node.beliefs = vector(F, [0])
            else:
                if node.beliefs <= 0:
                    node.beliefs = vector(F, [0])
                else:
                    node.beliefs = vector(F, [1])
            if tree.index(node) == len(tree)-1:
                done = True
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
            #breakpoint()
            node.beliefs = vector(F, list(map(lambda x: mod(x[0]+ x[1],2), zip(tree[node.l_child].beliefs, tree[node.r_child].beliefs))) + list(tree[node.r_child].beliefs))
            node.state = node_states[2]
            node, depth = tree[floor((tree.index(node)-1)/2)], depth-1# step u


    message = []
    for i in list(set(range(2**log(N, 2)-1, len(tree))) - set(frozen_set)):
        message.extend(tree[i].beliefs)
    print(message)
    for x in reversed(frozen_set):
        message.pop(x)
    print(f"sc output:= {message}")
    return message

node_states = ['l', 'r', 'u']
F = RealField(10)
r = vector(F, [1, 0, 1, 0, 1, 1, 0, 1])
p = 0.2
llr1 = p/(1-p)
initLLRS = (2*r-vector([1]*len(r)))*(1-p)
breakpoint()
message = list(BPSK_decoder(vector(F, [-1.38, -1.38, -1.38, 1.38, -1.38, 1.38, -1.38, 1.38]), 8, [0, 1, 2, 4]))


