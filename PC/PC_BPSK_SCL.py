from sage.all import *
import numpy as np
import HelperFunctions as HF
import pdb

""" Polar code initialization"""
n_int = 3
N = 2**n_int
K = 4

reliability_sequence = HF.get_realiability_sequence()
reliability_sequence_N = []
for index, elem in enumerate(reliability_sequence):
    if elem <=N:
        reliability_sequence_N.append(elem-1)
    if len(reliability_sequence_N) >= N:
        break

frozen_set = reliability_sequence_N[0:(N-K)]

#TODO: codeword should be generated from PC_Encoder.py
d = [1, 0, 1, 0, 0, 1, 0, 1]

"""communication initialization"""
EbNodB = 4
Rate = K/N
EbNo = n(10**(EbNodB/10), digits=4)
sigma = sqrt(1/(2*Rate*EbNodB))
d_modulated = HF.BPSK(d)
r = d_modulated + sigma * np.random.normal(0, 1, N)
print("D_MODULATED", d_modulated)
print("NORMAL", r)

""""Tree initialization"""
node_state = [0] * (2 * N - 1)
beliefs = np.zeros((n_int + 1, N), dtype=float)
beliefs = matrix(RR, beliefs)

for i in range(N):  # setting belief for root node
    beliefs[0, i] = r[i]

"""tree traverse variables"""
node, depth, done, npos = 0, 0, 0, 0
ucap = np.zeros((n_int+1,N), dtype=int)
ucap = matrix(ucap)

"""SCL initialization"""
#TODO: list size of L =8. From "The Development and Operations of...
n_decoders = 4
list_decoders = [HF.Decoder("", 0) for i in range(n_decoders)]

#breakpoint()
while done == 0:
    npos = (2 ** depth - 1) + node  # position of node in node state vector

    # check if leaf
    if depth == n_int:
        #breakpoint()
        node_state[npos] = 3
        decision_metric = HF.sign_rev(beliefs[n_int, node])
        list_decoders = HF.update_decoders((node in frozen_set), decision_metric, list_decoders, n_decoders)
        ucap[n_int, node] = decision_metric
        #breakpoint()
        if node == (N-1):
            done = 1
        else:
            node = floor(node/2); depth = depth-1

    # is non-leaf
    else:
        temp = 2 ** (n_int - depth)
        belief_node = list(beliefs[depth][temp * node:temp * (node + 1)])
        alpha_l = belief_node[0:int(temp / 2)]
        alpha_r = belief_node[int(temp / 2): N]

        if node_state[npos] == 0:   # node_state == 0 => f(a,b)
            node, depth, temp = node*2, depth+1, temp/2
            #breakpoint()
            beliefs[depth, int(temp * node): int(temp * (node + 1))] = HF.f_BPSK(alpha_l, alpha_r)
            node_state[npos] = 1

        elif node_state[npos] == 1:   # node_state == 1 => g(a,b,beta)
            left_node, left_depth, left_temp = 2*node, depth+1, temp/2
            node, depth, temp = node*2+1, depth+1, temp/2

            beta = list(ucap[left_depth][int(left_temp * left_node): int(left_temp * (left_node + 1))])
            beliefs[depth, int(temp * node): int(temp * (node + 1))] = HF.g_BPSK(alpha_l, alpha_r, beta)
            node_state[npos] = 2;

        # node_state == 2  => [a+b b]
        else:
            left_node, right_node = 2*node, 2*node+1
            cdepth, ctemp = depth + 1, int(temp/2)
            ucapnl = ucap[cdepth][int(ctemp * left_node): ctemp * (left_node + 1)]
            ucapnr = ucap[cdepth][int(ctemp * right_node): ctemp * (right_node + 1)]

            #breakpoint()
            result = list(mod(a+b, 2) for a, b in zip(ucapnr, ucapnl))
            result.extend(ucapnr)

            ucap[depth, int(temp * node): int(temp * (node + 1))] = vector(ZZ, result)
            node, depth, node_state[npos] = floor(node / 2), depth - 1, 3


breakpoint()
message = [x for x in list_decoders[0].cword]
for x in reversed(frozen_set):
    message.pop(x)


print("UCAP")
print(ucap)

print("BELIEFS")
print(beliefs)

print("decoder_input: ", r)
print("MESSAGE: ", message)

