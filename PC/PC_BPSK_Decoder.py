from sage.all import *
import numpy as np
import HelperFunctions as HF
import pdb

n_int = 3
N = 2**n_int
K = 4

EbNodB = 4
Rate = K/N
EbNo = n(10**(EbNodB/10), digits=4)
sigma = sqrt(1/(2*Rate*EbNodB))

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



def BPSK(codeword):
    return [(-1)**x for x in codeword]

# @arg float rv: real value digit
def sign(rv):
    if rv < 0:
        return -1
    return 1

def f(alpha1, alpha2):
    result = []
    for x,y in zip(alpha1,alpha2):
        result.append((sign(x)*sign(y))*min(abs(x), abs(y)))
    return vector(RR, result)



def g(alpha1, alpha2, beta):
    result = []
    for a1,a2,b in zip(alpha1, alpha2, beta):
        result.append( (a2 + (1-2*b)*a1) )
    return vector(RR, result)


d_modulated = BPSK(d)

#d_modulated = [ 1.1531403,  0.42099311,  1.04529288, -0.81021667,  0.96074431,  0.37959046, -0.86606408, -0.54630906]

r = d_modulated + sigma * np.random.normal(0, 1, N)
print("D_MODULATED", d_modulated)
print("NORMAL", r)


beliefs = np.zeros((n_int + 1, N), dtype=float)
beliefs = matrix(RR, beliefs)
node_state = [0] * (2 * N - 1)

# setting belief for root node
for i in range(N):
    beliefs[0, i] = r[i]

node = 0;
depth = 0;  # starting at root
done = 0

ucap = np.zeros((n_int+1,N), dtype=float)
ucap = matrix(ucap)

npos = 0
depth = 0

#breakpoint()
while done == 0:
    npos = (2 ** depth - 1) + node # position of node in node state vector

    # check if leaf
    if depth == n_int:
        node_state[npos] = 3
        if node in frozen_set:
            ucap[n_int, node] = 0
        else:
            #breakpoint()
            if beliefs[n_int, node] >= 0:
                ucap[n_int, node] = 0
            else:
                ucap[n_int, node] = 1
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


        # node_state == 0 => f(a,b)
        if node_state[npos] == 0:
            node = node * 2;
            depth = depth + 1
            temp = temp / 2
            #breakpoint()
            beliefs[depth, int(temp * node): int(temp * (node + 1))] = f(alpha_l, alpha_r)
            node_state[npos] = 1

        # node_state == 1 => g(a,b,beta)
        elif node_state[npos] == 1:
            left_node = 2 * node;
            left_depth = depth + 1
            left_temp = temp / 2

            node = node * 2 + 1;
            depth = depth + 1;
            temp = temp / 2

            beta = list(ucap[left_depth][int(left_temp * left_node): int(left_temp * (left_node + 1))])
            #breakpoint()
            beliefs[depth, int(temp * node): int(temp * (node + 1))] = g(alpha_l, alpha_r, beta)
            node_state[npos] = 2;

        # node_state == 2  => [a+b b]
        else:
            left_node = 2 * node;
            right_node = 2 * node + 1;
            cdepth = depth + 1;
            ctemp = int(temp / 2)
            ucapnl = ucap[cdepth][int(ctemp * left_node): ctemp * (left_node + 1)]
            ucapnr = ucap[cdepth][int(ctemp * right_node): ctemp * (right_node + 1)]

            #breakpoint()
            result = list( mod(a+b,2) for a, b in zip(ucapnr, ucapnl))
            result.extend(ucapnr)

            ucap[depth, int(temp * node): int(temp * (node + 1))] = vector(ZZ, result)
            node = floor(node / 2);
            depth = depth - 1
            node_state[npos] = 3


message = [int(x) for x in ucap[n_int]]
for x in reversed(frozen_set):
    message.pop(x)
#breakpoint()

print("UCAP")
print(ucap)

print("BELIEFS")
print(beliefs)

print("decoder_input: ", r)
print("MESSAGE: ", message)
