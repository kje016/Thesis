from sage.all import *
from numpy.random import default_rng
from numpy.random import uniform

import test_CRC


def standard_to_list(input_text):
    with open(input_text) as f:
        temp = f.read()
    temp = temp.replace('\n', ' ')
    temp = temp.split(' ')
    rel = {}
    for i in range(0, len(temp), 2):
        rel.update({int(temp[i]): int(temp[i+1])})
    return [rel.get(i) for i in range(len(rel))]


# return the reliability sequence as [ints] 1 indexed
def get_realiability_sequence():
    with open('Reliability_Sequence.txt') as f:
        reliability_sequence = f.read()
    reliability_sequence = reliability_sequence.split(",")
    reliability_sequence = [int(x) for x in reliability_sequence]
    return reliability_sequence


def get_inf_pos(K):
    output = []
    R = get_realiability_sequence()
    for i in R:
        if i <= K:
            output.append(i)
        if len(output) >= K:
            return output
    return output


# 2 is representing the erasure symbol
# returns the modulation of the codeword with added noise
def channel_noise(s, channel, p):
    F = RealField(10)
    if channel == 'BSC':
        noise = vector(F, [1 if x <= p else 0 for x in list(uniform(0, 1, size=len(s)))])
        r = vector(F, list(map(lambda y: (2 * y) - 1, (vector(F, s)+noise) % 2)))

    elif channel == 'AWGN':
        noise = vector(F, list(default_rng().normal(0, p, len(s))))
        r = 2*vector(F, s) - vector(F, [1]*len(s)) + noise

    else: # channel == 'BEC'
        s_mod = vector(F, list(map(lambda y: (2 * y) - 1, vector(F, s))))
        r = vector(F, [2 if list(uniform(0, 1, size=len(s)))[i] <= p else s_mod[i] for i, e in enumerate(s_mod)])
    return r


# @arg float rv: real value digit
def sign(rv):
    return -1 if rv<=0 else 1


# @arg float rv: real value digit
# reverse of sign(rv)
def sign_rev(rv):
    return 0 if rv >= 0 else 1


def BPSK(codeword):
    return [(-1)**x for x in codeword]


def ft(beliefs):
    F = RealField(10)
    result = []
    for a1, a2 in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)]):
        result.append((sign(a1)*sign(a2))*min(abs(a1), abs(a2)))
    return vector(F, result)


def gt(beliefs, beta):
    F = RealField(10)
    result = []
    for a1, a2, b in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)], beta):
        result.append(a2 + a1*(1-2*b))
    return vector(F, result)


def bec_xor(a1, a2):
    if 2 in [a1, a2]:
        return 2
    else:
        return ((sign_rev(a1)^sign_rev(a2))*2-1)*-oo


def uhat(belief, frozen):
    if frozen:
        return vector(RealField(10), [oo])
    else:
        return vector(RealField(10), [belief[0]])

    # return vector(RealField(10), list(map(lambda a: a*-oo * (1-frozen), belief)))


def f_bec(beliefs):
    F = RealField(10)
    result = []
    for x, y in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)]):
        result.append(bec_xor(x, y))
    return vector(F, result)


def lor(a1, a2):
    if a1 == 2:
        return a2
    elif a2 == 2:
        return a1
    else:
        return ((sign_rev(a1) or sign_rev(a2))*2-1)*-oo


def g_bec(beliefs, beta):
    F = RealField(10)
    result = []
    for a1, a2, b in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)], beta):
        result.append(lor(a2, bec_xor(a1, b)))
    return vector(F, result)


def g_BPSK(alpha1, alpha2, beta):
    F= RealField(10)
    result = []
    for a1, a2, b in zip(alpha1, alpha2, beta):
        result.append( (a2 + (1-2*b)*a1) )
    return vector(F, result)


class Decoder:
    def __init__(self, inf_bits, path_metric):
        self.inf_bits = inf_bits
        self.path_metric = path_metric

    def __repr__(self):
        pp = '('+ str(self.inf_bits) + ', ' + str(self.path_metric) + ')'
        return pp


class Node:
    def __init__(self, l_child, r_child, state):
        self.l_child = l_child
        self.r_child = r_child
        self.beliefs = []
        self.state = state

    def __str__(self):
        return f'beliefs:{str(self.beliefs)}, state:{self.state}, left_child: {str(self.l_child)}, right_child:{str(self.r_child)}'


def init_tree(N, r):
    tree, d, n = [], 0, 1
    while n < N:
        tree.extend([Node(i, i+1, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
        d, n = d+1, n*2
    tree.extend([Node(None, None, '') for i in range(2*(n-1)+1, 2*(n-1)+1+(n*2), 2)])
    tree[0].beliefs = r
    return tree


# @arg input_decoders list of current decoders and their respective path_metric
# @arg decoder_size the amount of decoders to keep after pruning
# @return returns a pruned list of decoders
def prune_decoders(input_decoders, decoder_size):
    if len(input_decoders) <= decoder_size:
        return input_decoders
    input_decoders.sort(key=lambda x: x.path_metric)
    input_decoders = input_decoders[0:decoder_size]
    return input_decoders


# @arg is_frozen_node Boolean of the current node is a frozen node
# @belief the belief for the current node
# @input_decoders list of current decoders
def update_decoders(is_frozen_node, belief, llr,  input_decoders, n_decoders, C_perm, crcbit, I_IL):
    if is_frozen_node:
        new_decoders = [Decoder(decoder.inf_bits, decoder.path_metric) for decoder in input_decoders]
    else:
        if len(input_decoders[0].inf_bits) in crcbit:
            new_decoders = []
            for decoder in input_decoders:
                iPI = crcbit[:crcbit.index(len(decoder.inf_bits))][::-1]
                cword = list(vector(GF(2), decoder.inf_bits))
                list(map(lambda x: cword.pop(x), iPI))
                check = (vector(GF(2), cword) * Matrix(GF(2), C_perm[:len(cword)]))[len(iPI)]
                if check == sign_rev(belief):
                    new_decoders.append(Decoder(decoder.inf_bits + str(check), decoder.path_metric))
        else:
            new_decoders = [Decoder(decoder.inf_bits+"1", decoder.path_metric + (1-sign_rev(belief))*abs(llr)) for decoder in input_decoders]
            new_decoders.extend([Decoder(decoder.inf_bits + "0",  decoder.path_metric + (sign_rev(belief))*abs(llr)) for decoder in input_decoders])
    new_decoders = prune_decoders(new_decoders, n_decoders)
    return new_decoders
