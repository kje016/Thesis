from sage.all import *
from numpy.random import default_rng
from numpy.random import uniform


def standard_to_list(input_text):
    temp = input_text.read()
    temp = temp.replace('\n', ' ')
    temp = temp.split(' ')
    rel = {}
    for i in range(0, len(temp), 2):
        rel.update({int(temp[i]): int(temp[i+1])})
    rel = [rel.get(i) for i in range(len(rel))]
    return rel


# return the reliability sequence as [ints] 1 indexed
def get_realiability_sequence():
    with open('Reliability_Sequence.txt') as f:
        reliability_sequence = f.read()
    reliability_sequence = reliability_sequence.split(",")
    reliability_sequence = [int(x) for x in reliability_sequence]
    return reliability_sequence


# 2 is representing the erasure symbol
# returns the modulation of the codeword with added noise
def channel_noise(s, channel, p):
    F = RealField(7)
    if channel == 'BSC':
        noise = vector(F, [1 if x <= p else 0 for x in list(uniform(0, 1, size=len(s)))])
        r = vector(F, list(map(lambda y: (2 * y) - 1, (vector(F, s)+noise) % 2)))

    elif channel == 'AWGN':
        noise = vector(F, list(default_rng().normal(0, p, len(s))))
        r = 2*vector(F, s) - vector(F, [1]*len(s)) + noise

    else: # channel == 'BEC'
        #noise = list(map(lambda lis, i: lis[i] = 1, [0]*len(s)) list(random.sample(range(0, len(s)), floor(len(s)*p))))
        noise = vector(F, [1 if x <= p else 0 for x in list(uniform(0, 1, size=len(s)))])
        s_mod = vector(F, list(map(lambda y: (2 * y) - 1, vector(F, s))))
        r = vector(F, [2 if noise[i] == 1 else s_mod[i] for i, e in enumerate(s_mod)])
    return r


# @arg float rv: real value digit
def sign(rv):
    return -1 if rv<=0 else 1


# @arg float rv: real value digit
# reverse of sign(rv)
def sign_rev(rv):
    return 0 if rv >= 0 else 1


def ft(beliefs, F):
    result = []
    for a1, a2 in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)]):
        result.append((sign(a1)*sign(a2))*min(abs(a1), abs(a2)))
    return vector(F, result)


def gt(beliefs, beta, F):
    result = []
    for a1, a2, b in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)], beta):
        result.append(a2 + a1*(1-2*b))
    return vector(F, result)


def bec_xor(a1, a2):
    if 2 in [a1, a2]:
        return 2
    else:
        return ((sign_rev(a1)^sign_rev(a2))*2-1)*-oo


def bec_uhat(belief, frozen):
    if frozen:
        return vector(RealField(7), [oo])
    else:
        return vector(RealField(7), [belief[0]*oo])


def uhat(belief, frozen, F):
    if frozen:
        return vector(F, [0])
    else:
        return vector(F, [sign_rev(belief[0])])


def f_bsc(beliefs):
    F = RealField(7)
    result = []
    for a1, a2 in zip(beliefs[0:len(beliefs) // 2], beliefs[len(beliefs) // 2: len(beliefs)]):
        result.append((sign(a1)*sign(a2))*min(abs(a1), abs(a2)))
    return vector(F, result)


def f_bec(beliefs):
    F = RealField(7)
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
    F = RealField(7)
    result = []
    for a1, a2, b in zip(beliefs[0:len(beliefs)//2], beliefs[len(beliefs)//2: len(beliefs)], beta):
        result.append(lor(a2, bec_xor(a1, b)))
    return vector(F, result)


class Decoder:
    def __init__(self, inf_bits, path_metric):
        self.inf_bits = inf_bits
        self.path_metric = path_metric

    def __repr__(self):
        pp = '('+ str(self.inf_bits) + ', ' + str(self.path_metric) + ')'
        return pp


class Node:
    def __init__(self, state, child):
        self.beliefs = []
        self.state = state
        self.child = child

    def __str__(self):
        return f'beliefs:{str(self.beliefs)}, state:{self.state}, state:{self.child}'


def init_tree(N, r):
    tree, d, n = [], 0, 1
    while d <= log(N, 2):
        tree.extend([Node('', '') for i in range(n)])
        d, n = d+1, n*2
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
def update_decoders(is_frozen_node, belief, input_decoders, n_decoders, C_perm, crcbit):
    if is_frozen_node:
        return input_decoders

    if len(input_decoders[0].inf_bits) in crcbit:
        new_decoders = []
        iPI = crcbit[:crcbit.index(len(input_decoders[0].inf_bits))]
        #breakpoint()
        for decoder in input_decoders:
            cword = list(vector(GF(2), decoder.inf_bits))
            for a in iPI[::-1]:
                cword.pop(a)
            #cword = vector(GF(2), list(map(lambda x: cword.pop(x), iPI[::-1]))
            check = (vector(GF(2), cword) * Matrix(GF(2), C_perm[:len(cword)]))[len(iPI)]
            if check == sign_rev(belief):
                new_decoders.append(Decoder(decoder.inf_bits + str(check), decoder.path_metric))
    else:
        new_decoders = [Decoder(decoder.inf_bits+"1", decoder.path_metric + (1-sign_rev(belief))*abs(belief)) for decoder in input_decoders]
        new_decoders.extend([Decoder(decoder.inf_bits + "0",  decoder.path_metric + sign_rev(belief)*abs(belief)) for decoder in input_decoders])
    new_decoders = prune_decoders(new_decoders, n_decoders)
    return new_decoders


def bec_update_decoders(is_frozen_node, belief, llr,  input_decoders, L, C_perm, crcbit):
    if is_frozen_node:
        return input_decoders

    if abs(belief) != oo:
        new_decoders = [Decoder(decoder.inf_bits + "1", decoder.path_metric + abs(llr)) for
                        decoder in input_decoders]
        new_decoders.extend([Decoder(decoder.inf_bits + "0", decoder.path_metric + abs(llr)) for decoder in
             input_decoders])

    elif len(input_decoders[0].inf_bits) in crcbit:
        new_decoders = []
        for decoder in input_decoders:
            iPI = crcbit[:crcbit.index(len(decoder.inf_bits))][::-1]
            cword = list(vector(GF(2), decoder.inf_bits))
            list(map(lambda x: cword.pop(x), iPI))
            check = (vector(GF(2), cword) * Matrix(GF(2), C_perm[:len(cword)]))[len(iPI)]
            if float(check) == uhat([belief], False, RealField(10))[0]:
                new_decoders.append(Decoder(decoder.inf_bits + str(check), decoder.path_metric))
    else:
        new_decoders = [Decoder(decoder.inf_bits + str(sign_rev(belief)), decoder.path_metric)for decoder in input_decoders]
    return prune_decoders(new_decoders, L)
