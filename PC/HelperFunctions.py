from sage.all import *
from numpy.random import default_rng
from numpy.random import uniform

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

# @arg float rv: real value digit
def sign(rv):
    return -1 if rv < 0 else 1


# @arg float rv: real value digit
# reverse of sign(rv)
def sign_rev(rv):
    return 0 if rv >= 0 else 1


def BPSK(codeword):
    return [(-1)**x for x in codeword]


def f_BPSK(alpha1, alpha2):
    result = []
    for x, y in zip(alpha1, alpha2):
        result.append((sign(x)*sign(y))*min(abs(x), abs(y)))
    return vector(RR, result)


def channel_noise(s, channel, p):
    if channel == 'BSC':
        noise = [1 if abs(x) <= p else 0 for x in list(default_rng().normal(0, 1, len(s)))]
        #noise = [1 if x >= p else 0 for x in list(default_rng().uniform(0, 1, len(s)))]
        r = vector(GF(2), s) + vector(GF(2), noise)
    else:
        print("channel_noise () not fully implemented")
        r = s
    return r


def g_BPSK(alpha1, alpha2, beta):
    result = []
    for a1,a2,b in zip(alpha1, alpha2, beta):
        result.append( (a2 + (1-2*b)*a1) )
    return vector(RR, result)


class Decoder:
    def __init__(self, inf_bits, path_metric):
        self.inf_bits = inf_bits
        self.path_metric = path_metric

    def __repr__(self):
        pp = '('+ str(self.inf_bits) + ', ' + str(self.path_metric) + ')'
        return pp


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
def update_decoders(is_frozen_node, belief, llr,  input_decoders, n_decoders):
    new_decoders = []
    if not is_frozen_node:
        penalty = sign_rev(belief)
        #expand decoders list
        new_decoders = [Decoder(decoder.inf_bits+"1", decoder.path_metric + (1-sign_rev(belief))*abs(llr)) for decoder in input_decoders]
    for decoder in input_decoders:
        new_decoders.append(Decoder(decoder.inf_bits + "0", decoder.path_metric + (sign_rev(belief))*abs(llr)))
    new_decoders = prune_decoders(new_decoders, n_decoders)
    return new_decoders
