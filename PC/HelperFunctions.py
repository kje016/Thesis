from sage.all import *


# return the reliability sequence as [ints] 1 indexed
def get_realiability_sequence():
    with open('Reliability_Sequence.txt') as f:
        reliability_sequence = f.read()
    reliability_sequence = reliability_sequence.replace(' ', '')
    reliability_sequence = reliability_sequence.split(",")
    reliability_sequence = [int(x) for x in reliability_sequence]
    return reliability_sequence

# @arg float rv: real value digit
def sign(rv):
    if rv < 0:
        return -1
    return 1

# @arg float rv: real value digit
# reverse of sign(rv)
def sign_rev(rv):
    if rv >= 0:
        return 0
    return 1

def BPSK(codeword):
    return [(-1)**x for x in codeword]

def f_BPSK(alpha1, alpha2):
    result = []
    for x,y in zip(alpha1,alpha2):
        result.append((sign(x)*sign(y))*min(abs(x), abs(y)))
    return vector(RR, result)

def g_BPSK(alpha1, alpha2, beta):
    result = []
    for a1,a2,b in zip(alpha1, alpha2, beta):
        result.append( (a2 + (1-2*b)*a1) )
    return vector(RR, result)




class Decoder:
    def __init__(self, cword, path_metric):
        self.cword = cword
        self.path_metric = path_metric

    def __repr__(self):
        pp = '('+ self.cword + ', ' + str(self.path_metric) + ')'
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
def update_decoders(is_frozen_node, belief, input_decoders, n_decoders):
    new_decoders = []
    if not is_frozen_node:
        #expand decoders list
        new_decoders = [Decoder(decoder.inf_bits, decoder.path_metric) for decoder in input_decoders]
        penalty = abs(1-belief)
        for decoder in new_decoders:
            decoder.cword = decoder.cword + "1"
            decoder.path_metric = decoder.path_metric + penalty

    penalty = abs(0 - belief)
    for decoder in input_decoders:
        decoder.inf_bits = decoder.inf_bits + "0"
        decoder.path_metric = decoder.path_metric + penalty
    input_decoders.extend(new_decoders)
    input_decoders = prune_decoders(input_decoders, n_decoders)
    return input_decoders