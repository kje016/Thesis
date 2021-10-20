n_decoders = 4
class Decoder:
    def __init__(self, cword, path_metric):
        self.cword = cword
        self.path_metric = path_metric

    def __repr__(self):
        pp = '('+ self.cword + ', ' + str(self.path_metric) + ')'
        return pp

"""brukt for å generere n antall decoders"""
#decoders = [Decoder("", 0) for i in range(n_decoders)]
"""Test for prune"""
dc0 = Decoder("", 0)
dc1 = Decoder("0000", 1)
dc2 = Decoder("0000", 1)
dc3 = Decoder("0000", 1)
dc4 = Decoder("0000", 1)

dc5 = Decoder("0001", 0)
dc6 = Decoder("0001", 0)
dc7 = Decoder("0001", 0)
dc8 = Decoder("0001", 0)

#decoders = [dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8]
decoders = [dc0]
#sort @decoders by 'path_metric' in the natural order



def print_decoders():
    for x in decoders:
        print(repr(x))

def prune_decoders(input_decoders):
    if len(input_decoders) <= n_decoders:
        return input_decoders
    input_decoders.sort(key=lambda x: x.path_metric)
    input_decoders = input_decoders[0:n_decoders]
    return input_decoders

def update_decoders(is_frozen_node, belief, input_decoders):
    new_decoders = []
    breakpoint()
    if not is_frozen_node:
        #expand decoders list
        new_decoders = [Decoder(decoder.inf_bits, decoder.path_metric) for decoder in input_decoders]
        penalty = abs(1-belief)
        for decoder in new_decoders:
            decoder.cword = decoder.cword + "1"
            decoder.path_metric = decoder.path_metric + penalty

    breakpoint()
    penalty = abs(0 - belief)
    for decoder in input_decoders:
        decoder.inf_bits = decoder.inf_bits + "0"
        decoder.path_metric = decoder.path_metric + penalty
    breakpoint()
    input_decoders.extend(new_decoders)
    return input_decoders



decoders = prune_decoders(decoders)
decoders = update_decoders(False, 0, decoders)
breakpoint()


