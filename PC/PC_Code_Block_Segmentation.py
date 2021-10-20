from sage.all import *


def main_block_segmentation(channel, a, A, G):
    I_seq = channel in ["PUCCH", "PUSCH"] and (A >= 1013 or A >= 360 or G >= 1088)
    C = 2 if I_seq else 1
    a_ap = a[:]
    if C == 2:
        if A%2 != 0:
            A_ap = (A+1)/2
            a_ap.inesert(0, 0)
        else:
            A_ap = A/2
        a_ap = [a_ap[i:A_ap] for i in range(0, A, A / 2)]
        return a_ap, C
    return [a_ap], C






