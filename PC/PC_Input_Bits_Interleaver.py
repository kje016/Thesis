# cd Desktop/Thesis/PySageMath/PC
from sage.all import *

K_max_il = 164
#pattern = "0 2 4 7 9 14 19 20 24 25 26 28 31 34 42 45 49 50 51 53 54 56 58 59 61 62 65 66 67 69 70 71 72 76 77 81 82 83 87 88 89 91 93 95 98 101 104 106 108 110 111 113 115 118 119 120 122 123 126 127 129 132 134 138 139 140 1 3 5 8 10 15 21 27 29 32 35 43 46 52 55 57 60 63 68 73 78 84 90 92 94 96 99 102 105 107 109 112 114 116 121 124 128 130 133 135 141 6 11 16 22 30 33 36 44 47 64 74 79 85 97 100 103 117 125 131 136 142 12 17 23 37 48 75 80 86 137 153 13 18 38 144 39 145 40 146 41 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163"
pattern = "0 2 4 7 9 14 19 20 24 25 26 28 31 34 42 45 49 50 51 53 54 56 58 59 61 62 65 66 67 69 70 71 72 76 77 81 82 83 87 88 89 91 93 95 98 101 104 106 108 110 111 113 115 118 119 120 122 123 126 127 129 132 134 138 139 140 1 3 5 8 10 15 21 27 29 32 35 43 46 52 55 57 60 63 68 73 78 84 90 92 94 96 99 102 105 107 109 112 114 116 121 124 128 130 133 135 141 6 11 16 22 30 33 36 44 47 64 74 79 85 97 100 103 117 125 131 136 142 12 17 23 37 48 75 80 86 137 143 13 18 38 144 39 145 40 146 41 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163"
interleaver_pattern = [int(x) for x in list(pattern.split(' '))]


def get_pi(K):
    PI = []
    for m in range(K_max_il):
        if interleaver_pattern[m] >= K_max_il-K:
            PI.append(interleaver_pattern[m]-K_max_il-K)
    return PI


def interleaver(flag_param, c_seq):
    if not flag_param:
        return c_seq

    PI = get_pi(len(c_seq))
    output = [c_seq[value] for value in PI]
    return output


def main_bit_interleaver(channel, c):
    flag = channel in ["PBCCH", "PDCCH"]    # interleaver is used in PBCH, PDCCH, bypassed for PUCCH & PUSCH
    c_ap = interleaver(flag, c)
    return flag, c_ap
