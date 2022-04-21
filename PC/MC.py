# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import sys
import csv
import datetime
import gc

import PC_Decoding
import PC_Encoder
import PC_Rate_Matching
import PC_Subchannel_Allocation
import HelperFunctions as HF
import PC_CRC
import test_CRC

# sage MC.py 12 0 BSC SC
R = [1/2] # [1/2, 2/5, 1/3, 1/4,  1/5]   # Rate of the code
A_min = 12
runs = 100
SNR = [0.1] #[0.3, 0.2, 0.1]   # really p_cross
# SNR = [0.5, 1, 2, 3, 4]
A = int(sys.argv[1])
I_IL = int(sys.argv[2])
channel, decoder = sys.argv[3].upper(), sys.argv[4].upper()
pol = PC_CRC.get_pol(A, I_IL)
K = A + pol.degree()
n_min, n_max = 5, 10 - I_IL  # n_max = 10 for uplink, 9 for downlink.

for rate in R:
    E = ceil(K / rate)
    n = min(ceil(log(E, 2)), n_max)     # TODO: hvordan velges egt 'n'?
    N = 2 ** n
    QN0 = PC_Subchannel_Allocation.get_Q_N0(N)
    npc, n_wm_pc = PC_Subchannel_Allocation.get_n_pc_bits(K, E, I_IL)
    QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, E, npc)
    GN = PC_Encoder.gen_G(n)
    QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)

    for i, snr in enumerate(SNR):
        # sigma = sqrt(1 / (2 * rate * 10 ** (snr / 10)))
        # N0 = 2 * sigma ** 2
        BLER, BER = 0, 0
        false_negative = 0  # counts of how many times the codeword is in scout, but is not the most likely output
        print(f"p_cross = {snr}")
        for iteration in range(runs):
            if iteration % 100 == 0:
                print(iteration)
            a = random_vector(GF(2), A)
            c, C = test_CRC.CRC(a, A, pol)
            u = PC_Subchannel_Allocation.calc_u(N, QNI, c, QNPC)
            d = vector(GF(2), u) * GN
            MS = PC_Rate_Matching.get_rm_set(N-E, matching_scheme, QN0)
            e = PC_Rate_Matching.circular_buffer(d, MS, matching_scheme)
            r = list(HF.channel_noise(s=e, channel=channel, p=snr)) # returns the modulated codeword with added noise
            scout = PC_Decoding.PC_Decoding(r=r, N=N, N0=None, QNF=QNF, ms=matching_scheme, MS=MS,
                                             p_cross=snr, channel=channel + '_' + decoder, I_IL=I_IL)
            BER = BER + (a + scout[:A]).hamming_weight()
            BLER = BLER + sign(test_CRC.CRC_check(scout, K, pol).hamming_weight())

            if decoder == 'SCL':
                if 1 * sign((a+scout[0].inf_bits).hamming_weight()):
                    for scout_i in range(1, len(scout)):
                        if scout[scout_i].inf_bits == a:
                            false_negative = false_negative + 1
                            break


        file_getter = channel + '_' + decoder
        with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\{file_getter}.csv', mode='a',
                  newline='') as file:
            result_writer = csv.writer(file)  # , delimeter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            result_writer.writerow(
                [A, rate, A, N, runs, BER / (runs * A), BLER / runs, false_negative, snr, datetime.datetime.now()])
            gc.collect()
