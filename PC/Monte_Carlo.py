# cd Desktop/Thesis/PySageMath/PC
from sage.all import *
import csv

import numpy as np
import PC_Decoding
import PC_Encoding
import HelperFunctions as HF

import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt
from scipy.stats import norm

import sys

# sage Monte_Carlo.py 12 0 AWGN
R = [1/2, 1/3, 1/4]    # Rate of the code
A_min = 12
# runs = 5
lim = 200  # number or errors to be reached before moving to the next SNR value
SNR = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]) #,  5, 5.5, 6])
A = int(sys.argv[1])
I_IL, channel = int(sys.argv[2]), sys.argv[3].upper()
#plist = list(map(lambda rate: sqrt(1 / (2 * rate * 10 ** (SNR[10] / 10))), R))

BEgraph = np.zeros((3, len(SNR)))   # Vector containing P_b
np.random.seed(314150304)   # Random seed. Use your own
for i, snr in enumerate(SNR):
    #sigma = sqrt(1 / (2 * R * 10 ** (snr / 10)))  # TODO: SNR[] hard-coded
    BE = np.zeros(3)
    runs = 0
    print(f"snr := {snr}")
    while np.min(BE) < lim:
        runs += 1
        a = random_vector(GF(2), A)
        sigma = sqrt(1 / (2 * R[0] * 10 ** (snr / 10)))
        N0 = 2 * sigma** 2
        N, e, QNF, ms, MS, pol = PC_Encoding.PC_Encoding(a=list(a), A=A, R=R[0], I_IL=I_IL)
        r = list(HF.channel_noise(s=e, channel=channel, p=sigma))
        sclout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=ms, MS=MS,
                                         p_cross=sigma, channel=channel, pol=pol)
        BE[0] = BE[0] + (a+sclout[0]).hamming_weight()

        sigma = sqrt(1 / (2 * R[1] * 10 ** (snr / 10)))
        N0 = 2 * sigma ** 2
        N, e, QNF, ms, MS, pol = PC_Encoding.PC_Encoding(a=list(a), A=A, R=R[1], I_IL=I_IL)
        r = list(HF.channel_noise(s=e, channel=channel, p=sigma))
        sclout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=ms, MS=MS,
                                         p_cross=sigma, channel=channel, pol=pol)
        BE[1] = BE[1] + (a + sclout[0]).hamming_weight()

        sigma = sqrt(1 / (2 * R[2] * 10 ** (snr / 10)))
        N0 = 2 * sigma** 2
        N, e, QNF, ms, MS, pol = PC_Encoding.PC_Encoding(a=list(a), A=A, R=R[2], I_IL=I_IL)
        r = list(HF.channel_noise(s=e, channel=channel, p=sigma))
        sclout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=ms, MS=MS,
                                         p_cross=sigma, channel=channel, pol=pol)
        BE[2] = BE[2] + (a + sclout[0]).hamming_weight()
        if runs > 500:
            break
    for j in range(3):
        BEgraph[j, i] = BE[j]/(runs*A)

    plt.plot(SNR, np.transpose(BEgraph))
    plt.yscale('log')
    plt.grid(axis='y')
    plt.legend(['1/2', '1/3', '1/4'])
    plt.xlabel('SNR (dB)')
    plt.ylabel('Pb')
plt.savefig(f'20_chart.png')




