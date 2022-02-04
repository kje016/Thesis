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
runs = 5
lim = 50  # number or errors to be reached before moving to the next SNR value
SNR = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6])
A = list(range(12, 40))
I_IL, channel = int(sys.argv[2]), sys.argv[3].upper()
plist = list(map(lambda rate: sqrt(1 / (2 * rate * 10 ** (SNR[10] / 10))), R))

BEgraph = np.zeros((3, len(A)))   # Vector containing P_b
np.random.seed(314150304)   # Random seed. Use your own
for i, Ai in enumerate(A):
    BE = np.zeros(len(R))
    while np.min(BE) < lim:
        a = random_vector(GF(2), Ai)
        for ri, rate in enumerate(R):
            sigma = sqrt(1 / (2 * rate * 10 ** (SNR[10] / 10)))  # TODO: SNR[] hard-coded
            N0 = 2 * sigma ** 2
            N, e, QNF, ms, MS, pol = PC_Encoding.PC_Encoding(a=list(a), A=Ai, R=rate, I_IL=I_IL)
            r = list(HF.channel_noise(s=e, channel=channel, p=sigma))
            sclout = PC_Decoding.PC_Decoding(r=r, N=N, N0=N0, QNF=QNF, ms=ms, MS=MS,
                                             p_cross=sigma, channel=channel, pol=pol)
            if not sclout[1]:
                BE[ri] += 1
            #BE[ri] = BE[ri]+(a+sclout[0]).hamming_weight()
    for j in range(len(R)):
        BEgraph[j, i] = BE[j]/(runs)   # Bit Error Probability
    #plt.plot(A, plist*len(A))
    plt.plot(A, np.transpose(BEgraph))
    plt.yscale('log')
    plt.grid(axis='y')
    plt.legend(['1/2', '1/3', '1/4'])
    plt.xlabel('A')
    plt.ylabel('Pb')
plt.savefig(f'A_chart.png')




