from sage.all import *
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Subchannel_Allocation
import CRC

# cd Desktop/Thesis/PySageMath/PC
run_vals = [[1/2, 0.14, 12, 'BSC'], [2/5, 0.12, 15, 'BSC'], [2/5, 0.1, 15, 'BSC'], [2/5, 0.08, 21, 'BSC'],
            ]
I_IL = 1
runs = 10
decoder = 'SCL'
for elem in run_vals:
    rate = elem[0]
    snr = elem[1]
    A = elem[2]
    channel = elem[3]

    pol = CRC.get_pol(A, I_IL)
    K = A + pol.degree()
    N0 = None
    false_negative = 'NaN'

    E = ceil(K / rate)
    n = min(ceil(log(E, 2)), 10)     # TODO: hvordan velges egt 'n'?
    N = 2 ** n
    QN0 = PC_Subchannel_Allocation.get_Q_N0(N)
    npc, n_wm_pc = PC_Subchannel_Allocation.get_n_pc_bits(K, E, I_IL)
    QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, E, npc, QN0, 0)
    GN = PC_Encoder.gen_G(n)
    QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)

    BLER, BER = 0, 0
    print(elem)
    for iteration in range(runs):
        a = random_vector(GF(2), A)
        #a = vector(GF(2), [1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0])
        #pol = x**8 + x**7 + x**4 + x**3 + x + x**0  #TODO: remove this line
        H = CRC.gen_CRC_mat(K, pol)
        Hperm = H.rows()
        c, G = CRC.CRC(a, A, pol)
        #ct = vector(GF(2), list(a) + (vector(GF(2), a) * H[-A:]).list())

        c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL=I_IL, c_seq=c)
        c_ap = vector(GF(2), c_ap)


        iPI = [PI.index(a) for a in PI if a >= A]
        coreperm = [a for a in PI if a < A]

        Gperm = Matrix(G[a] for a in [y for y in PI if y < A])
        Hperm = Matrix(H[a] for a in PI)

        end_perm = [c[-A:][a] for a in PI if a < A]
        ce = vector(GF(2), list(c[:-A]) + end_perm)
        ht = Matrix(GF(2), H.rows()[:-A] + Gperm.rows())



        breakpoint()
