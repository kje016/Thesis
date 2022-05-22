from sage.all import *
import PC_Encoder
import PC_Input_Bits_Interleaver
import PC_Subchannel_Allocation
import CRC

# cd Desktop/Thesis/PySageMath/PC
run_vals = [[2/5, 0.14, 15, 'BSC'], [2/5, 0.12, 15, 'BSC'], [2/5, 0.1, 15, 'BSC'], [2/5, 0.08, 21, 'BSC'],
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
    QNF, QNI, MS, matching_scheme = PC_Subchannel_Allocation.freeze(N, K, E, npc, QN0)
    GN = PC_Encoder.gen_G(n)
    QNPC = PC_Subchannel_Allocation.get_n_wm_pc(GN, n_wm_pc, QNI, npc)

    BLER, BER = 0, 0
    print(elem)
    for iteration in range(runs):
        a = random_vector(GF(2), A)
        H = CRC.gen_CRC_mat(K, pol)
        c, G = CRC.CRC(a, A, pol)
        ct = vector(GF(2), list(a) + (vector(GF(2), a) * H[-A:]).list())


        c_ap, PI = PC_Input_Bits_Interleaver.interleaver(I_IL=I_IL, c_seq=c)
        c_ap = vector(GF(2), c_ap)
        iPI = [PI.index(a) for a in PI if a >= A]

        gperm = Matrix(G[a] for a in [y for y in PI if y < A])

        cbits = [c_ap[a] for a in iPI]
        print(f'cbits:\n{cbits}')
        ca = c_ap[:iPI[1]]
        cc = c_ap[:iPI[1]+1]
        hperm = Matrix(GF(2), [H[a] for a in PI])
        print(ca*gperm[:iPI[1]])
        for i in range(hperm.nrows()-G.nrows()):
            hh = gperm[:iPI[0]].rows()
            hh.insert(0, hperm[i])
            res = cc*Matrix(GF(2), hh)
            if res[0]==0:
                print(i)
        #print(ca*hperm[-iPI[0]:])
        #print(cc*hperm[-(iPI[0]+1):])

        breakpoint()




"""
        g1 = gperm[:iPI[0]]
        h1 = g1.rows()
        h1.append(h[A+0])
        h1 = Matrix(GF(2), h1)
        #g1.insert(0, hperm[-14])
        g1 = Matrix(GF(2), g1)
        print((ca*g1)[0])
        print((cc*h1)[0])
        breakpoint()

        cc2 = c_ap[:iPI[1]+1]
        temp = vector(GF(2), [cc2[a] for a in iPI[:2]])
        c2ap = vector(GF(2), [cc2[i] for i, a in enumerate(PI[:iPI[1]]) if a < A])
        c2 = vector(GF(2), list(c2ap)+list(temp))
        g2 = gperm[:iPI[1]-1]
        h2 = hperm[:iPI[1]+1]
        print()
        print((c2ap*g2)[:2])
        print((cc2*h2)[:2])
        breakpoint()
"""