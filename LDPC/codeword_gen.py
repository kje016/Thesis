from sage.all import *
import csv

A = 21
snr = elem[1]
rate = elem[2]
channel = elem[3]
N0, sig = None, None
pol = CRC.get_pol(A)
B = A + pol.degree()

# sigma = vector(RealField(10), map(lambda z: sqrt(1 / (2 * R * 10 ** (z / 10))), SNR))
# N0 = 2 * sigma[0] ** 2
bg = PF.det_BG(A, rate)
L, C, B_ap = PF.get_code_block_param(bg=bg, B=B)
K_ap = B_ap // C
Kb = PF.determine_kb(B=B, bg=bg)
Zc, iLS, K = PF.det_Z(bg=bg, kb=Kb, lifting_set=lss, K_ap=K_ap)
BG = HF.get_base_matrix(bg, iLS, Zc)
BGB = BG.matrix_from_rows_and_columns(list(range(4)), list(range(10, 10 + 4)))
print(bg, iLS, Zc)
H = HF.Protograph(BG, Zc)
del BG;
gc.collect()

runs = 60000
breakpoint()
for iteration in range(runs):
    a = random_vector(GF(2), A)
    c, G = CRC.CRC(a, A, pol)
    # crk := padding the codeword
    crk = PF.calc_crk(C=C, K=K, K_ap=K_ap, L=L, b_bits=c)  # TODO: testing for C > 1 & need to split crk
    D = vector(GF(2), crk)
    u = LDPC_Encoding.Encoding(H=H, Zc=Zc, D=D, K=K, kb=Kb, BG=bg)
    with open(f'C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\LDPC\\u_codewords.txt', mode='a', newline='') as file:
        result_writer = csv.writer(file)
        result_writer.writerow([A, bg, u])