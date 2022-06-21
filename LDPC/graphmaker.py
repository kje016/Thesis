import csv
import os.path
import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt


f_ber = 'BER'
f_bler = 'BLER'

BSC_MS = {(0.3333333333333333, 'BER'): [[0.02157008153645553, 0.10551711652143796, 0.16775562044177916, 0.22925188019527643, 0.2642945622318367], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.5, 'BER'): [[0.01578608660426073, 0.06144631975967325, 0.09541695783834149, 0.12555683596773362, 0.1568727924894314], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.3333333333333333, 'BLER'): [[0.0924092409240924, 0.4171881518564873, 0.6939625260235948, 0.8741258741258742, 0.9514747859181731], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.5, 'BLER'): [[0.1358880282647099, 0.5022601707684581, 0.7716049382716049, 0.9115770282588879, 0.9699321047526673], [0.02, 0.04, 0.06, 0.08, 0.1]]}
BSC_IOMS = {(0.3333333333333333, 'BER'): [[0.008855602541386214, 0.05348350885805864, 0.11907733889020704, 0.17895535037245194, 0.19095699473057964], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.5, 'BER'): [[0.006204394024308091, 0.02515166694495464, 0.051917047142307056, 0.08211490626706304, 0.11761414307447497], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.3333333333333333, 'BLER'): [[0.04220422042204221, 0.2501876407305479, 0.5467468562055768, 0.7880220646178093, 0.9009009009009009], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.5, 'BLER'): [[0.0697069706970697, 0.2949852507374631, 0.5830903790087464, 0.8038585209003215, 0.9267840593141798], [0.02, 0.04, 0.06, 0.08, 0.1]]}
BSC_AMS = {(0.5, 'BER'): [[0.0071573195055354595, 0.028627195836044242, 0.06035591766723842, 0.08562070397195184, 0.11901035713657705], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.3333333333333333, 'BER'): [[0.012167254461295187, 0.06755656938490108, 0.12704625540835898, 0.19132792172382096, 0.22116118496955928], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.5, 'BLER'): [[0.08790879087908791, 0.33478406427854035, 0.6313131313131313, 0.8149959250203749, 0.9319664492078286], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.3333333333333333, 'BLER'): [[0.05670567056705671, 0.3037667071688943, 0.5910165484633569, 0.8257638315441783, 0.940733772342427], [0.02, 0.04, 0.06, 0.08, 0.1]]}

AWGN_MS = {(0.5, 'BER'): [[0.15605688978254653, 0.1061606483930589, 0.05374253820860558, 0.017203312508270823, 0.0031765739230357744, 0.00035850849146230895, 5.714857200005715e-05], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3, 'BER'): [[0.29753480197513416, 0.2461208333689547, 0.1609747771168971, 0.06849214087782257, 0.019578256282653195, 0.004168132934948634, 0.001162020964001162], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3333333333333333, 'BER'): [[0.05309745098742838, 0.017225532077017224, 0.003086022888003086, 0.0005667233390005668], [4.0, 5.0, 6.0, 7.0]],
(0.5, 'BLER'): [[0.9186954524575104, 0.6922810661128418, 0.3739715781600598, 0.12304663467454165, 0.02200110005500275, 0.0029001450072503624, 0.00030003000300030005], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3, 'BLER'): [[0.9857072449482503, 0.9062075215224287, 0.6958942240779401, 0.3551767004084532, 0.11138338159946536, 0.023651182559127956, 0.005200520052005201], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3333333333333333, 'BLER'): [[0.25119316754584275, 0.0784078407840784, 0.014801480148014802, 0.0027002700270027003], [4.0, 5.0, 6.0, 7.0]]}

AWGN_IOMS = {((0.3, 0.95, 0.4), 'BER'): [[0.21900815302884571, 0.1770873907420968, 0.10496929378283035, 0.039062890906282566, 0.010961320754716982, 0.0019622641509433963], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.3, 0.95, 0.2), 'BER'): [[0.2526393997486723, 0.1955661188533234, 0.12425648549364814, 0.05060713359035715, 0.013043396226415094, 0.0024528301886792454], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.5, 0.95, 0.4), 'BER'): [[0.09590442658963017, 0.05569365771492998, 0.022243297851660297, 0.006584905660377358, 0.0011232704402515724, 0.0001858490566037736], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.5, 0.95, 0.2), 'BER'): [[0.11075632962425415, 0.0661017076713504, 0.029037475002173723, 0.008665408805031447, 0.0012830188679245284, 0.00035943396226415094], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.3, 0.95, 0.4), 'BLER'): [[0.9528346831824679, 0.8257638315441783, 0.5474952094169176, 0.23060071486221606, 0.0672, 0.013], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.3, 0.95, 0.2), 'BLER'): [[0.9722897423432183, 0.8598452278589854, 0.6027727546714888, 0.27781636338380333, 0.078, 0.0159], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.5, 0.95, 0.4), 'BLER'): [[0.7974481658692185, 0.5399568034557235, 0.23415547923821417, 0.06986666666666666, 0.011733333333333333, 0.00155], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.5, 0.95, 0.2), 'BLER'): [[0.8547008547008547, 0.5983246908655764, 0.2764976958525346, 0.08026666666666667, 0.011666666666666667, 0.0031], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]}
AWGN_AMS = {(0.5, 'BER'): [[0.15605688978254653, 0.1061606483930589, 0.05374253820860558, 0.017203312508270823, 0.0031765739230357744, 0.00035850849146230895, 5.714857200005715e-05], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3, 'BER'): [[0.29753480197513416, 0.2461208333689547, 0.1609747771168971, 0.06849214087782257, 0.019578256282653195, 0.004168132934948634, 0.001162020964001162], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3333333333333333, 'BER'): [[0.05309745098742838, 0.017225532077017224, 0.003086022888003086, 0.0005667233390005668], [4.0, 5.0, 6.0, 7.0]],
(0.5, 'BLER'): [[0.9186954524575104, 0.6922810661128418, 0.3739715781600598, 0.12304663467454165, 0.02200110005500275, 0.0029001450072503624, 0.00030003000300030005], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3, 'BLER'): [[0.9857072449482503, 0.9062075215224287, 0.6958942240779401, 0.3551767004084532, 0.11138338159946536, 0.023651182559127956, 0.005200520052005201], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
(0.3333333333333333, 'BLER'): [[0.25119316754584275, 0.0784078407840784, 0.014801480148014802, 0.0027002700270027003], [4.0, 5.0, 6.0, 7.0]],}


def update_dicts(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        #next(myfile, None)
        test_results = {}
        breakpoint()
        for row in myfile:
            dict_getter = test_results.get((float(row[1]), float(row[8])), [0] * (4))
            dict_getter[0] += int(row[0]) * int(row[4])  # tot information bits sent
            dict_getter[1] += int(row[4])  # tot runs
            dict_getter[2] += int(row[5])  # tot bit errors
            dict_getter[3] += int(row[6])  # tot block errors
            test_results.update({(float(row[1]), float(row[7])): dict_getter})
        for key, value in test_results.items():
            test_results.update({key: [value[0], value[1], value[2] / (value[0]), value[3] / value[1]]})
        ber_plot, bler_plot = {}, {}
        for key, value in test_results.items():
            ber_get = ber_plot.get(key[0], [[], []])
            bler_get = bler_plot.get(key[0], [[], []])

            ber_get[0].append(value[2]);
            ber_get[1].append(key[1])
            bler_get[0].append(value[3]);
            bler_get[1].append(key[1])

            ber_plot.update({key[0]: ber_get})
            bler_plot.update({key[0]: bler_get})

        temp_dict = {}
        for key, value in ber_plot.items():
            plot_x = sorted(value[1])
            plot_y = [value[0][value[1].index(a)] for a in plot_x]
            temp_dict.update({(key, 'BER'): [plot_y, plot_x]})
            print(f'({key}, {f_ber!r}): [{plot_y}, {plot_x}],')
        for key, value in bler_plot.items():
            plot_x = sorted(value[1])
            plot_y = [value[0][value[1].index(a)] for a in plot_x]
            temp_dict.update({(key, 'BLER'): [plot_y, plot_x]})
            print(f'({key}, {f_bler!r}): [{plot_y}, {plot_x}],')
        return ber_plot, bler_plot


def plot_gam_var(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        #next(myfile, None)
        test_results = {}
        #breakpoint()
        for row in myfile:
            gam = float(row[-2].split(',')[1].split(':')[1])
            lam = float(row[-2].split(',')[-1][-3:])
            dict_getter = test_results.get((float(row[1]), gam, lam), [0] * (4))
            dict_getter[0] += int(row[0]) * (int(row[4])+1)  # tot information bits sent
            dict_getter[1] += int(row[4])+1  # tot runs
            dict_getter[2] += int(row[5])  # tot bit errors
            dict_getter[3] += int(row[6])  # tot block errors
            test_results.update({(float(row[1]), gam, lam, float(row[7])): dict_getter})
        #breakpoint()
        for key, value in test_results.items():
            test_results.update({key: [value[0], value[1], value[2] / (value[0]), value[3] / value[1]]})

        ber_plot, bler_plot = {}, {}
        #breakpoint()
        for key, value in test_results.items():
            ber_get = ber_plot.get((key[0], key[1], key[2]), [[], []])
            bler_get = bler_plot.get((key[0], key[1], key[2]), [[], []])

            ber_get[0].append(value[2]);
            ber_get[1].append(key[3])
            bler_get[0].append(value[3]);
            bler_get[1].append(key[3])
            ber_plot.update({(key[0], key[1], key[2]): ber_get})
            bler_plot.update({(key[0], key[1], key[2]): bler_get})

        temp_dict = {}
        #breakpoint()
        for key, value in ber_plot.items():
            plot_x = sorted(value[1])
            plot_y = [value[0][value[1].index(a)] for a in plot_x]
            temp_dict.update({(key, 'BER'): [plot_y, plot_x]})
            print(f'({key}, {f_ber!r}): [{plot_y}, {plot_x}],')
        for key, value in bler_plot.items():
            plot_x = sorted(value[1])
            plot_y = [value[0][value[1].index(a)] for a in plot_x]
            temp_dict.update({(key, 'BLER'): [plot_y, plot_x]})
            print(f'({key}, {f_bler!r}): [{plot_y}, {plot_x}],')
        return ber_plot, bler_plot


def plot_dict():
    pl1 = AWGN_MS
    pl2 = BSC_IOMS
    pl3 = AWGN_AMS
    name = 'AWGN'
    ber_val, bler_val = 'BER', 'BLER'
    fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True, facecolor='#F7F7F7')
    ax1.set_title('BER'); ax1.set_facecolor('#F7F7F7')
    ax2.set_title('BLER'); ax2.set_facecolor('#F7F7F7')

    ax1.semilogy(pl1.get((0.5, ber_val))[1], pl1.get((0.5, ber_val))[0], label='0.5_MS', linestyle='dotted', marker='|')
    #ax1.semilogy(pl2.get((0.5, ber_val))[1], pl2.get((0.5, ber_val))[0], label=f'0.5_IOMS', linestyle='dotted', marker='|')
    ax1.semilogy(pl3.get((0.5, ber_val))[1], pl3.get((0.5, ber_val))[0], label=f'0.5_AMS', linestyle='dotted', marker='|')

    ax2.semilogy(pl1.get((0.5, bler_val))[1], pl1.get((0.5, bler_val))[0], label=f'0.5_MS', linestyle='dotted', marker='|')
    #ax2.semilogy(pl2.get((0.5, bler_val))[1], pl2.get((0.5, bler_val))[0], label=f'0.5_IOMS', linestyle='dotted',marker='|')
    ax2.semilogy(pl3.get((0.5, bler_val))[1], pl3.get((0.5, bler_val))[0], label=f'0.5_AMS', linestyle='dotted', marker='|')

    ax1.legend()
    ax1.grid(True, linewidth=0.5)
    ax1.set_xlabel('SNR')
    if name != 'AWGN':
        ax1.invert_xaxis()

    ax2.legend()
    ax2.grid(True, linewidth=0.5)
    ax2.set_xlabel('SNR')
    if name != 'AWGN':
        ax2.invert_xaxis()
    fig.suptitle(f'{name}')
    fig.savefig(f'{name}.svg')

#update_dicts('Tests/AWGN_MS.csv')
plot_gam_var('Tests/AWGN_IOMS.csv')
#plot_dict()
