# cd Desktop/Thesis/PySageMath/LDPC
# cd PycharmProjects/Thesis/LDPC

import csv
import os.path
import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt


f_ber = 'BER'
f_bler = 'BLER'

BSC_MS = {((0.3333333333333333, 1.0, 0.0), 'BER'): [[0.021567924528301887, 0.10547311438777597, 0.1676392850600581, 0.22905166021257312, 0.264043331659373], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 1.0, 0.0), 'BER'): [[0.015783941755537326, 0.061415473213609156, 0.09534339040747153, 0.10362780209576868, 0.1567207839695773], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 1.0, 0.0), 'BLER'): [[0.0924, 0.4170141784820684, 0.6934812760055479, 0.8733624454148472, 0.9505703422053232], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 1.0, 0.0), 'BLER'): [[0.1358695652173913, 0.5020080321285141, 0.7710100231303006, 0.6963788300835655, 0.9689922480620154], [0.02, 0.04, 0.06, 0.08, 0.1]]}
BSC_IOMS = {((0.3333333333333333, 0.95, 0.4), 'BER'): [[0.008854716981132075, 0.05347013129206112, 0.11901226930611403, 0.1788144406477492, 0.1907851162474738], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.95, 0.4), 'BER'): [[0.006418867924528302, 0.023932701504217023, 0.05263328717327333, 0.08262325015216068, 0.11317276758542735], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.94, 0.4), 'BER'): [[0.005015094339622641, 0.023892291871590806, 0.04793801328602536, 0.08043371722617006, 0.11207721067733241], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.94, 0.4), 'BER'): [[0.009532075471698114, 0.05562684086203377, 0.10730473279224202, 0.1656217268797852, 0.20703770071325675], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.95, 0.4), 'BLER'): [[0.0422, 0.25012506253126565, 0.546448087431694, 0.7874015748031497, 0.9000900090009001], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.95, 0.4), 'BLER'): [[0.0779, 0.2880184331797235, 0.573394495412844, 0.8064516129032258, 0.9216589861751152], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.94, 0.4), 'BLER'): [[0.0657, 0.2962085308056872, 0.5858230814294083, 0.8116883116883117, 0.9216589861751152], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.94, 0.4), 'BLER'): [[0.0476, 0.2575991756826378, 0.5078720162519045, 0.7818608287724785, 0.931098696461825], [0.02, 0.04, 0.06, 0.08, 0.1]]}
BSC_AMS = {((0.5, 0.95, 0.4), 'BER'): [[0.007156603773584905, 0.028617615114546236, 0.060317838223915246, 0.08555098027164895, 0.11889954674818172], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.95, 0.4), 'BER'): [[0.012166037735849056, 0.06753605418010761, 0.12697121332010833, 0.19117006040226664, 0.22095332671300894], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.94, 0.4), 'BER'): [[0.0058358490566037735, 0.024923372418671253, 0.05253849952342126, 0.08086042460212006, 0.11522105153529708], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.94, 0.4), 'BER'): [[0.013050943396226415, 0.06519873069318993, 0.132227349728372, 0.19262981574539365, 0.21924276008722135], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.95, 0.4), 'BLER'): [[0.0879, 0.33467202141900937, 0.6309148264984227, 0.8143322475570033, 0.931098696461825], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.95, 0.4), 'BLER'): [[0.0567, 0.30367446097783174, 0.5906674542232723, 0.8250825082508251, 0.9398496240601504], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 0.94, 0.4), 'BLER'): [[0.074, 0.3047851264858275, 0.6086427267194157, 0.7824726134585289, 0.9199632014719411], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.3333333333333333, 0.94, 0.4), 'BLER'): [[0.0599, 0.293513354857646, 0.6191950464396285, 0.8375209380234506, 0.9551098376313276], [0.02, 0.04, 0.06, 0.08, 0.1]]}

AWGN_MS = {((0.5, 1.0, 0.0), 'BER'): [[0.1559852383179999, 0.1061239146046876, 0.053732490996377384, 0.017202254168528698, 0.0031764150943396226, 0.0003584905660377358, 5.714285714285714e-05], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.3333333333333333, 1.0, 0.0), 'BER'): [[0.053084116620028224, 0.017223809523809523, 0.003085714285714286, 0.0005666666666666667], [4.0, 5.0, 6.0, 7.0]],
((0.5, 1.0, 0.0), 'BLER'): [[0.9182736455463728, 0.6920415224913494, 0.3739016638624042, 0.12303906490310673, 0.022, 0.0029, 0.0003], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.3333333333333333, 1.0, 0.0), 'BLER'): [[0.25113008538422904, 0.0784, 0.0148, 0.0027], [4.0, 5.0, 6.0, 7.0]]}
AWGN_IOMS = {((0.5, 0.95, 0.4), 'BER'): [[0.09616796798896998, 0.052857862953042685, 0.022314038403420854, 0.006669811320754717, 0.0012283018867924528, 0.00017358490566037735, 5.660377358490566e-06], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.5, 0.95, 0.4), 'BLER'): [[0.7942811755361397, 0.5238344683080147, 0.2390628735357399, 0.0727, 0.0122, 0.0019, 0.0001], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]]}

AWGN_AMS = {((0.5, 0.95, 0.4), 'BER'): [[0.10232100462326628, 0.059493359303930904, 0.026740174972867486, 0.0070849056603773585, 0.0013622641509433962, 0.00013773584905660377, 2.4528301886792453e-05], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.3333333333333333, 0.95, 0.4), 'BER'): [[0.23713383745686809, 0.19580928360091188, 0.11654896173711395, 0.04817222836371812, 0.012964150943396226, 0.0025735849056603773, 0.0002849056603773585], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.5, 0.95, 0.4), 'BLER'): [[0.8278145695364238, 0.5672149744753261, 0.2700513097488523, 0.0743, 0.0133, 0.0012, 0.0003], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.3333333333333333, 0.95, 0.4), 'BLER'): [[0.9727626459143969, 0.856898029134533, 0.5586592178770949, 0.24826216484607747, 0.0675, 0.0133, 0.0016], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]]}
AWGN_OMS = {
((0.5, 1.0, 0.2), 'BER'): [[0.13049278676842255, 0.08331653285580447, 0.04029533912520389, 0.0125, 0.0024669811320754717, 0.00039433962264150943, 3.809523809523809e-05], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.5, 1.0, 0.4), 'BER'): [[0.11750850141891295, 0.07062580741742858, 0.032398727864263746, 0.010110377358490566, 0.0019622641509433963, 0.0003169811320754717], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
((0.5, 1.0, 0.2), 'BLER'): [[0.8783487044356609, 0.6285355122564424, 0.30656039239730226, 0.09715, 0.01815, 0.00275, 0.0002], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
((0.5, 1.0, 0.4), 'BLER'): [[0.8521516829995739, 0.5877167205406993, 0.2767783005812344, 0.0837, 0.01585, 0.0022], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
}

BEC_IOMS = {((0.5, 1.0, 0.0), 'BER'): [[0.0, 0.0, 1.8010291595197256e-05, 0.0003915984336062656, 0.030473034270488395, 0.17410383510183103], [0.1, 0.2, 0.25, 0.3, 0.4, 0.5]],
((0.5, 1.0, 0.0), 'BLER'): [[0.0, 0.0, 0.00042424242424242425, 0.005415094339622641, 0.2660310766833203, 0.9268221996812056], [0.1, 0.2, 0.25, 0.3, 0.4, 0.5]]}


def val_in_row(row, val):
    search_val = row[-2].split(',')
    ind = [i for i, elem in enumerate(search_val) if val in elem]
    return float(search_val[ind[0]].split(':')[1])


def plot_gam_var(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        #next(myfile, None)
        test_results = {}
        #breakpoint()
        for row in myfile:
            if not row[-2].split(',')[-1] == 'Hcore-check-CRC':
                continue
            else:
                gam = val_in_row(row, 'gamma')
                lam = val_in_row(row, 'lam')
                dict_getter = test_results.get((float(row[1]), gam, lam, float(row[7])), [0] * (4))
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


def plot_dicts():
    pl1 = BSC_MS
    pl2 = BSC_IOMS
    pl3 = BSC_AMS
    name = 'BSC'
    ber_val, bler_val = 'BER', 'BLER'
    fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True, facecolor='#F7F7F7')
    fig.set_figheight(7)
    ax1.set_title('BER'); ax1.set_facecolor('#F7F7F7')
    ax2.set_title('BLER'); ax2.set_facecolor('#F7F7F7')
    ax1.semilogy(pl1.get(((0.5, 1.0, 0.0), ber_val))[1], pl1.get(((0.5, 1.0, 0.0), ber_val))[0], label='0.5_MS', linestyle='dotted', marker='|', linewidth=2)
    ax1.semilogy(pl2.get(((0.5, 0.95, 0.4), ber_val))[1], pl2.get(((0.5, 0.95, 0.4), ber_val))[0], label=f'0.5_IOMS', linestyle='dotted', marker='|', linewidth=2)
    ax1.semilogy(pl3.get(((0.5, 0.95, 0.4), ber_val))[1], pl3.get(((0.5, 0.95, 0.4), ber_val))[0], label=f'0.5_AMS', linestyle='dotted', marker='|', linewidth=2)

    ax2.semilogy(pl1.get(((0.5, 1.0, 0.0), bler_val))[1], pl1.get(((0.5, 1.0, 0.0), bler_val))[0], label=f'0.5_MS', linestyle='dotted', marker='|', linewidth=2)
    ax2.semilogy(pl2.get(((0.5, 0.95, 0.4), bler_val))[1], pl2.get(((0.5, 0.95, 0.4), bler_val))[0], label=f'0.5_IOMS', linestyle='dotted',marker='|', linewidth=2)
    ax2.semilogy(pl3.get(((0.5, 0.95, 0.4), bler_val))[1], pl3.get(((0.5, 0.95, 0.4), bler_val))[0], label=f'0.5_AMS', linestyle='dotted', marker='|', linewidth=2)

    ax1.legend()
    ax1.grid(True, linewidth=0.5)
    ax1.set_xlabel('p')
    if name != 'AWGN':
        ax1.invert_xaxis()

    ax2.legend()
    ax2.grid(True, linewidth=0.5)
    ax2.set_xlabel('p')#alpha: '\u03B1'
    if name != 'AWGN':
        ax2.invert_xaxis()
    fig.suptitle(f'{name}')
    fig.savefig(f'LDPC_{name}.svg')


def plot_H(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        #next(myfile, None)
        He, Hc = {}, {}

        for row in myfile:
            if row[-2].split(',')[-1] == 'Hcore-check-CRC':
                dict_getter = Hc.get((float(row[1]), float(row[7])), [0]*5)
            elif row[-2].split(',')[-1] == 'Hextension-check':
                dict_getter = He.get((float(row[1]), float(row[7])), [0]*5)
            else:
                continue
            dict_getter[0] += int(row[0]) * (int(row[4]) + 1)  # tot information bits sent
            dict_getter[1] += int(row[4]) + 1  # tot runs
            dict_getter[2] += int(row[5])  # tot bit errors
            dict_getter[3] += int(row[6])  # tot block errors
            dict_getter[4] += int(row[8])   # tot iterations

            if row[-2].split(',')[-1] == 'Hcore-check-CRC':
                Hc.update({(float(row[1]), float(row[7])): dict_getter})
            elif row[-2].split(',')[-1] == 'Hextension-check':
                He.update({(float(row[1]), float(row[7])): dict_getter})
        for dict in [He, Hc]:
            for key, value in dict.items():
                dict.update({key: [value[0], value[1], value[2]/value[0], value[3]/value[1], value[4]/(value[1]-value[3])]})
        #breakpoint()
        ber_plot, bler_plot, iter_plot = {}, {}, {}
        for key, value in He.items():
            ber_get = ber_plot.get((key[0], 'He'), [[], []])
            bler_get = bler_plot.get((key[0], 'He'), [[], []])

            ber_get[0].append(value[2])
            ber_get[1].append(key[1])
            bler_get[0].append(value[3])
            bler_get[1].append(key[1])
            ber_plot.update({(key[0], 'He'): ber_get})
            bler_plot.update({(key[0], 'He'): bler_get})
        for key, value in Hc.items():
            ber_get = ber_plot.get((key[0], 'Hc'), [[], []])
            bler_get = bler_plot.get((key[0], 'Hc'), [[], []])

            ber_get[0].append(value[2])
            ber_get[1].append(key[1])
            bler_get[0].append(value[3])
            bler_get[1].append(key[1])
            ber_plot.update({(key[0], 'Hc'): ber_get})
            bler_plot.update({(key[0], 'Hc'): bler_get})

        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True, facecolor='#F7F7F7')
        fig.set_figheight(7)
        ax1.set_title('BER'); ax1.set_facecolor('#F7F7F7')
        ax2.set_title('BLER'); ax2.set_facecolor('#F7F7F7')
        #breakpoint()
        for key, value in ber_plot.items():
            plot_x = sorted(value[1], reverse=True)
            plot_y = [value[0][value[1].index(a)] for a in plot_x]
            ax1.semilogy(plot_x, plot_y, label=f'{key[1]}', linestyle='dotted', marker='|', linewidth=2)

        for key, value in bler_plot.items():
            plot_x = sorted(value[1], reverse=True)
            plot_y = [value[0][value[1].index(a)] for a in plot_x]
            ax2.semilogy(plot_x, plot_y, label=key[1], linestyle='dotted', marker='|', linewidth=2)

        ax1.legend()
        ax1.grid(True, linewidth=0.5)
        ax1.set_xlabel('SNR')
        #ax1.invert_xaxis()

        ax2.legend()
        ax2.grid(True, linewidth=0.5)
        ax2.set_xlabel('SNR') #alpha: '\u03B1'
        #ax2.invert_xaxis()

        fig.suptitle(f'He,Hc,AWGN')
        fig.savefig(f'Hc_AWGN.svg')





#update_dicts('Tests/BEC_IOMS.csv')
#plot_gam_var('Tests/AWGN_IOMS.csv')
plot_H('Tests/AWGN_IOMS.csv')
#plot_dicts()

#plot_dict(AWGN_IOMS, 'AWGN_IOMS')
