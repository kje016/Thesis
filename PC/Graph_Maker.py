import csv
import os.path
import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt


f_ber = 'BER'
f_bler = 'BLER'
BSC_SC = {
(0.5, 'BER'): [[4.761904761904762e-05, 0.0007952380952380953, 0.009461904761904762, 0.02034047619047619, 0.035965476190476194, 0.057539682539682536, 0.08086394557823129, 0.14262770562770563], [0.01, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.3333333333333333, 'BER'): [[1e-05, 0.00013, 0.0012666666666666666, 0.006613333333333333, 0.020613333333333334, 0.047516666666666665, 0.059513333333333335, 0.08566333333333333, 0.11824, 0.17614666666666667], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.4, 'BER'): [[6.31578947368421e-05, 0.0009763157894736842, 0.006573684210526316, 0.022594736842105265, 0.05565263157894737, 0.10365, 0.1629157894736842], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]],
(0.5, 'BLER'): [[0.0002, 0.0035, 0.03855, 0.08165, 0.14115, 0.20766666666666667, 0.29482142857142857, 0.4779090909090909], [0.01, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.3333333333333333, 'BLER'): [[0.0001, 0.0009, 0.00555, 0.0248, 0.0724, 0.1491, 0.192, 0.2618, 0.3427, 0.4971], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.4, 'BLER'): [[0.00045, 0.00435, 0.02575, 0.0769, 0.1717, 0.3087, 0.45625], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]]
}
BSC_SCL = {
(0.3333333333333333, 'BER'): [[6.666666666666667e-06, 3.3333333333333335e-05, 0.0006933333333333333, 0.00523, 0.017828387096774194, 0.02666666666666667, 0.04154202898550725, 0.05908205128205128, 0.08079393939393939, 0.14849463115325617, 0.15947333333333333], [0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.5, 'BER'): [[0.0006154761904761905, 0.008140476190476191, 0.01853809523809524, 0.031323809523809525, 0.04823809523809524, 0.07428792912513843, 0.10361635220125787], [0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.4, 'BER'): [[2.6315789473684212e-05, 0.0007302631578947369, 0.005422105263157895, 0.01856637168141593, 0.05019684210526316, 0.09898631578947369, 0.15724], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]],
(0.3333333333333333, 'BLER'): [[0.0001, 0.00025, 0.00355, 0.01905, 0.06813953488372093, 0.091, 0.1367391304347826, 0.1896923076923077, 0.25366666666666665, 0.31783610755441744, 0.466], [0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.5, 'BLER'): [[0.00315, 0.037425, 0.0811, 0.133225, 0.19966666666666666, 0.28967441860465115, 0.39158490566037735], [0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.4, 'BLER'): [[0.00015, 0.00385, 0.02244, 0.06706451612903226, 0.16092, 0.29812, 0.4464], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]]
}
AWGN_SC = {(0.5, 'BER'): [[0.2494591011962496, 0.02754283866795991, 0.0024148076301325574, 0.0018435374149659864, 0.0002227891156462585], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.3333333333333333, 'BER'): [[0.1063288, 0.0081584, 0.0008512, 0.0007666666666666667, 2e-05], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.4, 'BER'): [[0.0989578947368421, 0.03831578947368421, 0.009966623876765083, 0.0015712451861360718, 7.894736842105263e-05], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.5, 'BLER'): [[0.4102121212121212, 0.1390909090909091, 0.035666666666666666, 0.006821428571428571, 0.0009285714285714286], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.3333333333333333, 'BLER'): [[0.2465, 0.063575, 0.0153, 0.0028, 0.00015], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.4, 'BLER'): [[0.2626, 0.1017560975609756, 0.02746341463414634, 0.004585365853658536, 0.0003], [1.0, 2.0, 3.0, 4.0, 5.0]]}
AWGN_SCL = {(0.5, 'BER'): [[0.1082647619047619, 0.042403809523809524, 0.011325714285714286, 0.0017066666666666667, 0.00012476190476190475], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.4, 'BER'): [[0.0981778947368421, 0.038629473684210526, 0.009925263157894737, 0.0015157894736842106, 0.00011157894736842105], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.3333333333333333, 'BER'): [[0.07218666666666666, 0.023843333333333334, 0.00536, 0.00076, 3e-05], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.5, 'BLER'): [[0.34458, 0.14538, 0.04156, 0.00632, 0.00056], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.4, 'BLER'): [[0.25442, 0.10228, 0.02702, 0.00418, 0.00032], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.3333333333333333, 'BLER'): [[0.2039, 0.07135, 0.0169, 0.00265, 0.00015], [1.0, 2.0, 3.0, 4.0, 5.0]]
}
BEC_SC = {(0.5, 'BER'): [[5.476190476190476e-05, 0.0002785714285714286, 0.00185, 0.0075928571428571425, 0.007757142857142857, 0.026116666666666667, 0.06741904761904761, 0.14637380952380952, 0.26320238095238097], [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.4, 'BER'): [[0.0, 8.68421052631579e-05, 0.014089473684210527, 0.010726315789473684, 0.03676842105263158, 0.11582631578947368, 0.20513684210526315, 0.3052342105263158], [0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]],
(0.3333333333333333, 'BER'): [[0.0, 0.0, 0.0, 0.0022616666666666666, 0.043206666666666664, 0.13494666666666666, 0.24866333333333332], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65]],
(0.5, 'BLER'): [[0.0001, 0.00075, 0.0054, 0.02205, 0.0205, 0.076, 0.20215, 0.42395, 0.73405], [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.4, 'BLER'): [[0.0, 0.00025, 0.0317, 0.02405, 0.0811, 0.25865, 0.4708, 0.70775], [0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]],
(0.3333333333333333, 'BLER'): [[0.0, 0.0, 0.0, 0.006325, 0.1095, 0.34485, 0.61255], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65]]
}
BEC_SCL = {(0.3333333333333333, 'BER'): [[0.0, 0.0, 7.777777777777778e-05, 0.00019333333333333333, 0.035330434782608694, 0.07628888888888889, 0.15902898550724637, 0.2437072463768116], [0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65]],
(0.5, 'BER'): [[2.857142857142857e-05, 0.00022857142857142857, 0.0015492063492063492, 0.00680952380952381, 0.013442857142857144, 0.02263015873015873, 0.06831428571428572, 0.12428412698412698, 0.2592892857142857], [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.4, 'BER'): [[0.0, 3.8157894736842105e-05, 0.0016736842105263157, 0.013452631578947369, 0.021573684210526314, 0.12045263157894737, 0.2863894736842105], [0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6]],
(0.3333333333333333, 'BLER'): [[0.0, 0.0, 0.0002, 0.00055, 0.10091304347826087, 0.21733333333333332, 0.43430434782608696, 0.6375652173913043], [0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65]],
(0.5, 'BLER'): [[0.0001, 0.00065, 0.0048, 0.0224, 0.0428, 0.0721, 0.22261666666666666, 0.4004333333333333, 0.75465], [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.4, 'BLER'): [[0.0, 0.000125, 0.004225, 0.03325, 0.052975, 0.29985, 0.7064], [0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6]]
}

rates = [1/2, 2/5]

def update_dicts(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        next(myfile, None)
        test_results, y_axis = {}, {}
        #breakpoint()
        for row in myfile:
            dict_getter = test_results.get((float(row[1]), float(row[8]) ), [0]*(4))
            dict_getter[0]+=int(row[0])*int(row[4]) # tot information bits sent
            dict_getter[1]+= int(row[4])   # tot runs
            dict_getter[2]+= int(row[5])    # tot bit errors
            dict_getter[3]+= int(row[6])   # tot block errors
            test_results.update({(float(row[1]), float(row[8])): dict_getter})
    for key, value in test_results.items():
        test_results.update({key : [value[0], value[1], value[2]/(value[0]), value[3]/value[1]]})

    ber_plot, bler_plot = {}, {}
    for key, value in test_results.items():
        ber_get = ber_plot.get(key[0], [[], []])
        bler_get = bler_plot.get(key[0], [[], []])

        ber_get[0].append(value[2]); ber_get[1].append(key[1])
        bler_get[0].append(value[3]); bler_get[1].append(key[1])

        ber_plot.update({key[0] : ber_get})
        bler_plot.update({key[0] : bler_get})

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


def plot_dict(plotter, name):
    #breakpoint()
    fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True, facecolor='#F7F7F7')
    ax1.set_title('BER')
    ax2.set_title('BLER')
    for rate in rates:
        ax1.semilogy(plotter.get((rate, 'BER'))[1], plotter.get((rate, 'BER'))[0], label=f'{rate}', linestyle='dotted', marker='|', markerfacecolor='none')
        ax2.semilogy(plotter.get((rate, 'BLER'))[1], plotter.get((rate, 'BLER'))[0], label=f'{rate}', linestyle='dotted',
                     marker='|', markerfacecolor='none')
    ax1.legend()
    ax1.grid(True, linewidth=0.5)
    ax1.set_xlabel('SNR')

    if name.split('_')[0] != 'AWGN':
        ax1.invert_xaxis()

    ax2.legend()
    ax2.grid(True, linewidth=0.5)
    ax2.set_xlabel('SNR')
    if name.split('_')[0] != 'AWGN':
        ax2.invert_xaxis()

    fig.suptitle(f'Performance of {name}')
    fig.savefig(f'{name}.svg')


def plot_SC_SCL():
    plotter1 = BEC_SC
    plotter2 = BEC_SCL
    name = 'BEC'
    ber_val = 'BER'
    bler_val = 'BLER'
    fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True, facecolor='#F7F7F7')
    ax1.set_title('BER')
    ax1.set_facecolor('#F7F7F7')
    ax2.set_title('BLER')
    ax2.set_facecolor('#F7F7F7')
    #fig, ax1 = plt.subplots()
    fig.suptitle(f'Performance of {name}')
    #ax1.set_facecolor('#F7F7F7'); ax2.set_facecolor('#F7F7F7'); ax3.set_facecolor('#F7F7F7');


    #Plots comparison between SC and SCL, each ax is of same rate
    for rate in rates:
        ax1.semilogy(plotter1.get((rate, ber_val))[1], plotter1.get((rate, ber_val))[0], label=f'{rate}_SC', linestyle='dotted', marker='|', markerfacecolor='none')
        ax1.semilogy(plotter2.get((rate, ber_val))[1], plotter2.get((rate, ber_val))[0], label=f'{rate}_SCL', linestyle='dashed', marker='|', markerfacecolor='none')

    for rate in rates:
        ax2.semilogy(plotter1.get((rate, bler_val))[1], plotter1.get((rate, bler_val))[0], label=f'{rate}_SC', linestyle='dotted', marker='|', markerfacecolor='none')
        ax2.semilogy(plotter2.get((rate, bler_val))[1], plotter2.get((rate, bler_val))[0], label=f'{rate}_SCL', linestyle='dotted', marker='|', markerfacecolor='none')


    ax1.legend()
    ax1.grid(True, linewidth=0.5)
    ax1.set_xlabel('SNR')
    if name.split('_')[0] != 'AWGN':
        ax1.invert_xaxis()

    ax2.legend()
    ax2.grid(True, linewidth=0.5)
    ax2.set_xlabel('SNR')
    if name.split('_')[0] != 'AWGN':
        ax2.invert_xaxis()
    #ax3.legend()
    #ax3.grid(True, linewidth=0.5)
    #ax3.set_xlabel('SNR')
    #ax3.set_ylabel(f'{ber_val}')
    # BERax.set_title(f'Bit Error Rate {name}')
    #if name!= 'AWGN':
    #    ax3.invert_xaxis()
    #ax3.legend()

    fig.savefig(f'comparison_{name}.svg', facecolor=fig.get_facecolor())

#update_dicts('C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\BEC_SC.csv')
#plot_SC_SCL()
plot_dict(BEC_SC, 'BEC_SC')
