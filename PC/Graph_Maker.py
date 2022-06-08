import csv
import os.path
import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt


f_ber = 'BER'
f_bler = 'BLER'
BSC_SC = {
(0.5, 'BER'): [[4.761904761904762e-05, 0.0007952380952380953, 0.009609523809523809, 0.019896774193548388, 0.03598571428571429, 0.057539682539682536, 0.08086394557823129, 0.13613599274705349], [0.01, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.3333333333333333, 'BER'): [[1e-05, 0.00013, 0.0012666666666666666, 0.006613333333333333, 0.020613333333333334, 0.047516666666666665, 0.059513333333333335, 0.08566333333333333, 0.19947083333333335, 0.17614666666666667], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.4, 'BER'): [[6.31578947368421e-05, 0.0009763157894736842, 0.006573684210526316, 0.022594736842105265, 0.05565263157894737, 0.10365, 0.1629157894736842], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]],
(0.5, 'BLER'): [[0.0002, 0.0035, 0.0389, 0.0824, 0.14063333333333333, 0.20766666666666667, 0.29482142857142857, 0.4699622641509434], [0.01, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.3333333333333333, 'BLER'): [[0.0001, 0.0009, 0.00555, 0.0248, 0.0724, 0.1491, 0.192, 0.2618, 0.38236363636363635, 0.4971], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.4, 'BLER'): [[0.00045, 0.00435, 0.02575, 0.0769, 0.1717, 0.3087, 0.45625], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]]
}
BSC_SCL = {
(0.5, 'BER'): [[0.0011942622950819673, 0.020195901639344264, 0.01853809523809524, 0.09575901639344263, 0.04823809523809524, 0.24061808261886206, 0.2955250886974151], [0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.3333333333333333, 'BER'): [[6.666666666666667e-06, 3.3333333333333335e-05, 0.0006933333333333333, 0.00523, 0.01815157894736842, 0.02666666666666667, 0.042022222222222225, 0.05908205128205128, 0.07622222222222222, 0.14849463115325617, 0.15947333333333333], [0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.4, 'BER'): [[2.6315789473684212e-05, 0.001817391304347826, 0.011907608695652173, 0.04612735426008969, 0.1417358695652174, 0.2831489130434783, 0.1582078947368421], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]],
(0.5, 'BLER'): [[0.002425, 0.0318, 0.0811, 0.13015, 0.19966666666666666, 0.314046511627907, 0.4393717948717949], [0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]],
(0.3333333333333333, 'BLER'): [[0.0001, 0.00025, 0.00355, 0.01905, 0.07369565217391304, 0.091, 0.144, 0.1896923076923077, 0.25066666666666665, 0.31783610755441744, 0.466], [0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.17]],
(0.4, 'BLER'): [[0.00015, 0.003375, 0.019175, 0.07071698113207547, 0.183525, 0.354475, 0.4495], [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]]
}
AWGN_SC = {
(0.5, 'BER'): [[0.23858149130562922, 0.054572354211663066, 0.028334345581536594, 0.014788336933045357, 0.0028985727300334043, 0.002853131749460043, 0.0018147208121827411, 0.00034989200863930886, 0.00022081218274111675, 1.9438444924406046e-05], [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]],
(0.3333333333333333, 'BER'): [[0.10937623188405797, 0.006826400443704936, 0.0003547826086956522, 0.0010333333333333334], [1.0, 2.0, 3.0, 4.0]],
(0.4, 'BER'): [[0.0989578947368421, 0.06486052631578948, 0.03831578947368421, 0.021247368421052633, 0.009966623876765083, 0.0044631578947368425, 0.0015712451861360718, 0.0005052631578947368, 7.894736842105263e-05, 2.6315789473684212e-05], [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]],
(0.5, 'BLER'): [[0.41870833333333335, 0.1822608695652174, 0.14393023255813953, 0.05613043478260869, 0.03806976744186046, 0.011521739130434782, 0.007210526315789474, 0.0015217391304347826, 0.0009473684210526315, 8.695652173913044e-05], [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]],
(0.3333333333333333, 'BLER'): [[0.2873, 0.05161904761904762, 0.0104, 0.0038], [1.0, 2.0, 3.0, 4.0]],
(0.4, 'BLER'): [[0.2626, 0.17415, 0.1017560975609756, 0.05755, 0.02746341463414634, 0.01285, 0.004585365853658536, 0.0014, 0.0003, 0.0001], [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]]}
AWGN_SCL = {(0.5, 'BER'): [[0.1085904761904762, 0.04210714285714286, 0.012019047619047618, 0.0018071428571428572, 0.00010238095238095239], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.4, 'BER'): [[0.09893947368421052, 0.039086842105263156, 0.009494736842105262, 0.0014157894736842105, 8.157894736842105e-05], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.5, 'BLER'): [[0.34595, 0.14495, 0.0443, 0.00675, 0.00045], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.4, 'BLER'): [[0.2569, 0.1033, 0.0266, 0.00385, 0.0002], [1.0, 2.0, 3.0, 4.0, 5.0]]
}
BEC_SC = {(0.5, 'BER'): [[5e-06, 0.00142, 0.015523170731707317, 0.026116666666666667, 0.06741904761904761, 0.14637380952380952, 0.26320238095238097], [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.4, 'BER'): [[0.0, 8.68421052631579e-05, 0.010726315789473684, 0.03676842105263158, 0.11582631578947368, 0.20513684210526315, 0.3052342105263158], [0.1, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6]],
(0.25, 'BER'): [[0.00039523809523809526, 0.0028809523809523807, 0.01984285714285714, 0.08881428571428572, 0.21581666666666666], [0.5, 0.55, 0.6, 0.65, 0.7]],
(0.3333333333333333, 'BER'): [[0.0, 0.0, 0.0, 0.00018, 0.10839333333333333, 0.22676666666666667], [0.1, 0.2, 0.3, 0.4, 0.6, 0.65]],
(0.5, 'BLER'): [[0.0001, 0.0103, 0.0797, 0.076, 0.20215, 0.42395, 0.73405], [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.4, 'BLER'): [[0.0, 0.00025, 0.02405, 0.0811, 0.25865, 0.4708, 0.70775], [0.1, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6]],
(0.25, 'BLER'): [[0.00105, 0.00775, 0.04795, 0.2201, 0.531], [0.5, 0.55, 0.6, 0.65, 0.7]],
(0.3333333333333333, 'BLER'): [[0.0, 0.0, 0.0, 0.00115, 0.2796, 0.5751], [0.1, 0.2, 0.3, 0.4, 0.6, 0.65]]
}
BEC_SCL = {(0.5, 'BER'): [[0.0, 0.000925, 0.013380673499267935, 0.02545165945165945, 0.06235669362084456, 0.13284383954154727, 0.2713997477931904], [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.3333333333333333, 'BER'): [[0.0, 0.0, 0.0, 0.00019333333333333333, 0.03582222222222222, 0.07628888888888889, 0.15222222222222223, 0.2495111111111111], [0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65]],
(0.25, 'BER'): [[0.00016428571428571428, 0.0029261904761904763, 0.017945238095238094, 0.07381190476190476, 0.1888595238095238], [0.5, 0.55, 0.6, 0.65, 0.7]],
(0.4, 'BER'): [[0.0, 0.0, 4.736842105263158e-05, 0.00923421052631579, 0.09314473684210527, 0.26403947368421055], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]],
(0.5, 'BLER'): [[0.0, 0.0053, 0.05412121212121212, 0.08163636363636363, 0.20471698113207548, 0.41700906344410876, 0.768], [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5]],
(0.3333333333333333, 'BLER'): [[0.0, 0.0, 0.0, 0.00055, 0.10233333333333333, 0.21733333333333332, 0.41833333333333333, 0.643], [0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65]],
(0.25, 'BLER'): [[0.00065, 0.0081, 0.04745, 0.1984, 0.51045], [0.5, 0.55, 0.6, 0.65, 0.7]],
(0.4, 'BLER'): [[0.0, 0.0, 0.00025, 0.0221, 0.23415, 0.6667], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]]
}


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


def csv_to_plot(input_file):
    name = input_file.split('\\')[-1].replace('.csv', '')
    ber_plot, bler_plot = update_dicts(input_file)
    fig, (BERax, BLERax) = plt.subplots(2, constrained_layout = True)
    fig.suptitle(f'Performance of {name}')
    BERax.set_facecolor('#F7F7F7')
    #breakpoint()
    for key, value in ber_plot.items():
        plot_x = sorted(value[1])
        plot_y = [value[0][value[1].index(a)] for a in plot_x]
        #plot_x, plot_y = prep_axis(plot_x, plot_y, name.split('_')[0])
        BERax.semilogy(plot_x, plot_y, label=str(key), marker='d', markerfacecolor='none')

    BERax.grid(True, linewidth=0.5)
    BERax.set_xlabel('SNR')
    BERax.set_ylabel('BER')
    #BERax.set_title(f'Bit Error Rate {name}')
    if name.split('_')[0] != 'AWGN':
        BERax.invert_xaxis()
    BERax.legend()    # show a legend on the plot
    #ax.show() # function to show the plot
    #fig.savefig(f'{name}_BER.svg')
    #breakpoint()

    BLERax.set_facecolor('#F7F7F7')
    for key, value in bler_plot.items():
        plot_x = sorted(value[1])
        plot_y = [value[0][value[1].index(a)] for a in plot_x]
        #plot_x, plot_y = prep_axis(plot_x, plot_y, name.split('_')[0])
        BLERax.semilogy(plot_x, plot_y, label=str(key), marker='d', markerfacecolor='none')
        #BLERax.set(title=str(key))

    BLERax.grid(True, linewidth=0.5)
    BLERax.set_xlabel('SNR')
    BLERax.set_ylabel('BLER')
    if name.split('_')[0] != 'AWGN':
        BLERax.invert_xaxis()
    #BLERax.set_title(f'Block Error Rate {name}')
    BLERax.legend()    # show a legend on the plot
    #ax.show() # function to show the plot
    #fig.savefig(f'{name}_BLER.svg')
    #plt.legend()    # show a legend on the plot
    #plt.show() # function to show the plot
    fig.subplots_adjust(wspace=1/2)
    fig.savefig(f'{name}.svg')


def plot_dict():
    plotter1 = AWGN_SC
    plotter2 = AWGN_SCL
    name = 'AWGN'
    test_val = 'BER'
    fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
    #fig, ax1 = plt.subplots()
    fig.suptitle(f'Performance of {name}')
    #ax1.set_facecolor('#F7F7F7'); ax2.set_facecolor('#F7F7F7'); ax3.set_facecolor('#F7F7F7');
    ax1.set_facecolor('#F7F7F7')
    rates = [1/2, 2/5, 1/3]
    #breakpoint()
    ax1.semilogy(plotter1.get((rates[0], test_val))[1],plotter1.get((rates[0], test_val))[0], label=f'{rates[0]}_SC', linestyle='dotted', marker='|', markerfacecolor='none')
    ax1.semilogy(plotter2.get((rates[0], test_val))[1],plotter2.get((rates[0], test_val))[0], label=f'{rates[0]}_SCL', linestyle='dashed', marker='|', markerfacecolor='none')
    #ax1.semilogy(plotter1.get((rates[0], 'BLER')), label=f'{rates[0]}_BLER', marker='o', markerfacecolor='none')

    ax2.semilogy(plotter1.get((rates[1], test_val))[1],plotter1.get((rates[1], test_val))[0], label=f'{rates[1]}_SC', linestyle='dotted', marker='|', markerfacecolor='none')
    ax2.semilogy(plotter2.get((rates[1], test_val))[1], plotter2.get((rates[1], test_val))[0],label=f'{rates[1]}_SCL', linestyle='dashed', marker='|',markerfacecolor='none')

    #ax3.semilogy(plotter1.get((rates[2], test_val))[1], plotter1.get((rates[2], test_val))[0],label=f'{rates[2]}_SC', linestyle='dotted', marker='|', markerfacecolor='none')
    #ax3.semilogy(plotter2.get((rates[2], test_val))[1], plotter2.get((rates[2], test_val))[0],label=f'{rates[2]}_SCL', linestyle='dashed', marker='|', markerfacecolor='none')

    ax1.legend()
    ax1.grid(True, linewidth=0.5)
    ax1.set_xlabel('SNR')
    #ax1.set_ylabel(f'{test_val}')
    #BERax.set_title(f'Bit Error Rate {name}')
    if name.split('_')[0] != 'AWGN':
        ax1.invert_xaxis()

    ax2.legend()
    ax2.grid(True, linewidth=0.5)
    ax2.set_xlabel('SNR')
    #ax2.set_ylabel(f'{test_val}')
    # BERax.set_title(f'Bit Error Rate {name}')
    if name.split('_')[0] != 'AWGN':
        ax2.invert_xaxis()
    #ax3.legend()
    #ax3.grid(True, linewidth=0.5)
    #ax3.set_xlabel('SNR')
    #ax3.set_ylabel(f'{test_val}')
    # BERax.set_title(f'Bit Error Rate {name}')
    #if name!= 'AWGN':
    #    ax3.invert_xaxis()
    #ax3.legend()

    fig.savefig(f'comparison_{name}.svg')

update_dicts('C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\AWGN_SCL.csv')
plot_dict()

#csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/AWGN_SC.csv")))
#res = csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/BSC_SC.csv")))
#res = csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/BEC_SC.csv")))
#csv_to_plot('C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\BEC_SCL.csv')
#csv_to_plot(open(os.path.join("C:/Users/Kristian/Desktop/Thesis/PySageMath/PC/Tests", "AWGN_SCL.csv")))
