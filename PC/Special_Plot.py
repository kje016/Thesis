import csv
import os.path
import matplotlib as mpl

import Graph_Maker

mpl.use('template')
import matplotlib.pyplot as plt


f_ber = 'BER'
f_bler = 'BLER'

CA_BSC = {(0.5, 'BER'): [[0.0003875, 0.005459375, 0.02128125, 0.05409375, 0.1019546875], [0.02, 0.04, 0.06, 0.08, 0.1]],
(0.5, 'BLER'): [[0.001175, 0.014575, 0.0561, 0.1368, 0.2496625], [0.02, 0.04, 0.06, 0.08, 0.1]]}
CA_BEC = {(0.5, 'BER'): [[0.0, 0.0007203125, 0.017209375, 0.118646875, 0.325384375], [0.1, 0.2, 0.3, 0.4, 0.5]],
(0.5, 'BLER'): [[0.0, 0.0011625, 0.030033333333333332, 0.20846666666666666, 0.5795333333333333], [0.1, 0.2, 0.3, 0.4, 0.5]]}
CA_AWGN = {(0.5, 'BER'): [[0.08906354166666666, 0.0335890625, 0.008734375, 0.0014765625, 0.000140625], [1.0, 2.0, 3.0, 4.0, 5.0]],
(0.5, 'BLER'): [[0.20565833333333333, 0.0786875, 0.0203625, 0.0033875, 0.000325], [1.0, 2.0, 3.0, 4.0, 5.0]]}

RM_BSC = {((0.5, 40), 'BER'): [[0.0013856060606060607, 0.01623409090909091, 0.05296818181818182, 0.10413409090909091, 0.15980606060606062], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 56), 'BER'): [[0.001736, 0.016485, 0.052602, 0.102607, 0.157982], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 64), 'BER'): [[0.0011642857142857143, 0.012436904761904762, 0.03927380952380952, 0.08128333333333333, 0.13030952380952382], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 24), 'BER'): [[0.0005871951219512195, 0.007696951219512195, 0.033053658536585366, 0.08308963414634146, 0.1493939024390244], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 8), 'BER'): [[0.0002535714285714286, 0.0045448979591836735, 0.026365816326530612, 0.078575, 0.1568622448979592], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 40), 'BLER'): [[0.00995, 0.104025, 0.3018, 0.531475, 0.72265], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 56), 'BLER'): [[0.0122, 0.098875, 0.276225, 0.48725, 0.6749], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 64), 'BLER'): [[0.00805, 0.076825, 0.218575, 0.41055, 0.60005], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 24), 'BLER'): [[0.005775, 0.0613, 0.21125, 0.431775, 0.6527], [0.02, 0.04, 0.06, 0.08, 0.1]],
((0.5, 8), 'BLER'): [[0.002325, 0.0341, 0.1458, 0.354775, 0.60475], [0.02, 0.04, 0.06, 0.08, 0.1]]}
RM_BEC = {((0.5, 40), 'BER'): [[1.2121212121212122e-05, 0.003381060606060606, 0.03773030303030303, 0.14262424242424243, 0.2628530303030303], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 56), 'BER'): [[0.000103, 0.008637, 0.063618, 0.175329, 0.283428], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 64), 'BER'): [[8.214285714285714e-05, 0.006889285714285714, 0.04746071428571429, 0.13954880952380952, 0.24293452380952382], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 24), 'BER'): [[8.536585365853658e-06, 0.0010213414634146342, 0.022687804878048782, 0.12999268292682928, 0.2917768292682927], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 8), 'BER'): [[0.0, 0.0002260204081632653, 0.015720918367346938, 0.15950663265306123, 0.354740306122449], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 40), 'BLER'): [[5e-05, 0.01325, 0.167025, 0.60765, 0.94], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 56), 'BLER'): [[0.000325, 0.027675, 0.227525, 0.64515, 0.936925], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 64), 'BLER'): [[0.000225, 0.02075, 0.1709, 0.5381, 0.8845], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 24), 'BLER'): [[2.5e-05, 0.004725, 0.0919, 0.473025, 0.9039], [0.1, 0.2, 0.3, 0.4, 0.5]],
((0.5, 8), 'BLER'): [[0.0, 0.002025, 0.0901, 0.6053, 0.973275], [0.1, 0.2, 0.3, 0.4, 0.5]]}
RM_AWGN = {((0.5, 24), 'BER'): [[0.10669268292682926, 0.04072012195121951, 0.008742073170731708, 0.00105, 7.682926829268293e-05], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 40), 'BER'): [[0.1175189393939394, 0.05333636363636363, 0.016566666666666667, 0.0031015151515151516, 0.00027045454545454546], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 56), 'BER'): [[0.129582, 0.066066, 0.023976, 0.005652, 0.000783], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 64), 'BER'): [[0.09964880952380953, 0.04829047619047619, 0.016652380952380953, 0.003664285714285714, 0.0004321428571428571], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 8), 'BER'): [[0.12447448979591837, 0.0393, 0.006422448979591836, 0.000571938775510204, 3.724489795918367e-05], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 24), 'BLER'): [[0.479625, 0.21115, 0.053175, 0.007525, 0.00065], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 40), 'BLER'): [[0.556125, 0.2832, 0.094625, 0.01875, 0.00175], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 56), 'BLER'): [[0.569675, 0.323125, 0.127075, 0.031725, 0.005625], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 64), 'BLER'): [[0.461425, 0.237025, 0.084, 0.018825, 0.00225], [1.0, 2.0, 3.0, 4.0, 5.0]],
((0.5, 8), 'BLER'): [[0.45295, 0.163375, 0.03015, 0.0032, 0.000175], [1.0, 2.0, 3.0, 4.0, 5.0]]}


def update_CA_dicts(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        next(myfile, None)
        test_results = {}

        for row in myfile:
            if row[0] == '8' and row[7] == 'U=0, E:64, H-check, SCL':
                dict_getter = test_results.get((float(row[1]), float(row[8])), [0] * (4))
                dict_getter[0] += int(row[0]) * int(row[4])  # tot information bits sent
                dict_getter[1] += int(row[4])  # tot runs
                dict_getter[2] += int(row[5])  # tot bit errors
                dict_getter[3] += int(row[6])  # tot block errors
                test_results.update({(float(row[1]), float(row[8])): dict_getter})
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


def update_RM_dicts(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        next(myfile, None)
        test_results = {}

        for row in myfile:
            ue = row[7].replace(", ", "=").split("=")
            U, E = int(ue[1]), int(ue[3])

            dict_getter = test_results.get((float(row[1]), float(row[8]), U), [0]*(4))
            dict_getter[0] += int(row[0]) * int(row[4])  # tot information bits sent
            dict_getter[1] += int(row[4])  # tot runs
            dict_getter[2] += int(row[5])  # tot bit errors
            dict_getter[3] += int(row[6])  # tot block errors
            test_results.update({(float(row[1]), float(row[8]), U): dict_getter})
        for key, value in test_results.items():
            test_results.update({key: [value[0], value[1], value[2] / (value[0]), value[3] / value[1]]})
        ber_plot, bler_plot = {}, {}
        for key, value in test_results.items():
            ber_get = ber_plot.get((key[0], key[2]) , [[], []])
            bler_get = bler_plot.get((key[0], key[2]), [[], []])

            ber_get[0].append(value[2]);
            ber_get[1].append(key[1])
            bler_get[0].append(value[3]);
            bler_get[1].append(key[1])

            ber_plot.update({(key[0], key[2]): ber_get})
            bler_plot.update({(key[0], key[2]): bler_get})
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

#update_CA_dicts('Tests/CA_BSC.csv')

#update_RM_dicts('Tests/RM_BEC.csv')
Graph_Maker.plot_dict(CA_BEC, 'BEC_CA')