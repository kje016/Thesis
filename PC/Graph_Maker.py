import matplotlib.pyplot as plt
import csv
import os


def csv_to_plot(input_file):
    file = csv.reader(input_file)
    next(file, None)
    test_results, y_axis = {}, {}
    for row in csv.reader(input_file):
        dict_getter = test_results.get((float(row[1]),float(row[8]) ), [0]*(4))
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

    for key, value in ber_plot.items():
        plt.plot(value[1], value[0], label=str(key))

    name = input_file.name.split('/')[-1].replace('.csv', '')
    plt.xlabel('SNR')
    plt.ylabel('BER')
    plt.title(f'Bit Error Rate {name}')
    #plt.gca().invert_xaxis()
    plt.legend()    # show a legend on the plot
    plt.show() # function to show the plot

    for key, value in bler_plot.items():
        plt.plot(value[1], value[0], label=str(key))
    plt.xlabel('SNR')
    plt.ylabel('BLER')
    plt.title(f'Block Error Rate {name}')
    #plt.gca().invert_xaxis()
    plt.legend()    # show a legend on the plot
    plt.show() # function to show the plot



#csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/AWGN_SC.csv")))
#res = csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/BSC_SC.csv")))
res = csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/BEC_SC.csv")))