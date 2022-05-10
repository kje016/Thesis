import matplotlib.pyplot as plt
import csv


def csv_to_plot(input_file):
    test_results, y_axis = {}, {}
    for row in input_file:
        row_vals = row.split(',')[:-1]
        dict_getter = test_results.get((float(row_vals[1]),float(row_vals[8]) ), [0]*(4))
        dict_getter[0]+=int(row_vals[0]) # tot information bits sent
        dict_getter[1]+= int(row_vals[4])   # tot runs
        dict_getter[2]+= int(row_vals[5])    # tot bit errors
        dict_getter[3]+= int(row_vals[6])   # tot block errors
        test_results.update({(float(row_vals[1]), float(row_vals[8])): dict_getter})

    for key, value in test_results.items():
        test_results.update({key : [value[0], value[1], value[2]/(value[0]*value[1]), value[3]/value[1]]})

    ber_plot, bler_plot = {}, {}
    for key, value in test_results.items():
        ber_get = ber_plot.get(key[0], [[], []])
        bler_get = bler_plot.get(key[0], [[], []])

        ber_get[0].append(value[2]); ber_get[1].append(key[1])
        bler_get[0].append(value[3]); bler_get[1].append(key[1])

        ber_plot.update({key[0] : ber_get})
        bler_plot.update({key[0] : bler_get})

    for key, value in ber_plot.items():
        plt.plot(value[0], value[1], label=str(key))

    plt.xlabel('SNR')
    plt.ylabel('BER')
    plt.title('Here goes the title')
    plt.gca().invert_xaxis()
    plt.legend()    # show a legend on the plot
    plt.show() # function to show the plot

    for key, value in bler_plot.items():
        plt.plot(value[0], value[1], label=str(key))
    plt.xlabel('SNR')
    plt.ylabel('BLER')
    plt.title('Here goes the title')
    plt.legend()    # show a legend on the plot
    plt.show() # function to show the plot



res = csv_to_plot(open('test_csv.txt'))