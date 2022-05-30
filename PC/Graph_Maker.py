import csv
import os.path
import matplotlib as mpl
mpl.use('template')
import matplotlib.pyplot as plt


def prep_axis(x_axis, y_axis, channel):
    if x_axis[0] < x_axis[1] and channel != 'AWGN':
        return x_axis[::-1], y_axis[::-1]
    elif channel == 'AWGN' and x_axis[0] > x_axis[1]:
        return x_axis[::-1], y_axis[::-1]
    return x_axis, y_axis


def csv_to_plot(input_file):
    with open(input_file, mode='r',  newline='') as file:
        myfile = csv.reader(file)
        next(myfile, None)
        test_results, y_axis = {}, {}
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

    name = input_file.split('\\')[-1].replace('.csv', '')
    fig, ax = plt.subplots()
    ax.set_facecolor('#F7F7F7')
    #breakpoint()
    for key, value in ber_plot.items():
        plot_x = sorted(value[1])
        plot_y = [value[0][value[1].index(a)] for a in plot_x]
        plot_x, plot_y = prep_axis(plot_x, plot_y, name.split('_')[0])
        ax.semilogy(plot_x, plot_y, label=str(key), marker='d', markerfacecolor='none')
        ax.set(title=str(key))



    ax.grid(True, linewidth=0.5)
    ax.set_xlabel('SNR')
    ax.set_ylabel('BER')
    ax.set_title(f'Bit Error Rate {name}')
    if name.split('_')[0] != 'AWGN':
        ax.invert_xaxis()
    ax.legend()    # show a legend on the plot
    #ax.show() # function to show the plot
    fig.savefig(f'{name}_BER.svg')
    #breakpoint()

    fig, ax = plt.subplots()
    ax.set_facecolor('#F7F7F7')
    #breakpoint()
    for key, value in bler_plot.items():
        plot_x = sorted(value[1])
        plot_y = [value[0][value[1].index(a)] for a in plot_x]
        #plot_x, plot_y = prep_axis(plot_x, plot_y, name.split('_')[0])
        ax.semilogy(plot_x, plot_y, label=str(key), marker='d', markerfacecolor='none')
        ax.set(title=str(key))

    ax.grid(True, linewidth=0.5)
    ax.set_xlabel('SNR')
    ax.set_ylabel('BER')
    if name.split('_')[0] != 'AWGN':
        ax.invert_xaxis()
    ax.set_title(f'Block Error Rate {name}')
    ax.legend()    # show a legend on the plot
    #ax.show() # function to show the plot
    fig.savefig(f'{name}_BLER.svg')
    #plt.legend()    # show a legend on the plot
    #plt.show() # function to show the plot



#csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/AWGN_SC.csv")))
#res = csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/BSC_SC.csv")))
#res = csv_to_plot( open(os.path.expanduser("~/PycharmProjects/Thesis/PC/Tests/BEC_SC.csv")))
csv_to_plot('C:\\Users\\Kristian\\Desktop\\Thesis\\PySageMath\\PC\\Tests\\BSC_SC.csv')
#csv_to_plot(open(os.path.join("C:/Users/Kristian/Desktop/Thesis/PySageMath/PC/Tests", "AWGN_SCL.csv")))