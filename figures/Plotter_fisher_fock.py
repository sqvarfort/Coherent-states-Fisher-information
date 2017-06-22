import matplotlib.pyplot as plt
from numpy import *
import os

"""








filename1 = "data_processing/" + "GeneralQuantumFisher1photons"
filename2 = "data_processing/" + "GeneralQuantumFisher2photons"
filename3 = "data_processing/" + "GeneralQuantumFisher3photons"

with open(filename1) as f:
     fisher1 = f.readlines()
with open(filename2) as f:
     fisher2 = f.readlines()
with open(filename3) as f:
     fisher3 = f.readlines()

fisher1 = [float(i) for i in fisher1]
fisher2 = [float(i) for i in fisher2]
fisher3 = [float(i) for i in fisher3]

times = linspace(0, 2*pi, len(fisher1))

filename_pure = "data_processing/" + "FisherInfoN20Homodyne"
filename_deco01 = "data_processing/" + "FisherHomodyneN25Deco0.1"
filename_deco05 = "data_processing/" + "FisherHomodyneN25Deco0.5"
with open(filename_pure) as f:
     for line in f:
        fisher_pure = line.split(",")
with open(filename_deco01) as f:
     for line in f:
        fisher_deco01 = line.split(",")
with open(filename_deco05) as f:
     for line in f:
        fisher_deco05 = line.split(",")


fisher_pure = [float(i) for i in fisher_pure]
fisher_deco01 = [float(i) for i in fisher_deco01]
fisher_deco05 = [float(i) for i in fisher_deco05]


times_pure = list(load("data_processing/" + "FisherInfoN20HomodyneTimes.npy"))
times_deco01 = list(load("data_processing/" + "times2016-08-21-08.05.18.npy"))
times_deco05 = list(load("data_processing/" + "times2016-08-20-01.16.02.npy"))

for i in range(0,1):
    del fisher_pure[::2]
    del fisher_deco01[::2]
    del fisher_deco05[::2]
    del times_pure[::2]
    del times_deco01[::2]
    del times_deco05[::2]


filename_pure = "data_processing/" + "fisher_mirror_N30"
filename_deco = "data_processing/" + "FisherN25MirrorDeco0.1"

with open(filename_pure) as f:
     for line in f:
        fisher_pure = line.split(",")

with open(filename_deco) as f:
     for line in f:
        fisher_deco = line.split(",")

fisher_pure = [float(i) for i in fisher_pure]
fisher_deco = [float(i) for i in fisher_deco]
del fisher_pure[::2]
del fisher_pure[::2]
del fisher_pure[::2]


times_pure = linspace(0, 2*pi, len(fisher_pure))
times_deco = list(load("data_processing/" + "times2016-08-21-08.05.18.npy"))

del fisher_deco[::2]
del fisher_deco[::2]
del fisher_deco[::2]

del times_deco[::2]
del times_deco[::2]
del times_deco[::2]

"""

filename00 = "data/simulation2017-02-07-13.45.35/fock_fisher"
filename_times00 = "data/simulation2017-02-07-13.45.35/times"

filename005 = "data/simulation2017-02-08-09.36.08/fock_fisher"
filename_times005 = "data/simulation2017-02-08-09.36.08/times"



fisher00 = load(filename00 + '.npy')
times00 = load(filename_times00 + '.npy')

fisher005 = load(filename005 + '.npy')
times005 = load(filename_times005 + '.npy')


# Mask the decoherence
fisher005 = ma.masked_where(fisher005 > 3.8, fisher005)
fisher005 = ma.masked_where(fisher005 < 0.0, fisher005)

#plt.show()

def plot_fisher(times, data, chirange, filename):
    """
    function: plot_fisher
    - data/array: Fisher information to be plotted vs time.
    Output:
    - Plot/file
    """
    plt.figure(figsize=(13,11))
    # Use Latex
    params = {'backend': 'ps',
              'font.size': 12,
              'axes.labelsize': 12,
             # 'text.fontsize': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,    #benutze latex fuer schrift encoding -> automatisch selbe schriftart wie in latex
              'text.latex.unicode': True,
              'font.family': 'serif',
              'font.serif': 'cm',
              #'figure.figsize': fig_size,
              'text.latex.preamble': [r'\usepackage{physics}', r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    plt.xlabel('$t$', size = 40)
    plt.ylabel('$I_F$', size = 40)
    ax = plt.subplot(111)
    ax.tick_params(axis='both', which='major', pad=10)
    #plt.title('Fisher Information vs. time for ' + r'$\bar{g} = $' + str(self.args['gbar']) + ', $k = $' + str(self.k) + ', $N = $' + str(self.N) + ', $h = $' + str(self.h), size = 20, y=1.02)
    #plt.gca().grid(True, linewidth = 2)

    #plt.plot(times[0], data[0], '-o', color = 'k', label = '$1$ photon')
    #plt.plot(times[1], data[1], '-o', color = 'b', label = '$4$ photons')
    #plt.plot(times[2], data[2], '-o', color = 'r', label = '$9$ photons')

    #plt.plot(times[0], data[0], color = 'b', label = 'Analytic')
    plt.plot(times[0], data[0], '-o', color = 'k', markeredgewidth=0.0, label = '$\kappa = 0.0$')
    plt.plot(times[1], data[1], '-o', color = 'b', markeredgewidth=0.0, label = '$\kappa = 0.05$')
    plt.xticks([ 0, pi/2, pi, 3*pi/2,  2*pi], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'], size = 40)

    plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], [r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$', r'$2.5$', r'$3.0$', r'$3.5$'], size = 30)
    #plt.yticks([0.0, 1000, 2000, 3000, 4000, 5000, 6000], [r'$0.0$', r'$1000$', r'$2000$', r'$3000$', r'$4000$', r'$5000$', r'$6000$'], size = 40)
    #plt.yticks([0.0, 100, 200, 300, 400, 500, 600, 700, 800, 900], [r'$0.0$', r'$100$', r'$200$', r'$300$', r'$400$', r'$500$', r'$600$', r'$700$', r'$800$', r'$900$'], size = 40)
    #plt.yticks([0.0,  200,  400,  600,  800], [r'$0.0$', r'$200$',  r'$400$', r'$600$', r'$800$'], size = 30)
    plt.subplots_adjust(bottom=0.15)
    #plt.xlim([0, pi/2])
    #plt.ylim([0,300])

    plt.legend(loc = 1, fontsize = 30)
    path = os.path.join(filename)
    plt.savefig(path +  ".pdf",transparent=True, dpi=600)
    plt.show()


plot_fisher([times00, times005], [fisher00, fisher005], [0.0], "Fisher_Fock")
