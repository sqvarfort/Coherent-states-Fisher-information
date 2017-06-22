import matplotlib.pyplot as plt
from numpy import *
import os




filename1 = "data/GeneralQuantumFisher1photons"
filename2 = "data/GeneralQuantumFisher2photons"
filename3 = "data/GeneralQuantumFisher3photons"

fisher1 = []
fisher2 = []
fisher3 = []

with open(filename1) as f:
     for line in f:
        fisher1.append(line)

with open(filename2) as f:
     for line in f:
        fisher2.append(line)        

with open(filename3) as f:
     for line in f:
        fisher3.append(line)

times = linspace(0, 2*pi, len(fisher1))


#plt.figure(figsize=(13,11))
#plt.plot(times01, fisher01)
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
              'font.size': 50,
              'axes.labelsize': 50,
             # 'text.fontsize': 10,
              'legend.fontsize': 30,
              'xtick.labelsize': 50,
              'ytick.labelsize': 50,
              'text.usetex': True,    #benutze latex fuer schrift encoding -> automatisch selbe schriftart wie in latex
              'text.latex.unicode': True,
              'font.family': 'serif',
              'font.serif': 'cm',
              #'figure.figsize': fig_size,
              'text.latex.preamble': [r'\usepackage{physics}', r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    plt.xlabel('Time $t$')
    plt.ylabel('QFI '+ r'$H_Q(t)$' + ' ' + r'(s$^4$/m$^2$)')
    ax = plt.subplot(111)
    ax.tick_params(axis='both', which='major', pad=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    [i.set_linewidth(2) for i in ax.spines.itervalues()]

    #plt.title('Fisher Information vs. time for ' + r'$\bar{g} = $' + str(self.args['gbar']) + ', $k = $' + str(self.k) + ', $N = $' + str(self.N) + ', $h = $' + str(self.h), size = 20, y=1.02)
    #plt.gca().grid(True, linewidth = 2)

    #plt.plot(times[0], data[0], '-o', color = 'k', label = '$1$ photon')
    #plt.plot(times[1], data[1], '-o', color = 'b', label = '$4$ photons')
    #plt.plot(times[2], data[2], '-o', color = 'r', label = '$9$ photons')

    #plt.plot(times[0], data[0], color = 'b', label = 'Analytic')
    plt.plot(times[0], data[0], '-o', color = 'k', markeredgewidth=0.0,  label = r'$|\alpha|^2 = 1$', markersize = 2, linewidth=3)
    plt.plot(times[0], data[1], '-o', color = 'b', markeredgewidth=0.0, label = r'$|\alpha|^2= 4$', markersize = 2, linewidth=3)
    plt.plot(times[0], data[2], '-o', color = 'r', markeredgewidth=0.0, label = r'$|\alpha|^2 = 9$', markersize = 2, linewidth=3)

    #plt.xticks([ 0, pi/2, pi, 3*pi/2,  2*pi], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'], size = 40)
    plt.xticks([ 0, pi/2, pi, 3*pi/2,  2*pi], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

    #plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], [r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$', r'$2.5$', r'$3.0$', r'$3.5$'], size = 30)
    #plt.yticks([-100, 0.0, 1000, 2000, 3000, 4000, 5000, 6000], [r'$0.0$', r'$1000$', r'$2000$', r'$3000$', r'$4000$', r'$5000$', r'$6000$'], size = 40)
    #plt.yticks([-20, 0.0, 100, 200, 250], [r'', r'$0$', r'$100$', r'$200$', r'$ $'])
    #plt.yticks([0.0,  200,  400,  600,  800], [r'$0.0$', r'$200$',  r'$400$', r'$600$', r'$800$'], size = 30)
    plt.subplots_adjust(bottom=0.15,left = 0.15) 
    #plt.xlim([0, pi/2])
    #plt.ylim([0,300])
    plt.ylim(ymin=-500)

    plt.legend(loc = 2)
    path = os.path.join("fisher_quantum2")
    plt.savefig(path +  ".jpg",transparent=True, dpi=600)
    plt.show()

plot_fisher([times], [fisher1, fisher2, fisher3], [0.0], "fisher_quantum2")
