import matplotlib.pyplot as plt
from numpy import *
import os





filename00 = "data/simulation2017-05-05-15.34.15/fisher_result"
filename00 = "data/simulation2017-05-05-14.03.12/fisher_result"

filename00 = "data/simulation2017-06-14-20.34.34/fisher_leaky_momentum"
filename00_times = "data/simulation2017-06-14-20.34.34/times"
fisher = load(filename00 + '.npy')
times = load(filename00_times + '.npy')



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
    plt.ylabel('Leaking photons '+ r'$\bar{I}_F(t)$')
    ax = plt.subplot(111)
    ax.tick_params(axis='both', which='major', pad=10)
    [i.set_linewidth(2) for i in ax.spines.itervalues()]

    #plt.title('Fisher Information vs. time for ' + r'$\bar{g} = $' + str(self.args['gbar']) + ', $k = $' + str(self.k) + ', $N = $' + str(self.N) + ', $h = $' + str(self.h), size = 20, y=1.02)
    #plt.gca().grid(True, linewidth = 2)

    #plt.plot(times[0], data[0], '-o', color = 'k', label = '$1$ photon')
    #plt.plot(times[1], data[1], '-o', color = 'b', label = '$4$ photons')
    #plt.plot(times[2], data[2], '-o', color = 'r', label = '$9$ photons')

    #plt.plot(times[0], data[0], color = 'b', label = 'Analytic')
    plt.plot(times[0], data[0], '-o', color = 'k', markeredgewidth=0.0,  label = '$\gamma = 0.5$', markersize = 2, linewidth=3)
    #plt.plot(times[1], data[1], '-o', color = 'b', markeredgewidth=0.0, label = '$\kappa = 0.05$', markersize = 2)
    #plt.plot(times[2], data[2], '-o', color = 'r', markeredgewidth=0.0, label = '$\kappa = 0.1$', markersize = 2)

    #plt.xticks([ 0, pi/2, pi, 3*pi/2,  2*pi], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'], size = 40)
    plt.xticks([ 0, pi/2, pi, 3*pi/2,  2*pi], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

    plt.yticks([-1, 0.0, 2, 4, 6, 8, 10, 12], [' ', r'$0$', r'$2$', r'$4$', r'$6$', r'$8$', r'$10$', r'$12$'])
    #plt.yticks([0.0, 1000, 2000, 3000, 4000, 5000, 6000], [r'$0.0$', r'$1000$', r'$2000$', r'$3000$', r'$4000$', r'$5000$', r'$6000$'], size = 40)
    #plt.yticks([-20, 0.0, 100, 200, 250], [r'', r'$0$', r'$100$', r'$200$', r'$ $'])
    #plt.yticks([0.0,  200,  400,  600,  800], [r'$0.0$', r'$200$',  r'$400$', r'$600$', r'$800$'], size = 30)
    plt.subplots_adjust(bottom=0.15,left = 0.15) 
    #plt.xlim([0, pi/2])
    #plt.ylim([0,300])

    plt.legend(loc = 2)
    path = os.path.join("Fisher_leaky")
    plt.savefig(path +  ".jpg",transparent=True, dpi=600)
    plt.show()

plot_fisher([times], [fisher], [0.0], "fisher_leaky")
