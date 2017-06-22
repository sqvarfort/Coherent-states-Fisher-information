import matplotlib.pyplot as plt
from numpy import *
import os

filename = "data/x_expg22016-08-21-00.17.42"
with open(filename) as f:
    x_exp = f.readlines()
filename = "data/p_expg22016-08-21-00.17.42"
with open(filename) as f:
    p_exp = f.readlines()
print "Data read"




x_exp = [float(i) for i in x_exp]
p_exp = [float(i) for i in p_exp]


def plot_quads(x_exp, p_exp, filename):
    """
    function: plot_fisher
    - data/array: Fisher information to be plotted vs time.
    Output:
    - Plot/file
    """
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
              'text.latex.preamble': [r'\usepackage{physics}', r'\usepackage{amsmath}', r'\usepackage[dvips]{graphicx}\usepackage{xfrac}']}
    plt.rcParams.update(params)
    plt.figure(figsize=(13,11))
    # Use Latex
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
#    plt.gca().grid(True, linewidth = 2)

    #plt.xlabel(r'$\langle \hat{x} \rangle$', size = 40)
    #plt.ylabel(r'$\langle \hat{p} \rangle$', size = 40)
    plt.xlabel(r'Position quadrature ($\langle \hat{x} \rangle$)', size = 40)
    plt.ylabel(r'Momentum quadrature ($\langle \hat{p} \rangle$)', size = 40)
    ax = plt.subplot(111)
    [i.set_linewidth(2) for i in ax.spines.itervalues()]
    ax.tick_params(axis='both', which='major', pad=5, width=2)
    plt.subplots_adjust(bottom=0.15, left = 0.15)
    #plt.title('Fisher Information vs. time for ' + r'$\bar{g} = $' + str(self.args['gbar']) + ', $k = $' + str(self.k) + ', $N = $' + str(self.N) + ', $h = $' + str(self.h), size = 20, y=1.02)
    line = plt.plot(x_exp, p_exp, linewidth = 3,  color = 'black', label = r'$\bar{g} = 2$')
    #plt.plot(times[1], data[1], color = 'b', label = 'Analytic')
    #plt.xticks([0.0, pi/2, pi, 3*pi/2,  2.*pi], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'], size = 40)
    plt.xticks([-0.7,-0.5,0.0,0.5,1.0,1.5,2.0, 2.3], [r' ', r'$-0.5$', r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'], size = 40)

    #plt.xticks([-0.7,-0.5,0.0,0.5,1.0,1.5,2.0, 2.3], [r' ', r'$-0.5$', r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'], size = 40)
    #plt.xticks([0.0, 2*pi, 4*pi, 6*pi,  8*pi], [r'$0$', r'$2\pi$', r'$4\pi$', r'$6\pi$', r'$8\pi$'], size = 40)
    plt.yticks([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0], [r' ', r'$-1.5$', r'$-1.0$', r'$-0.5$', r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$'], size = 40)
    #plt.xlim([0, 8*pi])
    #plt.xlim([0, 2*pi])
    # Adding arrows
    #add_arrow(line[0], size = 50, position = 2, direction = 'right', color = 'k')
    #add_arrow(line[0], size = 50, position = 1.5, direction = 'right', color = 'k')

    i = 0
    line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

    i = 650
    line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

    #i = 1750
    #line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

    #i = 4580
    #line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

    i = 2980
    line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

    i = 3400
    line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

    i = 630
    line[0].axes.annotate('',xytext = (x_exp[-i],p_exp[-i]),xy = (x_exp[-i + 1],p_exp[-i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)


    plt.legend(loc = 1, fontsize = 40)

    path = os.path.join(filename)
    plt.savefig(path +  ".jpg",transparent=False, dpi=600)
    plt.show()

plot_quads(x_exp, p_exp, "quadsg2")
