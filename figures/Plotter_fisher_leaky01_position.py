import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import os


filename01pos = "data/simulation2017-06-14-20.34.34/fisher_leaky_position"
filename01 = "data/simulation2017-06-14-20.34.34/fisher_leaky_momentum"
filename01_times = "data/simulation2017-06-14-20.34.34/times"
fisher01pos = np.load(filename01pos + '.npy')
fisher01 = np.load(filename01 + '.npy')
times01 = np.load(filename01_times + '.npy')


filename00 = "data/simulation2017-05-05-15.34.15/fisher_result"

filename05pos = "data/simulation2017-06-24-14.23.30/fisher_leaky_position"
filename05 = "data/simulation2017-06-24-14.23.30/fisher_leaky_momentum"
filename_times05 = "data/simulation2017-06-24-14.23.30/times"
fisher05 = np.load(filename05 + '.npy')
fisher05pos = np.load(filename05pos + '.npy')
times05 = np.load(filename_times05 + '.npy')

for n in range(0,2):
  times01 = times01[::2]
  fisher01 = fisher01[::2]
  fisher01pos = fisher01pos[::2]
  times05 = times05[::2]
  fisher05 = fisher05[::2]
  fisher05pos = fisher05pos[::2]


""" Latex Parameters """
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
params = {'font.size': 50,
              'axes.labelsize': 50,
              'legend.fontsize': 40,
              'xtick.labelsize': 50,
              'ytick.labelsize': 50,
              'text.latex.preamble': [r'\usepackage{physics}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)

""" Create plot """
plt.figure(figsize=(13,13))

""" Axes """
ax = plt.subplot(111)
ax.tick_params(axis='both', which='major')
plt.xticks([ 0, np.pi/2, np.pi, 3*np.pi/2,  2*np.pi], [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

""" Plotting """
plt.plot(times01, fisher01pos, '-o', color = 'k', markeredgewidth=0.0,  label = r'$\bar{\gamma}= 0.1$', markersize = 12, linewidth=2)

""" Limit data range """ 
plt.xlim([0, 2*np.pi])

""" Labels """
plt.xlabel(r'\textbf{time} (t)')
plt.ylabel(r'\textbf{CFI} ($\bar{I}_F$)')
plt.subplots_adjust(top=0.8)
plt.legend(loc = 2, frameon = False)
ax.annotate(r'$\lambda = 0$', xy=(1, 1), xytext=(-100, -36), fontsize=40,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0, boxstyle="round,pad=0.4"),
    horizontalalignment='right', verticalalignment='top')


plt.savefig("Fisher_leaky01_position.pdf",transparent=True, dpi=600,bbox_inches='tight')
plt.show()

