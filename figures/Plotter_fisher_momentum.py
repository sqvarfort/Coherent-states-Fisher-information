"""
Demo of TeX rendering.

You can use TeX to render all of your matplotlib text if the rc
parameter text.usetex is set.  This works currently on the agg and ps
backends, and requires that you have tex and the other dependencies
described at http://matplotlib.org/users/usetex.html
properly installed on your system.  The first time you run a script
you will see a lot of output from tex and associated tools.  The next
time, the run may be silent, as a lot of the information is cached in
~/.tex.cache

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

filename00 = "data/fisher_gbar1"

filename001 = "data/simulation2017-06-12-11.34.26/fisher_momentum"
filename_times001 = "data/simulation2017-06-12-11.34.26/times"

filename001 = "data/simulation2017-06-22-11.36.53/fisher_momentum"
filename_times001 = "data/simulation2017-06-22-11.36.53/times"

filename005 = "data/simulation2017-06-13-15.38.09/fisher_momentum"
filename_times005 = "data/simulation2017-06-13-15.38.09/times"

filename01 = "data/simulation2018-03-07-13.58.46/fisher_momentum"
filename_times01 = "data/simulation2018-03-07-13.58.46/times"

fisher00 = []
with open(filename00) as f:
     for line in f:
        fisher00.append(float(line))

times00 = np.linspace(np.pi, 3.*np.pi, len(fisher00))

fisher001 = np.load(filename001 + '.npy')
times001 = np.load(filename_times001 + '.npy')
fisher005 = np.load(filename005 + '.npy')
times005 = np.load(filename_times005 + '.npy')
times01 = np.load(filename_times01 + '.npy')
fisher01 = np.load(filename01 + '.npy')

times001_short = []
fisher001_short = []
fisher005_short=[]
times005_short = []
times01_short = []
fisher01_short = []

# print fisher01

# Remove all times that are less than pi
for i in range(0, len(times001)):
  if times001[i] > np.pi:
    times001_short.append(times001[i])
    fisher001_short.append(fisher001[i])

for i in range(0, len(times005)):
  if times005[i] > np.pi:
    times005_short.append(times005[i])
    fisher005_short.append(fisher005[i])

for i in range(0, len(times01)):
  if times01[i] > np.pi:
    times01_short.append(times01[i])
    fisher01_short.append(fisher01[i])

for i in range(0,2):
  times00 = times00[1::2]
  fisher00 = fisher00[1::2]
  times001_short = times001_short[1::2]
  fisher001_short = fisher001_short[1::2]
  times005_short = times005_short[1::2]
  fisher005_short = fisher005_short[1::2]
  times01_short = times01_short[1::2]
  fisher01_short = fisher01_short[1::2]

times01_short = times01_short[1::2]
fisher01_short = fisher01_short[1::2]
times01_short = times01_short[1::2]
fisher01_short = fisher01_short[1::2]

times00 = np.asarray(times00)
fisher00 = np.asarray(fisher00)
fisher00 = np.squeeze(fisher00)

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
plt.xticks([np.pi-0.5, np.pi, 3*np.pi/2,  2*np.pi, 5*np.pi/2, 3*np.pi], [r' ', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$5\pi/2$', r'$3\pi$'])

""" Plotting """
plt.plot(times00, fisher00, '-o', color = 'k', markeredgewidth=0.0,  label = r'$\bar{\kappa} = 0$', markersize = 12, linewidth=2)
plt.plot(times001_short, fisher001_short, '-^', color = 'b', markeredgewidth=0.0, label = r'$\bar{\kappa} = 0.01$', markersize = 12, linewidth=2)
plt.plot(times005_short, fisher005_short, '-v', color = 'm', markeredgewidth=0.0, label = r'$\bar{\kappa} = 0.05$', markersize = 12, linewidth=2)
plt.plot(times01_short, fisher01_short, '-s', color = 'g', markeredgewidth=0.0, label = r'$\bar{\kappa} = 0.1$', markersize = 12, linewidth=2)

""" Limit data range """ 
plt.xlim([np.pi, 3*np.pi])

""" Labels """
plt.xlabel(r'\textbf{time} (t)')
plt.ylabel(r'\textbf{CFI} ($\bar{I}_F$)')
plt.subplots_adjust(top=0.8)
plt.legend(loc = 2,frameon=False)
ax.annotate(r'$\lambda = 1/2$', xy=(1, 1), xytext=(-33, -33), fontsize=40,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0, boxstyle="round,pad=0.3"),
    horizontalalignment='right', verticalalignment='top')

plt.savefig("Fisher_momentum.pdf",transparent=True, dpi=600,bbox_inches='tight')
plt.show()