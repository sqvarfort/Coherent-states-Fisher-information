"""
Fisher information for position measurement

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

filename00 = "data/simulation2017-02-16-10.00.08/homodyne_fisher"
filename_times00 = "data/simulation2017-02-16-10.00.08/times"

filename00 = "data/simulation2017-04-25-16.09.58/homodyne_fisher"
filename_times00 = "data/simulation2017-04-25-16.09.58/times"

filename001 = "data/simulation2017-06-12-11.34.26/fisher_position"
filename_times001 = "data/simulation2017-06-12-11.34.26/times"

filename001 = "data/simulation2017-06-22-11.36.53/fisher_position"
filename_times001 = "data/simulation2017-06-22-11.36.53/times"

filename005 = "data/simulation2017-06-13-15.38.09/fisher_position"
filename_times005 = "data/simulation2017-06-13-15.38.09/times"

filename005 = "data/simulation2018-03-19-11.42.17/fisher_position"
filename_times005 = "data/simulation2018-03-19-11.42.17/times"

filename01 = "data/simulation2018-03-08-12.09.50/fisher_position"
filename_times01 = "data/simulation2018-03-08-12.09.50/times"


fisher00 = np.load(filename00 + '.npy')
times00 = np.load(filename_times00 + '.npy')
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

for i in range(0, len(fisher001)):
  if times001[i] < 2*np.pi:
    times001_short.append(times001[i])
    fisher001_short.append(fisher001[i])

for i in range(0, len(fisher005)):
  if times005[i] < 2*np.pi:
    times005_short.append(times005[i])
    fisher005_short.append(fisher005[i])

for i in range(0, len(times01)):
  if times01[i] < 2*np.pi:
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

times00 = times00[1::2]
fisher00 = fisher00[1::2]
times005_short = times005_short[1::2]
fisher005_short = fisher005_short[1::2]
times01_short = times01_short[1::2]
fisher01_short = fisher01_short[1::2]

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
plt.plot(times00, fisher00, '-o', color = 'k', markeredgewidth=0.0,  label = r'$\bar{\kappa} = 0$', markersize = 12, linewidth=2)
plt.plot(times001_short, fisher001_short, '-^', color = 'b', markeredgewidth=0.0, label = r'$\bar{\kappa} = 0.01$', markersize = 12, linewidth=2)
plt.plot(times005_short, fisher005_short, '-v', color = 'm', markeredgewidth=0.0, label = r'$\bar{\kappa} = 0.05$', markersize = 12, linewidth=2)
plt.plot(times01_short, fisher01_short, '-s', color = 'g', markeredgewidth=0.0, label = r'$\bar{\kappa} = 0.1$', markersize = 12, linewidth=2)

""" Limit data range """ 
plt.xlim([0, 2*np.pi])

""" Labels """
plt.xlabel(r'\textbf{time} (t)')
plt.ylabel(r'\textbf{Fisher information} ($\bar{I}_F$)')
plt.subplots_adjust(top=0.8)
plt.legend(loc = 2)
ax.annotate(r'$\lambda = 0$', xy=(1, 1), xytext=(-120, -33), fontsize=40,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.2, boxstyle="round,pad=0.3"),
    horizontalalignment='right', verticalalignment='top')


plt.savefig("Fisher_position.pdf",transparent=True, dpi=600,bbox_inches='tight')
plt.show()