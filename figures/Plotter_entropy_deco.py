import matplotlib.pyplot as plt
import numpy as np
import os

filename = "data/EntropyDecoherence0.2"
with open(filename) as f:
     for line in f:
        readout = line.split(",")

entropies = [float(i) for i in readout]
del entropies[::2]
del entropies[150:]

times = np.linspace(0., 6.*np.pi, len(entropies))

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


""" Plotting """
plt.plot(times, entropies, '-o', color = 'k', markeredgewidth=0.0,  label = r'$\bar{\gamma} = 0.2$',  markersize = 12, linewidth=2)


""" Axes """
ax = plt.subplot(111)
ax.tick_params(axis='both', which='major')
plt.xticks([0.0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi, 5*np.pi, 6*np.pi], [r'$0$', r'$\pi$',r'$2\pi$',r'$3\pi$', r'$4\pi$', r'$5\pi$', r'$6\pi$'])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8], [r'$0.0$',  r'$0.2$', r'$0.4$', r'$0.6$', r'$ 0.8$'])


""" Limit data range """ 
plt.xlim([0, 6*np.pi])

""" Labels """
plt.xlabel(r'\textbf{time} $t$')
plt.ylabel(r'\textbf{entropy} '+ r'$S(t)$')
plt.subplots_adjust(top=0.8)
plt.legend(loc = 1,frameon=False)



plt.savefig("entropy_deco.pdf",transparent=True, dpi=600,bbox_inches='tight')
plt.show()