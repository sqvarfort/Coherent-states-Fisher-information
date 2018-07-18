import matplotlib.pyplot as plt
from numpy import *
import os

filename = "data/x_expg12016-08-20-15.02.42"

x_exp = []
p_exp = []
with open(filename) as f:
     for line in f:
        x_exp.append(float(line))


filename = "data/p_expg12016-08-20-15.02.42"
with open(filename) as f:
    for line in f:
        p_exp.append(float(line))

x_exp = [float(i) for i in x_exp]/sqrt(2.)
p_exp = [float(i) for i in p_exp]/sqrt(2.)


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

""" Plotting """
line = plt.plot(x_exp, p_exp, linewidth = 3,  color = 'black', label = r'$\bar{g} = 1$')



""" Labels """
plt.xlabel( r'$\langle \hat{X}_{\mathrm{C}} \rangle$')
plt.ylabel(r'$\langle \hat{P}_{\mathrm{C}} \rangle$')
#plt.subplots_adjust(top=0.8)
plt.legend(loc = 1, frameon = False)
plt.xticks([-0.5, 0.0, 0.5, 1.0, 1.5], [r'-$0.5$', r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$', r'$2.0$'])
plt.yticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5], [r'-$1.5$', r'-$1.0$', r'-$0.5$', r'$0.0$', r'$0.5$', r'$1.0$', r'$1.5$'])

""" Limit data range """ 
plt.xlim([-0.7, 1.7])
plt.ylim([-1.7, 1.7])

""" Arrows """
i = 0
line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

i = 650
line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

i = 2750
line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

i = 3650
line[0].axes.annotate('',xytext = (x_exp[i],p_exp[i]),xy = (x_exp[i + 1],p_exp[i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)

i = 630
line[0].axes.annotate('',xytext = (x_exp[-i],p_exp[-i]),xy = (x_exp[-i + 1],p_exp[-i +1]),arrowprops=dict(arrowstyle="fancy",color = 'k'),size = 50)


""" Save image """
plt.savefig("quadg1.pdf",transparent=True, dpi=600,bbox_inches='tight')
plt.show()



