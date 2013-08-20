import matplotlib
from matplotlib import rc
rc('text',usetex=True)

#matplotlib.use("PDF")

import matplotlib.pyplot as plot
from numpy import loadtxt                 
import sys

tick_label_fontsize=14
axis_label_fontsize=14
matplotlib.rc('xtick', labelsize=tick_label_fontsize )
matplotlib.rc(('xtick.major','xtick.minor'),  pad=10)
matplotlib.rc('ytick', labelsize=tick_label_fontsize)
matplotlib.rc('text.latex', preamble=[r"\usepackage{amsmath}"])
matplotlib.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                  'monospace':['Computer Modern Typewriter']})

if( len(sys.argv) < 2 ):
    print "Error: Must supply kinetic rates data file name!"
    sys.exit(1)

datafile = sys.argv[1]

file = open(datafile, "r")
file.readline()
species = file.readline().split()
file.close()


n_species = len(species)

data = loadtxt( datafile, comments="#", skiprows=3)

T = data[:,0]

h_cea_data       = data[:,             1:  n_species+1]
h_stat_mech_data = data[:,   n_species+1:2*n_species+1]
e_tr             = data[:, 2*n_species+1:3*n_species+1]
e_vib            = data[:, 3*n_species+1:4*n_species+1]
e_el             = data[:, 4*n_species+1:5*n_species+1]
e_0              = data[:, 5*n_species+1:6*n_species+1]
cv_trans         = data[:, 6*n_species+1:7*n_species+1]
cv_rot           = data[:, 7*n_species+1:8*n_species+1]

fig = plot.figure()
axes = fig.add_subplot(111)

axes.yaxis.major.formatter.set_powerlimits((0,0)) 

axes.set_xlabel(r"$T$ [K]")
axes.set_ylabel(r"$h$ [J/kg]")

cea_plot = ["b-", "r-", "g-", "c-", "k-" ]
stat_mech_plot = ["b--", "r--", "g--", "c--", "k--" ]
e_vib_plot = ["b:", "r:", "g:", "c:", "k:" ]
e_el_plot = ["b-.", "r-.", "g-.", "c-.", "k-." ]
e_0_plot = ["b.", "r.", "g.", "c.", "k." ]

for i,s in enumerate(species):
    axes.plot(T, h_cea_data[:,i], cea_plot[i], label=r"CEA "+s)

for i,s in enumerate(species):
    axes.plot(T, h_stat_mech_data[:,i], stat_mech_plot[i], label=r"StatMech "+s)

axes.legend(loc="upper left")
axes.grid(True)

plot.savefig("h.png")

fig2 = plot.figure()
axes2 = fig2.add_subplot(111)

axes2.yaxis.major.formatter.set_powerlimits((0,0)) 

axes2.set_xlabel(r"$T$ [K]")
axes2.set_ylabel(r"$\frac{h_{CEA} - h_{SM}}{h_{SM}}$")

for i,s in enumerate(species):
    axes2.plot(T, (h_cea_data[:,i]-h_stat_mech_data[:,i])/h_stat_mech_data[:,i], cea_plot[i], label=r""+s)

axes2.legend(loc="lower right")
axes2.grid(True)

plot.savefig("h_diff.png")

fig3 = plot.figure()
axes3 = fig3.add_subplot(111)

axes3.yaxis.major.formatter.set_powerlimits((0,0)) 

axes3.set_xlabel(r"$T$ [K]")
axes3.set_ylabel(r"$e$ [J/kg]")

for i,s in enumerate(species):
    axes3.plot(T, cv_trans[:,i]*T, cea_plot[i], label=r"$e_{trans}$ "+s)
    axes3.plot(T, cv_rot[:,i]*T, stat_mech_plot[i], label=r"$e_{rot}$ "+s)
    axes3.plot(T, e_vib[:,i], e_vib_plot[i], label=r"$e_{vib}$ "+s)
    axes3.plot(T, e_el[:,i], e_el_plot[i], label=r"$e_{el}$ "+s)
    axes3.plot(T, e_0[:,i], e_0_plot[i], label=r"$e_0$ "+s)

axes3.legend(loc="upper left")
axes3.grid(True)

plot.savefig("e.png")

plot.show()



