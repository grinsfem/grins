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

data = loadtxt(datafile,comments="#",skiprows=3)

T = data[:,0]
omega_dot = data[:,1:n_species+1]

fig = plot.figure()
axes = fig.add_subplot(111)

axes.yaxis.major.formatter.set_powerlimits((0,0)) 

axes.set_xlabel(r"$T$ [K]")
# The "r" at the beginning is important. Not sure why yet.
axes.set_ylabel(r"$\dot{\omega}$ $[kg/m^3-s]$")

line_formats = [ "b-", "r-", "g-", "c-", "k-", \
                 "b--", "r--", "g--", "c--", "k--", \
                 "b:", "r:", "g:", "c:", "k:", \
                 "b-.", "r-.", "g-.", "c-.", "k-.", \
                 "b.", "r.", "g.", "c.", "k." ]

for i,s in enumerate(species):
    axes.plot(T, omega_dot[:,i], line_formats[i], label=r""+s)

axes.legend(loc="upper left")
axes.grid(True)

plot.savefig("omega_dot.png")

plot.show()

