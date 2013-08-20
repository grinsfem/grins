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
    print "Error: Must supply transport data file name!"
    sys.exit(1)

datafile = sys.argv[1]

file = open(datafile, "r")
file.readline()
species = file.readline().split()
file.close()


n_species = len(species)

data = loadtxt( datafile, comments="#", skiprows=3)

T = data[:,0]

mu       = data[:,1]
k = data[:,2]
D = data[:,3:3+n_species+1]

line_formats = [ "b-", "r-", "g-", "c-", "k-", \
                 "b--", "r--", "g--", "c--", "k--", \
                 "b:", "r:", "g:", "c:", "k:", \
                 "b-.", "r-.", "g-.", "c-.", "k-.", \
                 "b.", "r.", "g.", "c.", "k." ]

fig = plot.figure()
axes = fig.add_subplot(111)

axes.yaxis.major.formatter.set_powerlimits((0,0)) 

axes.set_xlabel(r"$T$ [K]")
axes.set_ylabel(r"$mu$")

axes.plot(T,mu,"b-")

axes.grid(True)

plot.savefig("mu.png")


fig2 = plot.figure()
axes2 = fig2.add_subplot(111)

axes2.yaxis.major.formatter.set_powerlimits((0,0)) 

axes2.set_xlabel(r"$T$ [K]")
axes2.set_ylabel(r"$k$")

axes2.plot(T,k,"b-")

axes2.grid(True)

plot.savefig("k.png")



fig3 = plot.figure()
axes3 = fig3.add_subplot(111)

axes3.yaxis.major.formatter.set_powerlimits((0,0)) 

axes3.set_xlabel(r"$T$ [K]")
axes3.set_ylabel(r"$D$")

for i,s in enumerate(species):
    axes3.plot(T, D[:,i], line_formats[i], label=r"D "+s)

axes3.legend(loc="upper left")
axes3.grid(True)

plot.savefig("D.png")

plot.show()



