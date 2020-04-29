#
# Example of using larch to reproduce the second part of the
# XAFS processing tutorial 
#

# libraries/functions used in this script are obtained using import

# larch function for reading ascii data (usually txt)
from larch.io import read_ascii
# larch function for converting a group to a dictionary and
# a dictionary into a group
from larch.utils import group2dict, dict2group
# larch function for normalisation and flattening
from larch.xafs import pre_edge
# larch function for post-edge background substraction
from larch.xafs import autobk
# larch plot labels to save us from writing 
from larch.wxlib import plotlabels as plab
# larch function for fourier transform
from larch.xafs import xftf
# logarithm function from numpy
from numpy import log
# ploting library
import matplotlib.pyplot as plt

# these two lines define constants to help formating
# subscripta and supperscripts for printing
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

# custom function for calculating mu,
# could need indicating which arrays to use
# for calculation
def get_mu(xafs_group):
    # calculate mu
    mu_e = log(abs(xafs_group.i0/xafs_group.it))
    # get a dictionary from te group
    xafs_dict = group2dict(xafs_group)
    # add mu to the dictionary
    xafs_dict['mu'] = mu_e
    xafs_group = dict2group(xafs_dict)
    return xafs_group

# custom function for creating a copy of a
# group
def copy_group(xafs_group):
    # get a dictionary from te group
    xafs_dict = group2dict(xafs_group).copy()
    #xafs_dict['filename']=xafs_dict['filename']+"_copy"
    new_group = dict2group(xafs_dict)
    new_group.filename = new_group.filename+"_copy" 
    return new_group

# custom function for printing group properties
# needs fixing as groups can change depending
# on the source data
def print_contents(xafs_group):
  # the following prints show the contents of each
    print("path:\t\t", xafs_group.path)
    print("filename:\t", xafs_group.filename)
    print(xafs_group.header)
    print(xafs_group.data)
    print(xafs_group.attrs)
    print(xafs_group.energy)
    print(xafs_group.i0)
    print(xafs_group.it)
    print(xafs_group.inttime)
    print(xafs_group.array_labels)


 #######################################################
# |       Restore state of previous session           | #
# V                                                   V #
 #######################################################

# read the data    
fe_xafs = read_ascii('XAFSExamples/Fe_standards/Fe_lepidocrocite.000')

# calculate mu and normalise with background extraction
fe_xafs = get_mu(fe_xafs)

# calculate pre-edge and post edge and add them to group
pre_edge(fe_xafs, group=fe_xafs)

# perform background removal
autobk(fe_xafs) # using defaults so no additional parameters are passed

# calculate fourier transform
xftf(fe_xafs, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

# create copy of initial group
fe_xafs_copy = copy_group(fe_xafs)

# redo calculations with modified bacground parameter
autobk(fe_xafs_copy, rbkg=0.2)
xftf(fe_xafs_copy, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

# plot magnitudes in r-space to see we are at the same point where we left in
# session 1
plt.plot(fe_xafs_copy.r, fe_xafs_copy.chir_mag,label=fe_xafs_copy.filename)
plt.plot(fe_xafs.r, fe_xafs.chir_mag,label=fe_xafs.filename)
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
plt.show()

# Background tweaking
# Making the background parameter small caused some issues, but so does changing
# it to a larger value. The following example shows what happens with a larger
# rbkg parameter

# https://vimeo.com/340215763 00:50
# change rbkg to a larger number and redo calculations
# with modified bacground parameter
autobk(fe_xafs_copy, rbkg=2.0)
xftf(fe_xafs_copy, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

# plot magnitude in r-space to see results of change
plt.plot(fe_xafs_copy.r, fe_xafs_copy.chir_mag,label=fe_xafs_copy.filename)
plt.plot(fe_xafs.r, fe_xafs.chir_mag,label=fe_xafs.filename)
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
plt.show()

# When looking again in energy, the result is that the background and the signal
# are now closer than before changing the parameter
# https://vimeo.com/340215763 02:19
# visualise and interpret in energy
plt.plot(fe_xafs_copy.energy, fe_xafs_copy.bkg, label='background')
plt.plot(fe_xafs_copy.energy, fe_xafs_copy.mu, label=fe_xafs_copy.filename) # plot mu 
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlabel(plab.energy)
plt.ylabel(r'x$\mu$(E)')
plt.title(fe_xafs_copy.filename+" in Energy")
plt.legend()
plt.show()

# varying the rbkg parameter it is possible to find an interval of values in which
# the changes in the r-space are safe (don't produce spurious peaks).
# with larch the safe margin is from 0.7 to 1.0

# https://vimeo.com/340215763 02:50
# change rbkg to a larger number and redo calculations
# with modified bacground parameter
autobk(fe_xafs_copy, rbkg=0.7)
xftf(fe_xafs_copy, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

# plot magnitude in r-space to see results of change
plt.plot(fe_xafs_copy.r, fe_xafs_copy.chir_mag,label=fe_xafs_copy.filename)
plt.plot(fe_xafs.r, fe_xafs.chir_mag,label=fe_xafs.filename)
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
plt.show()


# Now see what happens after 1.0, example 1.2, in this case the peak afer 1
# becomes clearer, while the values before 1 are again closer to 0. 
# However, according to Ravel, this is not good because the peak seems too high.
# Notice that in the tutorial this effect is reached at 1.4 not 1.2

# https://vimeo.com/340215763 05:50
# change rbkg to a larger number and redo calculations
# with modified bacground parameter
autobk(fe_xafs_copy, rbkg=1.2)
xftf(fe_xafs_copy, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

# plot magnitude in r-space to see results of change
plt.plot(fe_xafs_copy.r, fe_xafs_copy.chir_mag,label=fe_xafs_copy.filename)
plt.plot(fe_xafs.r, fe_xafs.chir_mag,label=fe_xafs.filename)
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
