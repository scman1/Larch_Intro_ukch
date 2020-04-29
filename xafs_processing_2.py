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

# import the larch.io libraries for managing athena files
from larch.io import create_athena, read_athena, extract_athenagroup

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
# |     By default athena recalculates everything     | #
# |     so we can create a function that calculates   | #
# V    mu, pre_edge, autobk and xftf for a group      V #
 #######################################################

def calc_with_defaults(xafs_group):
    # calculate mu and normalise with background extraction
    # should let the user specify the colums for i0, it, mu, iR. 
    if not hasattr(xafs_group, 'mu'):
        xafs_group = get_mu(xafs_group)    
    # calculate pre-edge and post edge and add them to group
    pre_edge(xafs_group, group=fe_xafs)
    # perform background removal
    autobk(xafs_group) # using defaults so no additional parameters are passed
    # calculate fourier transform
    xftf(xafs_group, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
    return xafs_group

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
plt.show()

# k-weight
# All the parameters can also influence data analysis. For instance looking 
# again at the initial group, we can analise how kweight affects the analyis. First
# take a look at backgroun in energy. In this case, the values after 7200 are
# so close to the background that they seem to disappear. However, the data
# still contains useful values which can be enhanced to reveal meaningful
# information.

# https://vimeo.com/340215763 12:08
plt.plot(fe_xafs.energy, fe_xafs.bkg, label='background')
plt.plot(fe_xafs.energy, fe_xafs.mu, label=fe_xafs.filename) # plot mu 
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlabel(r'Energy (eV)')
plt.ylabel(r'x$\mu$(E)')
plt.title(fe_xafs.filename+" in Energy")
plt.legend()
plt.show()

# This can be observed by plotting in the k space whith kweight 2. This enhances
# the signals after 6, making them easier to see.

# https://vimeo.com/340215763 12:37
plt.plot(fe_xafs.k, fe_xafs.chi*fe_xafs.k**2, label=fe_xafs.filename)
plt.xlabel(r'Wavenumber ($A^{-1}$)')
plt.ylabel(plab.chikw.format(2))
plt.title(fe_xafs.filename+" in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,14.5)
plt.legend()
plt.show()

# This same signal before amplification can be seen in the following diagram. in
# this case, whitout amplification is difficult to see values above 6.

# https://vimeo.com/340215763 14:16
plt.plot(fe_xafs.k, fe_xafs.chi)
plt.ylabel(r'$\chi(k)$')
plt.xlabel(r'Wavenumber ($A^{-1}$)')
plt.title(fe_xafs.filename+" in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,14.5)
plt.show()

# Signal amplified by a k-weight of 1 showing more data

# https://vimeo.com/340215763 14:23
plt.plot(fe_xafs.k, fe_xafs.chi*fe_xafs.k, label=fe_xafs.filename)
plt.xlabel(r'Wavenumber ($A^{-1}$)')
plt.ylabel(plab.chikw.format(1))
plt.title(fe_xafs.filename+" in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,14.5)
plt.legend()
plt.show()

# Changing the kweight also affects the shape of the fourier transform. The
# following diagram shows the fourier transforms when k-weight is 1 and 2

xftf(fe_xafs, kweight=1.0, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
# https://vimeo.com/340215763 15:40
# plot magnitude in r-space
plt.plot(fe_xafs.r, fe_xafs.chir_mag,label=fe_xafs.filename+" kw=1")
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
xftf(fe_xafs, kweight=2.0, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
# https://vimeo.com/340215763 15:40
# plot magnitude in r-space
plt.plot(fe_xafs.r, fe_xafs.chir_mag,label=fe_xafs.filename+" kw=2")
plt.title(fe_xafs.filename+" at k = 1 and 2")
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
plt.show()

# The function of athena that allows visualising with different k-weights can
# be replicated as follows
plt.subplot(3, 1, 1)
plt.xlim(0,6)
xftf(fe_xafs, kweight=1.0, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
plt.plot(fe_xafs.r, fe_xafs.chir_mag, 'b', label=fe_xafs.filename+" kw=1")
plt.legend()
plt.subplot(3, 1, 2)
plt.xlim(0,6)
xftf(fe_xafs, kweight=2.0, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
plt.plot(fe_xafs.r, fe_xafs.chir_mag, 'r', label=fe_xafs.filename+" kw=2")
plt.legend()
plt.subplot(3, 1, 3)
plt.xlim(0,6)
xftf(fe_xafs, kweight=3.0, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
plt.plot(fe_xafs.r, fe_xafs.chir_mag, 'g', label=fe_xafs.filename+" kw=3")
plt.legend()
plt.show()

# Import more data
# Now import the next two files Fe_lepidocrocite.100 and Fe_lepidocrocite.200
# using the same defaults for calculating mu

fe_100 = read_ascii('XAFSExamples/Fe_standards/Fe_lepidocrocite.100')
fe_200 = read_ascii('XAFSExamples/Fe_standards/Fe_lepidocrocite.200')

# calculate mu and normalise with background extraction
# using defaults
fe_100 = calc_with_defaults(fe_100)
fe_200 = calc_with_defaults(fe_200)

# save as an athena project

project_name = 'XAFSExamples/Fe_standards/lepidocrocite.prj'
fe_project = create_athena(project_name)
fe_project.add_group(fe_xafs)
fe_project.add_group(fe_100)
fe_project.add_group(fe_200)
fe_project.save()
vars(fe_project)

# read athena project

fe_project = read_athena(project_name)
vars(fe_project)
gr_0 = extract_athenagroup(fe_project.Fe_lepidocrocite_000)
gr_0 = calc_with_defaults(gr_0)
#vars(gr_0)
gr_1 = fe_project.Fe_lepidocrocite_100
gr_1 = calc_with_defaults(gr_1)
gr_2 = fe_project.Fe_lepidocrocite_200
gr_2 = calc_with_defaults(gr_2)
plt.plot(gr_0.energy, gr_0.mu, label= gr_0.label + ' $\mu$')
plt.plot(gr_1.energy, gr_1.mu, label= gr_1.label + ' $\mu$')
plt.plot(gr_2.energy, gr_2.mu, label= gr_2.label + ' $\mu$')
plt.xlabel('Energy')
plt.ylabel('$\mu$')
plt.legend() # include the leyend in the plot
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.show()
