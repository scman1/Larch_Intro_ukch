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
# import the larch.io function for merging groups interpolating if necessary
from larch.io import merge_groups


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

# https://vimeo.com/340215763 35:00
plt.plot(gr_0.energy, gr_0.flat, label=gr_0.label) # plot flattened and normalised energy
plt.plot(gr_1.energy, gr_1.flat, label=gr_1.label) # plot flattened and normalised energy
plt.plot(gr_2.energy, gr_2.flat, label=gr_2.label) # plot flattened and normalised energy
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel(r'normalised x$\mu$(E)') # label y axis
plt.title('lepidocrocite')
plt.legend() # show legend
plt.show()

# https://vimeo.com/340215763 35:00
plt.plot(gr_0.energy, gr_0.flat, 'b', label=gr_0.label) # plot flattened and normalised energy
plt.plot(gr_1.energy, gr_1.flat,'g', label=gr_1.label) # plot flattened and normalised energy
plt.plot(gr_2.energy, gr_2.flat, 'y', label=gr_2.label) # plot flattened and normalised energy
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel(r'normalised x$\mu$(E)') # label y axis
plt.title('lepidocrocite')
plt.xlim(7100,7220)
plt.legend() # show legend
plt.show()

# https://vimeo.com/340215763 36:40
plt.plot(gr_0.k, gr_0.chi*gr_0.k**2, label=gr_0.label)
plt.plot(gr_1.k, gr_1.chi*gr_1.k**2, label=gr_1.label)
plt.plot(gr_2.k, gr_2.chi*gr_2.k**2, label=gr_2.label)
plt.xlabel(plab.r)
plt.ylabel(plab.chir.format(2))
plt.title("lepidocrocite in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,14.5)
plt.legend()
plt.show()

# https://vimeo.com/340215763 36:46
# plot magnitude in r-space
plt.title("lepidocrocite in R space")
plt.plot(gr_0.r, gr_0.chir_mag,label=gr_0.label)
plt.plot(gr_1.r, gr_1.chir_mag,label=gr_1.label)
plt.plot(gr_2.r, gr_2.chir_mag,label=gr_2.label)
plt.xlabel(plab.r)
plt.ylabel(plab.chirmag.format(3))
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.legend()
plt.xlim(0,6)
plt.show()

# E0
# The value for E0 is calculated during normalisation, but it can be set to an
# specific value, but his will alter the results of some of the calculations for
# instance changing the value for E0 in the second group (Fe_lepidocrocite_001).
# The change is almost imperceptible on the normalised graph.

# https://vimeo.com/340215763 39:20

# change e0, recalculate and plot
autobk(gr_1, e0 = 7124.74)
xftf(gr_1, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

plt.plot(gr_0.energy, gr_0.flat, 'y', label=gr_0.label) # plot flattened and normalised energy
plt.plot(gr_1.energy, gr_1.flat,'r', label=gr_1.label) # plot flattened and normalised energy
plt.plot(gr_2.energy, gr_2.flat, 'c', label=gr_2.label) # plot flattened and normalised energy
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel(r'normalised x$\mu$(E)') # label y axis
plt.title('lepidocrocite')
plt.xlim(7120,7140)
plt.legend() # show legend
plt.show()


# However, the plot of the k space shows several differences for the second
# group when compared to the other two. This causes a shift in the lower
# energies which decreases on higher energies.

# https://vimeo.com/340215763 36:40
plt.plot(gr_0.k, gr_0.chi*gr_0.k**2, label=gr_0.label)
plt.plot(gr_1.k, gr_1.chi*gr_1.k**2, 'b', label=gr_1.label)
plt.plot(gr_2.k, gr_2.chi*gr_2.k**2, label=gr_2.label)
plt.xlabel(plab.r)
plt.ylabel(plab.chir.format(2))
plt.title("lepidocrocite in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,14.5)
plt.legend()
plt.show()

# Alignment 
# Pending to find how alignment of various readings is performed in Larch

# Setting values for a group

# Using the group_names property of the project and a loop allows setting values
# accross groups easily, for instance to set E0 to the value of the first group
# can be done as follows.

# https://vimeo.com/340215763 47:03
common_e0 = 0.0
for group_x in fe_project._athena_groups:
    if common_e0 == 0.0:
        common_e0 = fe_project._athena_groups[group_x].e0
    else:
        fe_project._athena_groups[group_x].e0 = common_e0
    autobk(fe_project._athena_groups[group_x], e0 = common_e0)
    xftf(fe_project._athena_groups[group_x], kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')

# Plot again to show that the change has removed the shift in the k space

# https://vimeo.com/340215763 47:25
plt.plot(gr_0.k, gr_0.chi*gr_0.k**2, label=gr_0.label)
plt.plot(gr_1.k, gr_1.chi*gr_1.k**2, 'b', label=gr_1.label)
plt.plot(gr_2.k, gr_2.chi*gr_2.k**2, label=gr_2.label)
plt.xlabel(plab.r)
plt.ylabel(plab.chir.format(2))
plt.title("lepidocrocite in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,15)
plt.legend()
plt.show()

# Merge groups 
# Once all signals look similar to each other as in the last example above, the
# next step in a typical analysis is to merge the signals in an average and
# continue working on the average.

# https://vimeo.com/340215763 54:00
# use the larch.io function for merging groups 

fe_merge = merge_groups(list(fe_project._athena_groups.values()))
autobk(fe_merge)
xftf(fe_merge, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')


# The result is a new group wich can be plotted over all the other groups.
# In the following image, the components are all in blue while the merge is
# shown on top as a yellow line.

plt.plot(gr_0.k, gr_0.chi*gr_0.k**2, 'b', label=gr_0.label)
plt.plot(gr_1.k, gr_1.chi*gr_1.k**2, 'b', label=gr_1.label)
plt.plot(gr_2.k, gr_2.chi*gr_2.k**2, 'b', label=gr_2.label)
plt.plot(fe_merge.k, fe_merge.chi*fe_merge.k**2, 'y', label="Merged")
plt.xlabel(plab.r)
plt.ylabel(plab.chir.format(2))
plt.title("lepidocrocite in k space")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlim(0,15)
plt.legend()
plt.show()

# https://vimeo.com/340215763 55:00
plt.plot(gr_0.energy, gr_0.flat, 'b', label=gr_0.label) # plot flattened and normalised energy
plt.plot(gr_1.energy, gr_1.flat, 'b', label=gr_1.label) # plot flattened and normalised energy
plt.plot(gr_2.energy, gr_2.flat, 'b', label=gr_2.label) # plot flattened and normalised energy
plt.plot(fe_merge.energy, fe_merge.flat, 'y', label="Merged")
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel(r'normalised x$\mu$(E)') # label y axis
plt.title('lepidocrocite')
plt.xlim(7100,7200)
plt.legend() # show legend
plt.show()

# name the group
fe_merge.label = "Lepidocrocite"
# add to the project
fe_project2 = create_athena(project_name)
for a_group in list(fe_project._athena_groups.values()):
    fe_project2.add_group(a_group)
# save the project
fe_project2.add_group(fe_merge)
fe_project2.save()

# https://vimeo.com/340215763 55:00
plt.plot(gr_0.energy, gr_0.flat, 'b', label=gr_0.label) # plot flattened and normalised energy
plt.plot(gr_1.energy, gr_1.flat, 'b', label=gr_1.label) # plot flattened and normalised energy
plt.plot(gr_2.energy, gr_2.flat, 'b', label=gr_2.label) # plot flattened and normalised energy
plt.plot(fe_merge.energy, fe_merge.flat, 'y', label=fe_merge.label)
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel(r'normalised x$\mu$(E)') # label y axis
plt.title('lepidocrocite')
plt.xlim(7100,7200)
plt.legend() # show legend
plt.show()

