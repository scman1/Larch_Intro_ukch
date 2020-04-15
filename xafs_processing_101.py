#
# Example of using larch to reproduce the first part of the
# XAFS tutorial
#

# libraries/functions used in this script are obtained using import

# larch function for reading ascii data (usually txt)
from larch.io import read_ascii
# larch function for converting a group to a dictionary and
# a dictionary into a group
from larch.utils import group2dict, dict2group
# larch function for  normalisation and flattening
from larch.xafs import pre_edge
# larch function for post-edge background substraction
from larch.xafs import autobk

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


# https://vimeo.com/340207346 07:01
fe_xafs = read_ascii('XAFSExamples/Fe_standards/Fe_lepidocrocite.000')
# using vars(fe_data) we see that the object has the following properties: 
# path, filename, header, data, attrs, energy, i0, it, ir, inttime, and labels
# vars(fe_xafs)

# the following prints show the contents of each
print_contents(fe_xafs)

# https://vimeo.com/340207346 15:06
fe_xafs = get_mu(fe_xafs)

plt.plot(fe_xafs.energy, fe_xafs.mu, label='$\mu$')
plt.xlabel('Energy')
plt.ylabel('$\mu$')
plt.legend() # include the leyend in the plot
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.title(fe_xafs.filename+" $\mu$")
plt.show()

# https://vimeo.com/340207346 25:00
# calculate pre-edge and post edge and add them to group
pre_edge(fe_xafs, group=fe_xafs)

#show pre-edge and post-edge fitting to mu
plt.plot(fe_xafs.energy, fe_xafs.pre_edge, 'g', label='pre-edge') # plot pre-edge in green
plt.plot(fe_xafs.energy, fe_xafs.post_edge, 'r', label='post-edge')# plot post-edge in green
plt.plot(fe_xafs.energy, fe_xafs.mu, 'b', label=fe_xafs.filename) # plot mu in blue
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel('x$\mu$(E)') # label y axis
plt.legend() # show legend
plt.title("pre-edge and post_edge fitting to $\mu$")
plt.show()

# https://vimeo.com/340207346 27:00
plt.plot(fe_xafs.energy, fe_xafs.flat, label=fe_xafs.filename) # plot flattened and normalised energy
plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
plt.xlabel('Energy (eV)') # label y graph
plt.ylabel(r'normalised x$\mu$(E)') # label y axis
plt.title("normalised and flattened $\mu$")
plt.legend() # show legend
plt.show()

print( 'Element Symbol', "\t", fe_xafs.atsym)
print ('Edge',"\t\t", fe_xafs.edge)
print ('E0'.translate(SUB),"\t\t", fe_xafs.e0)
print ('Edge Step',"\t", fe_xafs.edge_step)

# https://vimeo.com/340207346 30:00
autobk(fe_xafs) # using defaults so no additional parameters are passed
plt.plot(fe_xafs.energy, fe_xafs.bkg, label='background')
plt.plot(fe_xafs.energy, fe_xafs.mu, label=fe_xafs.filename) # plot mu 
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.xlabel(r'Energy (eV)')
plt.ylabel(r'x$\mu$(E)')
plt.title(fe_xafs.filename+" in Energy")
plt.legend()
plt.show()

plt.plot(fe_xafs.k, fe_xafs.chi)
plt.xlabel(r'$k\, ({\rm\AA})^{-1}$')
plt.title(fe_xafs.filename+" K-$\chi$ ")
plt.grid(linestyle=':', linewidth=1) #show and format grid
plt.show()
