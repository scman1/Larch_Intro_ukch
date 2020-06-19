# required larch modules
# calculate pre-edge and post edge for normalisation
from larch.xafs import pre_edge
# perform background removal
from larch.xafs import autobk
# calculate fourier transform
from larch.xafs import xftf
# managing athena files
from larch.io import create_athena, read_athena, extract_athenagroup
# linear combination fitting
from larch.math import lincombo_fit

# additional libraries
from numpy import log
import matplotlib.pyplot as plt

 #######################################################
# |             Athena recalculates everything        | #
# |     so we can create a function that calculates   | #
# V               all for each new group              V #
 #######################################################

def calc_with_defaults(xafs_group):
    # calculate mu and normalise with background extraction
    # should let the user specify the colums for i0, it, mu, iR. 
    if not hasattr(xafs_group, 'mu'):
        xafs_group = get_mu(xafs_group)    
    # calculate pre-edge and post edge and add them to group
    pre_edge(xafs_group)
    # perform background removal
    autobk(xafs_group) # using defaults so no additional parameters are passed
    # calculate fourier transform
    xftf(xafs_group, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
    return xafs_group

# plot mu vs flat normalised mu for selected groups
def plot_NxmuE_E_athena_prj(athena_project, group_keys, group_names, title = "Normalised Mu vs E"):    
    # plot mu vs flat normalised mu for selected groups
    for group_key in group_keys:
        gr_0 = extract_athenagroup(athena_project._athena_groups[group_key])
        # recalculate normalisation
        calc_with_defaults(gr_0)
        plt.plot(gr_0.energy, gr_0.flat, label=group_names[group_key])

    # set plot format
    plt.xlabel("Energy")
    plt.ylabel("normalised xmuE" )
    plt.title(title)
    plt.grid(linestyle=':', linewidth=1) #show and format grid    
    plt.xlim(11860,12000)
    plt.legend()
    return plt


# https://vimeo.com/340216087 08:29 open cyanobacteria project 
project_name = 'XAFSExamples/Au+Cyanobacteria/cyanobacteria.prj'
cianobacteria_project = read_athena(project_name)
#vars(cianobacteria_project)

# https://vimeo.com/340216087 10:40 plot readings of sample 
# set group names
group_names = {'ozun':"0.12",'wtnk':"2.42", 'd_4_73':"4.73",
               'd_7_03':"7.03",'d_9_33':"9.33", 'qvxh':"20",'rjbc':"33",
               'lshy':"720",'hqlr':"Au Foil",'Au1_Cl':"Au1 Cl",
               'tscd':"Au3 Cl aq",'ixde':"Au hydroxide",'gcpx' :"Au cyanide",
               'nyux' :"Au thiocyanide",'qhxp':"Au sulphide",
               'Au_thiosulphate_aq':"Au thiosulphate aq",
               'ryzf':"Au thiomalate aq"}

# get the group keys for first 8 groups
group_keys = list(cianobacteria_project._athena_groups.keys())[0:8]

# plot mu vs flat normalised mu for selected groups
plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names, title = "Cianobacteria Readings")
plt.show()

# https://vimeo.com/340216087 12:10 plot standards 
# get the group keys for first last 9 groups
group_keys = list(cianobacteria_project._athena_groups.keys())[8:17]

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names, title = "Cianobacteria Standards")
plt.show()

# https://vimeo.com/340216087 15:05 compare readings to standars 
group_keys = ['ozun', 'lshy','hqlr','tscd']

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Compare Readings to Standards")
plt.show()

# https://vimeo.com/340216087 15:05 compare readings to standars 
# Earliest measure compared to Au3Cl
group_keys = ['ozun', 'tscd']

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Compare AuCl to 0.12")
plt.show()

# https://vimeo.com/340216087 15:05 compare readings to standars 
# final measure compared to Au foil
group_keys = ['lshy','hqlr']

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Compare Au Foil to 720")
plt.show()
