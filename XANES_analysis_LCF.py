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
# |         Athena recalculates everything so we      | #
# |      need to create a function that calculates    | #
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

 #######################################################
# |       The code for plotting Nmu vs E repeats      | #
# |   so it is useful to have a plotting function     | #
# V            to reduce duplicated code              V #
 #######################################################
# plot mu vs flat normalised mu for selected groups
def plot_NxmuE_E_athena_prj(athena_project, group_keys, group_names,
                            title = "Normalised Mu vs E", xlimits = None,
                            ylimits = None):    
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
    if xlimits != None:
        plt.xlim(xlimits[0],xlimits[1])
    if ylimits != None:
        plt.ylim(ylimits[0],ylimits[1])
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
plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Cianobacteria Readings",
                              xlimits = [11860,12000])
plt.show()

# https://vimeo.com/340216087 12:10 plot standards 
# get the group keys for first last 9 groups
group_keys = list(cianobacteria_project._athena_groups.keys())[8:17]

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Cianobacteria Standards",
                              xlimits = [11860,12000])
plt.show()

# https://vimeo.com/340216087 15:05 compare readings to standars 
group_keys = ['ozun', 'lshy','hqlr','tscd']

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Compare Readings to Standards",
                              xlimits = [11860,12000])
plt.show()

# https://vimeo.com/340216087 15:05 compare readings to standars 
# Earliest measure compared to Au3Cl
group_keys = ['ozun', 'tscd']

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Compare AuCl to 0.12",
                              xlimits = [11860,12000])
plt.show()

# https://vimeo.com/340216087 15:05 compare readings to standars 
# final measure compared to Au foil
group_keys = ['lshy','hqlr']

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, group_keys, group_names,
                              title = "Compare Au Foil to 720",
                              xlimits = [11860,12000])
plt.show()

#Linear Combination Fitting

# https://vimeo.com/340216087 19:58 linear combination for inner state 
standard_keys = ['hqlr','tscd'] #Au Foil and Au3Cl
intermidate_state_key = 'd_7_03'

# get the intermediate group
mid_group = gr_0 = extract_athenagroup(cianobacteria_project._athena_groups[intermidate_state_key])
# recalculate normalisation
calc_with_defaults(mid_group)

# get the list of standard groups
components = {}
for group_key in standard_keys:
    components[group_key] = extract_athenagroup(cianobacteria_project._athena_groups[group_key])
    # recalculate normalisation
    calc_with_defaults(components[group_key])
    
# perform linear combination fitting
comb = lincombo_fit(mid_group,list(components.values()),[0.5,0.5])

plt = plot_NxmuE_E_athena_prj(cianobacteria_project, standard_keys, group_names,
                              title = "Compare to LCF 2 standards",
                              xlimits = [11900,12000],
                              ylimits = [0.0,1.2])
#add mid_group and LCF result to plot
plt.plot(mid_group.energy, mid_group.flat, label="7.03",color="blue")
plt.plot(comb.xdata, comb.ydata, label="LCF result",color="r")
plt.legend() # needed for showing legends of last two lines
plt.show()

# https://vimeo.com/340216087 23:50 add another group 
standard_keys = ['hqlr','tscd', 'qhxp'] #Au Foil, Au3Cl, and Au sulphide
components['qhxp'] = extract_athenagroup(cianobacteria_project._athena_groups['qhxp'])
calc_with_defaults(components['qhxp'])
   
# perform linear combination fitting
comb = lincombo_fit(mid_group,list(components.values()),[0.333,0.333,0.333])


plt = plot_NxmuE_E_athena_prj(cianobacteria_project, standard_keys, group_names,
                              title = "Compare to LCF 3 standards",
                              xlimits = [11900,12000],
                              ylimits = [0.0,1.2])
#add mid_group and LCF result to plot
plt.plot(mid_group.energy, mid_group.flat, label="7.03",color="blue")
plt.plot(comb.xdata, comb.ydata, label="LCF result",color="r")
plt.legend() # needed for showing legends of last two lines
plt.show()
