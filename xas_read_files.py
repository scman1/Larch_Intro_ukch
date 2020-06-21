# script for bulk processing
# read xas text files from a directory
# create basic plots of energy and fitting
# merge groups and save as athena

# import library for managing files
from pathlib import Path
import sys
from difflib import SequenceMatcher

# import library for managing csv files
import csv

# get the data from the csv_file, assuming first column is integer id
def get_csv_data(input_file, id_field):
    csv_data = {}
    fieldnames=[]
    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if fieldnames==[]:
                fieldnames=list(row.keys())
            csv_data[int(row[id_field])]=row
    return csv_data, fieldnames

# writes data to the given file name
def write_csv_data(values, filename):
    fieldnames = []
    for item in values.keys():
        for key in values[item].keys():
            if not key in fieldnames:
                fieldnames.append(key)
    #write back to a new csv file
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key in values.keys():
            writer.writerow(values[key])

# add logging
# save processing steps in log file
import logging

file_dir = Path("log_dir")
if not file_dir.exists():
    file_dir.mkdir(parents=True)

logging.basicConfig(format='[%(asctime)s] %(message)s', filename='log_dir/processing.log', filemode='w', level=logging.INFO)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

# add a configuration file to read parameters
import configparser

# files for processing xas data
import numpy as np
import pylab
# ploting library
import matplotlib.pyplot as plt
import larch

# larch-xas processing functions
from larch_plugins.xafs import autobk, xftf

# import the larch.io function for merging groups interpolating if necessary
from larch.io import merge_groups

# import the larch.io libraries for managing athena files
from larch.io import create_athena, read_athena, extract_athenagroup

# create a larch interpreter, passed to most functions
my_larch = larch.Interpreter()

# get_file_groups:
# groups files in a directory according to common text patterns assuming that 
# patterns correspond to files which must be processed together. 
# input:
#  - the directory path where the files to be processed are placed
#  - the string which is used to filter the files to process (wildcad caracheters)
# output:
#  - an indexed list of files, which uses common patterns found as keys

def get_file_groups(source_dir, filename_pattern, group_files):
    #i_counter = 0
    files_list = []
    for filepath in sorted(source_dir.glob(filename_pattern)):
        #i_counter += 1
        files_list.append(filepath)
    log_message = "Found " + str(len(files_list)) + \
        " files to process with pattern: " + \
        filename_pattern + \
        "in dir: " + str(source_dir) + " Group: " +str(group_files)
    logging.info(log_message)
    pattern_current = ""
    file_groups = {}

    if group_files:
        for index, a_file in enumerate(files_list):
            file_name = a_file.name
            if index < len(files_list)-1:
                another_file = files_list[index+1].name
                common_pattern = get_common(file_name, another_file)
                if pattern_current == "":
                    pattern_current = common_pattern
                    if pattern_current in file_groups:
                        file_groups[pattern_current] += [file_name, another_file]
                    else:
                        file_groups[pattern_current] = [file_name, another_file]
                elif pattern_current == common_pattern:
                    file_groups[pattern_current].append(another_file)
                else:
                    pattern_current = ""
            else:
                another_file = ""

        patterns_found = file_groups.keys()
        patterns_remove = []

        # mark redundant patterns to remove
        for pattern_check in file_groups:
            for pattern_other in patterns_found:
                if pattern_other != pattern_check:
                    # use sets to ignore duplicates
                    set_check = set(file_groups[pattern_check])
                    set_other = set(file_groups[pattern_other])
                    if pattern_other in pattern_check:
                        #merge check into other and remove check
                        file_groups[pattern_other] = list(set_other.union(set_check))
                        patterns_remove.append(pattern_check)
                    elif pattern_other in pattern_check:
                        #merge other into check and remove other
                        patterns_remove.append(pattern_other)
                        file_groups[pattern_check] = list(set_other.union(set_check))
    
        # remove redundant patterns
        if len(patterns_remove) > 0:
            for pattern_rem in patterns_remove:
                file_groups.pop(pattern_rem)
    else:
        file_groups['unique'] = []
        for file in files_list:
            file_groups['unique'].append(file.name)
        
    for pattern in file_groups: 
        log_message = str(len(file_groups[pattern])) + " files to process with pattern " + pattern + "\nFiles: " + str(file_groups[pattern])
        logging.info(log_message)
    return file_groups


# get the common pattern between the two strings
def get_common(file_1, file_2):
    common_pattern = ""
    seqMatch = SequenceMatcher(None, file_1, file_2)
    match = seqMatch.find_longest_match(0, len(file_1), 0, len(file_2))
    if (match.size!=0): 
        common_pattern = file_1[match.a: match.a + match.size]
    else: 
         log_message = 'No longest common sub-string found'
         logging.info(log_message)
    return common_pattern

# basic plot of a group
# input:
#   - a larch xas group
#   - the detination dir (where to save)
def basic_plot(xas_group, dest_dir, show_plot = False):
    fig=plt.figure(figsize=(10,8))
    plt.tick_params(axis='both', labelsize=6)
    
    # plot grid of results:
    # mu + bkg
    plt.subplot(2, 2, 1)
    plt.title('$\mu$ and background')
    plt.plot(xas_group.energy, xas_group.bkg, 'r--', label = 'background')
    plt.plot(xas_group.energy, xas_group.mu, label = "$\mu$")
    plt.xlabel('Energy (eV)')
    plt.grid(linestyle=':', linewidth=1)
    plt.legend()

    # normalized XANES
    # find array bounds for normalized mu(E) for [e0 - 25: e0 + 75]
    j0 = np.abs(xas_group.energy-(xas_group.e0 - 25.0)).argmin()
    j1 = np.abs(xas_group.energy-(xas_group.e0 + 75.0)).argmin()

    
    plt.subplot(2, 2, 2)
    plt.title('normalized $\mu$')
    plt.plot(xas_group.energy[j0:j1], xas_group.norm[j0:j1], label="$\mu$ Normalised")
    plt.xlabel('Energy (eV)')
    plt.grid(linestyle=':', linewidth=1)
    plt.legend()
    
    # chi(k)
    plt.subplot(2, 2, 3)
    plt.title(r"$\chi(k)$")
    plt.plot(xas_group.k, xas_group.chi*xas_group.k**2, label= r'$ \chi(k^2)$')
    plt.plot(xas_group.k, xas_group.kwin, 'r--', label= r'$k$ window')
    plt.xlabel(r'$ k (\AA^{-1}) $', fontsize='small')
    plt.ylabel(r'$ k^2 \chi(\AA^{-2}) $', fontsize='small')
    plt.grid(linestyle=':', linewidth=1)
    plt.legend()

    # chi(R)
    plt.subplot(2, 2, 4)
    plt.title(r"$\chi(R)$")
    plt.plot(xas_group.r, xas_group.chir_mag, label = r"$\chi(R)$ magnitude")
    plt.plot(xas_group.r, xas_group.chir_re, 'r--', label = r"$\chi(R)$ re")
    plt.xlabel(r'$ R (\AA) $',fontsize='small')
    plt.ylabel(r'$ \chi(R) (\AA^{-3}) $', fontsize='small')
    plt.grid(linestyle=':', linewidth=1)
    plt.legend()

    
    save_as = dest_dir / (xas_group.label + "_01.jpg")

    if not save_as.parent.exists():
        save_as.parent.mkdir(parents=True)
        
    fig.tight_layout(pad=3.0)
    fig.suptitle(xas_group.label)    
    plt.savefig(str(save_as))
    if show_plot:
        plt.show()
    plt.clf()

# save energy and normalised mu
def save_e_nmu(xafsgroup, save_dir):
    export = {}
    for n_index, value in enumerate(xafsgroup.energy):
        export[n_index] = {'energy':value, 'norm':xafsgroup.norm[n_index]}
    write_csv_data(export, save_dir/(xafsgroup.label+"_EvNm.csv"))

# xas_read_files
# groups files in a directory according to common text patterns
# process groups of files using larch with defaults:
#   get mu (calculate if needed)
#   get normal, pre-edge and post-edge E0
#   plot groups
#   merge groups
#   plot merge
#   save diagrams
#   save all as athena project

def xas_read_files(argv):
    try:
        # required
        file_path = argv[0]
        name_pattern = argv[1]
        if len(argv) == 3:
            if argv[2] == 'T':
                group_files = True
        else:
            group_files = False
            
    except:
        print("missing arguments"+
              "\n -string files path (eg: ../documents/ascii_path)"+
              "\n -string file pattern (eg: *experiment_FeO2_sample*)"+
              "\n -character T or F to indicate if grouping")
        return
    
    file_dir= Path(file_path)
    # initialisation file
    ini_file = file_dir / "xas_processing.ini"
    data_columns = {}
    if ini_file.exists():
        log_message = 'reading ini file: '+ str(ini_file)
        logging.info(log_message) 
        xas_config = configparser.ConfigParser()
        xas_config.read(ini_file)
        # especify data columns
        if 'col_energy' in xas_config['Data_Columns']:
            data_columns['energy'] = int(xas_config['Data_Columns']["col_energy"])
        if 'col_time' in xas_config['Data_Columns']:
            data_columns['time'] = int(xas_config['Data_Columns']["col_time"])
        if 'col_i0' in xas_config['Data_Columns']:
            data_columns['i0'] = int(xas_config['Data_Columns']["col_i0"])
        if 'col_it' in xas_config['Data_Columns']:
            data_columns['it'] = int(xas_config['Data_Columns']["col_it"])
        if 'col_ir' in xas_config['Data_Columns']:        
            data_columns['ir'] = int(xas_config['Data_Columns']["col_ir"])
        if 'col_mu' in xas_config['Data_Columns']:
            data_columns['mu'] = int(xas_config['Data_Columns']["col_mu"])
        if 'col_mur' in xas_config['Data_Columns']:
            data_columns['mur'] = int(xas_config['Data_Columns']["col_mur"])
    else:
        log_message = 'missing ini file: '+ str(ini_file)
        logging.info(log_message)
        log_message = 'using defaults'
        logging.info(log_message)
        # especify data columns
        # energy 0
        # time   1
        # i0     2
        # it     3
        # ir     4
        # mu     5
        # mur    6
        data_columns = {'energy':0, 'time':1,'i0':2,'it':3,'ir':4, 'mu':5,
                        'mur': 6}
    log_message = "data_columns:" + str(data_columns)
    logging.info(log_message)
      
    file_groups = get_file_groups(file_dir, name_pattern, group_files)
    
    # process file groups
    for pattern in file_groups:
        groups = []
        for file in file_groups[pattern]:
            save_dir = file_dir / 'result' / pattern[1:][:-1]
            file_path = file_dir / file
            xafsdat = larch.io.read_ascii(file_path)
            # get data columns specified in ini_file
            if 'energy' in data_columns:
                xafsdat.energy = xafsdat.data[data_columns['energy']]
            if 'time' in data_columns:
                xafsdat.time = xafsdat.data[data_columns['time']]
            if 'i0' in data_columns:    
                xafsdat.i0 = xafsdat.data[data_columns['i0']]
            if 'it' in data_columns:
                xafsdat.it = xafsdat.data[data_columns['it']]
            if 'ir' in data_columns:
                xafsdat.ir = xafsdat.data[data_columns['ir']]
            if 'mu' in data_columns:
                xafsdat.mu = xafsdat.data[data_columns['mu']]
            if 'mur' in data_columns:
                xafsdat.mue = xafsdat.data[data_columns['mur']]
            # run autobk on the xafsdat Group, including a larch Interpreter....
            # note that this expects 'energy' and 'mu' to be in xafsdat, and will
            # write data for 'k', 'chi', 'kwin', 'e0', ... into xafsdat
            autobk(xafsdat, rbkg=1.0, kweight=2, _larch=my_larch)
            
            # Fourier transform to R space, again passing in a Group (here,
            # 'k' and 'chi' are expected, and writitng out 'r', 'chir_mag',
            # and so on
            xftf(xafsdat, kmin=2, kmax=15, dk=3, kweight=2, _larch=my_larch)

            xafsdat.label = xafsdat.filename[:-4]
            
            # add group to list
            groups.append(xafsdat)

            # plot and save each file in group 
            basic_plot(xafsdat, save_dir)

            # save energy v normalised mu
            save_e_nmu(xafsdat, save_dir)
        # merge groups
        merged_group = merge_groups(groups)
        merged_group.label = pattern[1:][:-1] + "_merge"
        autobk(merged_group, rbkg=1.0, kweight=2, _larch=my_larch)
        xftf(merged_group, kmin=2, kmax=15, dk=3, kweight=2, _larch=my_larch)
        # plot and save for merge
        basic_plot(merged_group, save_dir)
        # save energy v normalised mu for merge
        save_e_nmu(merged_group, save_dir)
        
        groups.append(merged_group)

        log_message = "Processed groups (including merge): " + str(len(groups)) + " for pattern " + (pattern[1:][:-1])
        logging.info(log_message)
        
        # save as an athena project

        project_name = save_dir / (pattern[1:][:-1] + '.prj')
        athena_project = create_athena(project_name)
        for a_group in groups:
            athena_project.add_group(a_group)
        athena_project.save()
        log_message = "Saved athena project " + str(project_name)
    

if __name__ == "__main__":
   xas_read_files(sys.argv[1:])

