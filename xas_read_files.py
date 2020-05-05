# read xas text files from a directory

# import library for managing files
from pathlib import Path
import sys
from difflib import SequenceMatcher

# files for processing xas data
import numpy as np
import pylab
# ploting library
import matplotlib.pyplot as plt
import larch

# now import larch-specific Python code
from larch_plugins.xafs import autobk, xftf

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

def get_file_groups(source_dir, filename_pattern):
    #i_counter = 0
    files_list = []
    
    
    for filepath in sorted(source_dir.glob(filename_pattern)):
        #i_counter += 1
        files_list.append(filepath)
    
    print("Found ", len(files_list), "to process with pattern:", filename_pattern)
    pattern_current = ""
    file_groups = {}
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
    
    for pattern in file_groups:
        print(len(file_groups[pattern]), "files to process with pattern" , pattern, "\nFiles: ", file_groups[pattern])
        
    return file_groups


# get the common pattern between the two strings
def get_common(file_1, file_2):
    common_pattern = ""
    seqMatch = SequenceMatcher(None, file_1, file_2)
    match = seqMatch.find_longest_match(0, len(file_1), 0, len(file_2))
    if (match.size!=0): 
        #print (match)
        common_pattern = file_1[match.a: match.a + match.size]
    else: 
         print ('No longest common sub-string found') 
    return common_pattern

# xas_read_files
# groups files in a directory according to common text patterns
# process groups of files using larch with defaults:
#   get mu (calculate if needed)
#   get normal, pre-edge and post-edge E0
#   merge groups
#   save diagrams

def xas_read_files(argv):
    try:
        # required
        file_path = argv[0]
        name_pattern = argv[1]
    except:
        print("missing arguments"+
              "\n -string files path (eg: ../documents/ascii_path)"+
              "\n -string file pattern (eg: *experiment_FeO2_sample*)")
        return
    file_dir= Path(file_path)
    file_groups = get_file_groups(file_dir, name_pattern)
    print(file_groups)
    groups = []
    for pattern in file_groups:
        for file in file_groups[pattern]:
            file_path = file_dir / file
            xafsdat = larch.io.read_ascii(file_path)
            # especify data columns
            # energy 0
            # time   1
            # i0     2
            # it     3
            # ir     4
            # mu     5
            # mur    6
            xafsdat.energy = xafsdat.data[0]
            xafsdat.time =   xafsdat.data[1]
            xafsdat.i0 =     xafsdat.data[2]
            xafsdat.it =     xafsdat.data[3]
            xafsdat.ir =     xafsdat.data[4]
            xafsdat.mu =     xafsdat.data[5]
            xafsdat.mue =    xafsdat.data[6]
            # run autobk on the xafsdat Group, including a larch Interpreter....
            # note that this expects 'energy' and 'mu' to be in xafsdat, and will
            # write data for 'k', 'chi', 'kwin', 'e0', ... into xafsdat
            autobk(xafsdat, rbkg=1.0, kweight=2, _larch=my_larch)
            
            # Fourier transform to R space, again passing in a Group (here,
            # 'k' and 'chi' are expected, and writitng out 'r', 'chir_mag',
            # and so on
            xftf(xafsdat, kmin=2, kmax=15, dk=3, kweight=2, _larch=my_larch)

            # add group to list
            groups.append(xafsdat)
            
            #
            # plot grid of results:
            # mu + bkg
            plt.subplot(2, 2, 1)
            plt.plot(xafsdat.energy, xafsdat.bkg, 'r--')
            plt.plot(xafsdat.energy, xafsdat.mu)
            plt.xlabel('Energy (eV)')

            # normalized XANES
            # find array bounds for normalized mu(E) for [e0 - 25: e0 + 75]
            j0 = np.abs(xafsdat.energy-(xafsdat.e0 - 25.0)).argmin()
            j1 = np.abs(xafsdat.energy-(xafsdat.e0 + 75.0)).argmin()

            plt.subplot(2, 2, 2)
            plt.plot(xafsdat.energy[j0:j1], xafsdat.norm[j0:j1])
            plt.xlabel('Energy (eV)')

            # chi(k)
            plt.subplot(2, 2, 3)
            plt.plot(xafsdat.k, xafsdat.chi*xafsdat.k**2)
            plt.plot(xafsdat.k, xafsdat.kwin, 'r--')
            plt.xlabel(r'$ k (\AA^{-1}) $')
            plt.ylabel(r'$ k^2 \chi(\AA^{-2}) $')

            # chi(R)
            plt.subplot(2, 2, 4)
            plt.plot(xafsdat.r, xafsdat.chir_mag)
            plt.plot(xafsdat.r, xafsdat.chir_re, 'r--')
            plt.xlabel(r'$ R (\AA) $')
            plt.ylabel(r'$ \chi(R) (\AA^{-3}) $')
 
            save_as = file_dir / 'processed' / pattern[1:][:-1] / (file[:-4] + ".jpg")
            if not save_as.parent.exists():
                save_as.parent.mkdir()
                
            plt.savefig(str(save_as))
            #pylab.show()

    print("Processed groups for pattern:", len(groups), groups)
    

if __name__ == "__main__":
   xas_read_files(sys.argv[1:])

