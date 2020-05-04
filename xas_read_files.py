# read xas text files from a directory

# import library for managing files
from pathlib import Path
import sys
from difflib import SequenceMatcher 

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
        pdf_path = argv[0]
        name_pattern = argv[1]
    except:
        print("missing arguments"+
              "\n -string files path (eg: ../documents/ascii_path)"+
              "\n -string file pattern (eg: *experiment_FeO2_sample*)")
        return
    pdf_dir= Path(pdf_path)
    files_list = get_file_groups(pdf_dir, name_pattern)


if __name__ == "__main__":
   xas_read_files(sys.argv[1:])

