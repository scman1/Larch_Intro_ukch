# import Document library from Chem Data Extractor
from chemdataextractor import Document
# import library for managing files
from pathlib import Path
import sys
from difflib import SequenceMatcher 


def get_files_list(source_dir, filename_pattern):
    i_counter = 0
    files_list = []
    for filepath in sorted(source_dir.glob(filename_pattern)):
        i_counter += 1
        files_list.append(filepath)
    return files_list

def get_uniques(cde_doc):
    uniques={}
    for chement in cde_doc.cems:
        if not chement.text in uniques:
            uniques[chement.text] = 1
        else:
            uniques[chement.text] += 1
    return uniques

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

def xas_read_files(argv, pdf_path = "ascii"):
    try:
        pdf_path = argv[0]
        name_pattern = argv[1]
    except:
        print("missing arguments"+
              "\n -string files path")
        return
    pdf_dir= Path(pdf_path)
    files_list = get_files_list(pdf_dir, name_pattern)
    pattern_current = ""
    common_files = {}
    for index, a_file in enumerate(files_list):
        file_name = a_file.name
        if index < len(files_list)-1:
            another_file = files_list[index+1].name
            common_pattern = get_common(file_name, another_file)
            #common_pattern = file_name[match.a: match.a + match.size]
            if pattern_current == "":
                pattern_current = common_pattern
                common_files[pattern_current] = [file_name, another_file]
            elif pattern_current == common_pattern:
                common_files[pattern_current].append(another_file)
            else:
                pattern_current = ""
        else:
            another_file = ""
            #pattern_current = get_common(pattern_current, file_name)
        print(file_name, pattern_current)
    for pattern in common_files:
        print(len(common_files[pattern]), "files with pattern" , pattern, "\nFiles: ", common_files[pattern])
    
if __name__ == "__main__":
   xas_read_files(sys.argv[1:])

