# import Document library from Chem Data Extractor
from chemdataextractor import Document
# import library for managing files
from pathlib import Path
import sys

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
def get_common(current_pattern, file_name):
    return ""

def get_max(uniques):
    max_val = 0
    max_lbl = ""
    for chement in uniques:
        if uniques[chement] > max_val:
            max_val = uniques[chement]
            max_lbl = chement.replace('\n',' ')
    return max_lbl, max_val

def cde_read_pdfs(argv, pdf_path = "./pdfs"):
    try:
        pdf_path = argv[0]
        name_pattern = argv[1]
    except:
        print("missing arguments"+
              "\n -string pdf files path")
        return
    pdf_dir= Path(pdf_path)
    files_list = get_files_list(pdf_dir, name_pattern)
    for a_file in files_list:
        file_name = a_file.name
        pattern_current = ""
        if pattern_current == "":
            pattern_current = file_name
        else:
            pattern_current = common(pattern_current, file_name)
        print(file_name)
        
if __name__ == "__main__":
   cde_read_pdfs(sys.argv[1:])

