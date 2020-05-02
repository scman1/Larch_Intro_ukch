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
    print(files_list)

if __name__ == "__main__":
   cde_read_pdfs(sys.argv[1:])

