import os
import glob

def find_ext_files(directory, ext_file):
    csv_files = []
    for file in glob.glob(os.path.join(directory, '*.' + ext_file)):
        csv_files.append(os.path.abspath(file))
        
    return csv_files