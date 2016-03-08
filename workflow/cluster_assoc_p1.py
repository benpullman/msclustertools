#!/usr/bin/python

import sys
#import cluster_base
import glob
import os

if len(sys.argv) == 3:
    [program_name,path_to_data,datalist_output] = sys.argv
else:
    sys.exit('Error!  Wrong number of inputs specified.')


def generate_datalist(path_to_data):
    path_to_datalist = ""
    datalist = [
        filename
        for filename in glob.glob(path_to_data + "/*")
    ]
    with open(datalist_output, 'w') as datalist_file:
        datalist_file.write('\n'.join(datalist))
        path_to_datalist = datalist_file.name

generate_datalist(path_to_data)
