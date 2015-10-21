#!/usr/bin/python

import sys
import cluster_base
import glob
import os

if len(sys.argv) == 3:
    [program_name,phase1_result_mgf,phase1_result_mgf_altered] = sys.argv
else:
    sys.exit('Error!  Wrong number of inputs specified.')

def process_phase1_result_mgf(phase1_result_mgf,phase1_result_mgf_altered):
    for filename in glob.glob(phase1_result_mgf + "/*.mgf"):
        with open(filename, 'r') as readfile:
            with open(phase1_result_mgf_altered + "/" + os.path.split(filename)[1], 'w') as writefile:
                for line in readfile:
                    if "CLUSTER_SIZE=" in line:
                        pass
                    else:
                        writefile.write(line.replace("TITLE=out_0_0.","SCANS="))

process_phase1_result_mgf(phase1_result_mgf,phase1_result_mgf_altered)
