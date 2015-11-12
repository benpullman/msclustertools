#!/usr/bin/python

import sys
import cluster_base
import glob
import os


if len(sys.argv) == 6:
    [program_name,phase1_result_clust,phase2_result_clust,phase1_result_mgf,datalist,output_mgf] = sys.argv
else:
    print(sys.argv)
    sys.exit('Error!  Wrong number of inputs specified.')

data_map = cluster_base.datalist_to_map(datalist)
clusters = cluster_base.initialize_unidentified_clusters(phase1_result_clust,phase2_result_clust,"out_0_0.",data_map)

print(data_map)

def process_final_result_mgf(output_mgf, phase1_result_mgf, clusters):
    # print(data_map)
    # for c in clusters:
    #     print(c)
    for filename in glob.glob(phase1_result_mgf + "/*.mgf"):
        with open(filename, 'r') as readfile:
            is_ions = False
            current_cluster = None
            with open(output_mgf + "/" + os.path.split(filename)[1], 'w') as writefile:
                for line in readfile:
                    if "BEGIN" in line:
                        is_ions = False
                    if "TITLE" in line:
                        current_cluster_id = int(line.replace("TITLE=out_0_0.",""))
                        current_cluster = clusters[str(current_cluster_id)]
                    elif line[0].isdigit() and not is_ions:
                        is_ions = True
                        writefile.write('\n'.join([
                            'SCANS=' + str(current_cluster.scan_number),
                            #'ASSIGNED_CHARGE=' + str(current_cluster.charge),
                            'CLUSTER_SQS=' + str(current_cluster.sqs),#cluster sqs
                            'AVG_SPECTRUM_SQS=' + str(current_cluster.avg_spectra_sqs()),#average spectrum sqs
                            'SPECTRA=' + str(','.join([
                                spectrum.output_pair()
                                for spectrum in current_cluster.spectra
                            ])),#spectrum list (as string)
                            line
                        ]))
                    else:
                        writefile.write(line)

process_final_result_mgf(output_mgf, phase1_result_mgf, clusters)
