#!/usr/bin/python

import sys
import glob
import os
import bisect

import cluster_base

if len(sys.argv) == 7:
    [program_name,clust_in,tsv_clust_in,tsv_spec_in,params,clustOrSpec,tsv_out_file] = sys.argv
else:
    sys.exit('Error!  Wrong number of inputs specified.')

file_id_col = 2

print(str(sys.argv))

isClust = clustOrSpec=="--clust"
parameters = cluster_base.parse_xml_file(params)
# class MiniSpectrum(object):
#     def __init__(self,file_name,scan,cluster,sqs,charge,precursor_mz):
#         self.file_name = file_name
#         self.scan = scan
#         self.cluster = cluster
#         self.sqs = sqs
#         self.charge = charge
#         self.precursor_mz = precursor_mz
h_ion = 1.007825035
mass_mh = 18.010564686 + h_ion
aa_weights = {
    'A':71.037113787,
    'R':156.101111026,
    'D':115.026943031,
    'N':114.042927446,
    'C':160.030648200,
    'E':129.042593095,
    'Q':128.058577510,
    'G':57.021463723,
    'H':137.058911861,
    'I':113.084063979,
    'L':113.084063979,
    'K':128.094963016,
    'M':131.040484605,
    'F':147.068413915,
    'P':97.052763851,
    'S':87.032028409,
    'T':101.047678473,
    'W':186.079312952,
    'Y':163.063328537,
    'V':99.068413915
}

def extract_middle(peptide_string):
    first_letter = peptide_string.split(".",1)[0]
    last_letter = peptide_string.rsplit(".",1)[1]
    return ''.join([peptide_string.rsplit(".",1)[0].split(".",1)[1]])

def calculate_mass(peptide_string, charge):
    stripped = extract_middle(peptide_string)
    weights = [aa_weights[aa] for aa in stripped if aa.isalpha()]
    mods = [float(mod) for mod in ''.join([aa for aa in stripped if not aa.isalpha()]).split('+') if mod != '']
    return (sum(weights) + sum(mods) + (charge-1)*h_ion + mass_mh)/charge

def isotopes(peptide_string,charge,n=4):
    peptide_mass = calculate_mass(peptide_string, charge)
    return [
        peptide_mass + i*(1/charge)
        for i in range(0,n+1)
    ]

def calculate_precursor_mix(cluster,spectra,tolerance):
    tolerance_percent = tolerance/1000000
    pep_rep = cluster.purity.representative_spectrum
    mix = 0
    if pep_rep != 'PEPTIDE':
        rep_isotopes = isotopes(pep_rep,int(cluster.charge))
        similar = 0
        total = len(spectra)
        for spectrum in spectra:
            for isotope in rep_isotopes:
                if (float(spectrum.precursor_mz) > (isotope - tolerance_percent*isotope) and float(spectrum.precursor_mz) < (isotope + tolerance_percent*isotope)):
                    similar += 1
        mix = similar/total
    return mix

def parseSpectra(cluster, spectra_string):
    spectra_list = []
    spectra = spectra_string.split(",")
    for spectrum in spectra:
        parts = spectrum.split(":")
        file_name = parts[0]
        scan_number = parts[1]
        sqs = parts[2]
        charge = parts[3]
        precursor_mz = parts[4]
        spectra_list.append(cluster_base.Spectrum(scan_number,file_name,precursor_mz,charge,None,None,sqs,None))
    return spectra_list

def read_tsv_result_file(tsv_file_name, tsv_file_id_col, is_cluster):
    header = ""
    tsv_lines = {}
    peptides = {}
    with open(tsv_file_name, 'r') as tsv_file:
        header = tsv_file.readline()
        for line in tsv_file:
            split_line = line.replace('\n','').split('\t')
            scan_number = int(split_line[tsv_file_id_col])
            peptide = split_line[7]
            file_name = split_line[0]
            if is_cluster:
                id_for_dict = file_name + ":" + str(scan_number)
            else:
                id_for_dict = file_name + ":" + str(scan_number)
                peptides[id_for_dict] = peptide
                # print(id_for_dict)
            tsv_lines[id_for_dict] = split_line
    return tsv_lines, header, peptides

clustered_lines, clustered_header, clustered_peptides = read_tsv_result_file(tsv_clust_in, 1, True)
unclustered_lines, unclustered_header, unclustered_peptides = read_tsv_result_file(tsv_spec_in, 2, False)

# print(unclustered_lines)

if isClust:
    modified_header = clustered_header.replace("\n","") + "\tCluster Size\tCluster Purity\tRepresentative Spectrum\tCluster SQS\tAverage Spectra SQS\tMix Score\tSpectra\tPrecursor Mix"
else:
    modified_header = unclustered_header.replace("\n","") + "\tCluster\tSpectra SQS"

class parseMGF(object):

    def __init__(self, file_name):
        self.file_name = file_name
        self.cluster = None
        self.cluster_sqs = None
        self.avg_spectra_sqs = None
        self.spectra = None
        self.pepmass = None
        self.charge = None
        self.spec_index = None
        self.cluster_size = None
        self.purity = None
        self.mix_score = None
        self.precursor_mix = None
        #self.assigned_charge = None
    def as_array(self):
        return [str(self.cluster_size), str(self.purity.purity_no_undentified), str(self.purity.representative_spectrum), str(self.cluster_sqs),str(self.avg_spectra_sqs),str(self.mix_score),str(self.spectra),str(self.precursor_mix)]
    def as_no_id_tsv(self):
        return [
            self.file_name, # #SpecFile
            str(self.spec_index), # SpecIndex
            str(self.cluster), # Scan#
            parameters['fragmentation.fragmentation'][0], # FragMethod
            str(self.pepmass), # Precursor
            "0", # PMError(ppm)
            str(self.charge), # Charge
            "", # Peptide
            "",# Protein
            "0", # DeNovoScore
            "0", # MSGFScore
            "1", # SpecProb
            "1", # P-value
            "1", # FDR
            "1"  # PepFDR
        ]
    def as_no_id_tsv_unclustered(self, scan, file_name, precursor_mz, charge):
        return [
            file_name, # #SpecFile
            str(scan), # SpecIndex
            str(scan), # Scan#
            parameters['fragmentation.fragmentation'][0], # FragMethod
            str(precursor_mz), # Precursor
            "0", # PMError(ppm)
            str(charge), # Charge
            "", # Peptide
            "",# Protein
            "0", # DeNovoScore
            "0", # MSGFScore
            "1", # SpecProb
            "1", # P-value
            "1", # FDR
            "1"  # PepFDR
        ]


mgf = []
mgf_files = glob.glob(clust_in + "/*.mgf")
server_file_name = cluster_base.upload_file_mapping(params = parameters)
print("server_file_name: " + str(server_file_name))
server_file_name_swap = cluster_base.swap(server_file_name)
with open(tsv_out_file, 'w') as tsv_out:
    tsv_out.write(modified_header + "\n")
    for mgf_file_name in mgf_files:
        index = 1
        with open(mgf_file_name, 'r') as mgf_file:
            for line in mgf_file:
                if "BEGIN" in line:
                    current_mgf = parseMGF(file_name = os.path.split(mgf_file_name)[1])
                elif "SCANS" in line:
                    current_mgf.cluster = int(line.replace("SCANS=","").replace("\n",""))
                elif "CLUSTER_SQS" in line:
                    current_mgf.cluster_sqs = line.replace("CLUSTER_SQS=","").replace("\n","")
                elif "AVG_SPECTRUM_SQS" in line:
                    current_mgf.avg_spectra_sqs = line.replace("AVG_SPECTRUM_SQS=","").replace("\n","")
                elif "PEPMASS" in line:
                    current_mgf.pepmass = line.replace("PEPMASS=","").replace("\n","")
#                elif "ASSIGNED_CHARGE" in line:
#                    current_mgf.assigned_charge = line.replace("ASSIGNED_CHARGE=","").replace("\n","")
                elif "CHARGE" in line:
                    current_mgf.charge = line.replace("CHARGE=","").replace("\n","")
                elif "SPECTRA" in line:
                    current_mgf.spectra = line.replace("SPECTRA=","").replace("\n","")
                elif "END" in line:
                    current_mgf.spec_index = index

                    ## add info to
                    spectra = parseSpectra(current_mgf.cluster, current_mgf.spectra)
                    current_mgf.cluster_size = len(spectra)

                    peptides = []


                    for spectrum in parseSpectra(current_mgf.cluster, current_mgf.spectra):
                        spec_id = spectrum.file_id + ":" + spectrum.scan_number
                        # try:
                        #     spec_id = server_file_name[spectrum.file_id] + ":" + spectrum.scan_number
                        # except:
                        #     pass
                        spectrum.cluster_number = current_mgf.cluster

                        try:
                            peptides.append(unclustered_peptides[spec_id])
                        except:
                            peptides.append(None)
                        if not isClust:
                            tsv_start_unclustered = unclustered_lines.get(spec_id, current_mgf.as_no_id_tsv_unclustered(spectrum.scan_number,spectrum.file_id, spectrum.precursor_mz, spectrum.charge))
                            output_unclustered = tsv_start_unclustered + [str(spectrum.cluster_number), str(spectrum.sqs)]
                            tsv_out.write('\t'.join(output_unclustered) + '\n')
                    if isClust:
                        current_mgf.purity = cluster_base.purity_from_peptide_list(peptides)
                        current_mgf.mix_score = cluster_base.mix_score(spectra)
                        alt_mixture = calculate_precursor_mix(current_mgf, spectra, 30)
                        current_mgf.precursor_mix = alt_mixture
                        tsv_start_clustered = clustered_lines.get(os.path.split(mgf_file_name)[1] + ":" + str(index),current_mgf.as_no_id_tsv())
                        output_clustered = tsv_start_clustered + current_mgf.as_array()
                        tsv_out.write('\t'.join(output_clustered) + '\n')
                    index += 1
