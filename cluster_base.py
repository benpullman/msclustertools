import sys
import glob
import os
import bisect
from itertools import groupby

class SpectrumBase(object):

    """Base class for all spectrum-like objects"""

    def __init__(self, scan_number, file_id, precursor_mz, charge):
        self.scan_number = scan_number
        self.file_id = file_id
        self.precursor_mz = precursor_mz
        self.charge = charge

    def __eq__(self, other):
        equality = False
        if issubclass(self.__class__, SpectrumBase) and issubclass(other.__class__, SpectrumBase):
            equality = (
                str(self.scan_number) == str(other.scan_number) and
                self.file_id == other.file_id
            )
        return equality

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return "Base Spectrum (Scan number: {base.scan_number})".format(base=self)

class Cluster(SpectrumBase):

    """Class for individual cluster."""

    def __init__(self, scan_number, file_id, precursor_mz, charge, spectra,
                 purity, consensus_peptide, sqs = 0, ms_cluster_charge = None,
                 spec_prob = None, p_value = None, fdr = None, pep_fdr = None):
        self.spectra = spectra
        self.purity = purity
        self.consensus_peptide = consensus_peptide
        self.sqs = 0
        self.spectrum_weights = self.set_spectrum_weights()
        self.ms_cluster_charge = ms_cluster_charge
        self.spec_prob = spec_prob
        self.p_value = p_value
        self.fdr = fdr
        self.pep_fdr = pep_fdr
        super(Cluster, self).__init__(scan_number, file_id, precursor_mz, charge)
    def avg_spectra_sqs(self):
            output = 0.0
            all_sqs = [
                spectrum.sqs
                for spectrum in self.spectra
            ]
            try:
                 output = sum(all_sqs)/float(len(all_sqs))
            except:
                pass
            return output
    def set_spectrum_weights(self):
        return [
            float(spectrum.precursor_mz)
            for spectrum in self.spectra
        ]
    def guess_charge(self):
        return [
            spectrum.charge
            for spectrum in self.spectra
        ][0]


class Spectrum(SpectrumBase):

    """Class for individual spectrum (or cluster solely for SQS)."""

    def __init__(self, scan_number, file_id, precursor_mz, charge, similarity,
                 p_value, sqs, peptide, cluster_number = None, spec_prob = None,
                 fdr = None, pep_fdr = None):
        self.similarity = similarity
        self.p_value = p_value
        self.sqs = sqs
        self.peptide = peptide
        self.cluster_number = cluster_number
        self.spec_prob = spec_prob
        self.fdr = fdr
        self.pep_fdr = pep_fdr
        super(Spectrum, self).__init__(scan_number, file_id, precursor_mz, charge)
    def output_pair(self):
        return str(self.file_id) + ":" + str(self.scan_number) + ":" + str(self.sqs) + ":" + str(self.charge) + ":" + str(self.precursor_mz)

class DatabaseSpectrum(SpectrumBase):

    """Class for database results."""

    def __init__(self, scan_number, file_id, precursor_mz, charge, peptide, is_cluster):
        self.peptide = peptide
        self.is_cluster = is_cluster
        super(DatabaseSpectrum, self).__init__(scan_number, file_id, precursor_mz, charge)

class Purity(object):

    """Class for purity measures"""

    def __init__(self, representative_spectrum, purity_incl_undentified,
                 purity_no_undentified, peptide_counts):
        self.representative_spectrum = representative_spectrum
        self.purity_incl_undentified = purity_incl_undentified
        self.purity_no_undentified = purity_no_undentified
        self.peptide_counts = peptide_counts

    def distinct_peptides(self):
        return len(peptide_counts.keys())

    def __str__(self):
        return "%.1f" % self.purity_incl_undentified + ": " + str(self.representative_spectrum)

##Core functions

def countby(func, items):
    sorted_items = sorted(items)
    return dict(
        (key,len(list(value)))
        for key, value in groupby(sorted_items, func)
    )

def safe_ratio(item1, item2):

    """
    Mainly for calculating purity,
    sometimes we want to ignore a ZeroDivisionError
    """

    ratio = 0
    try:
        ratio = (float(item1)/float(item2))
    except ZeroDivisionError:
        pass
    return ratio

def swap(dictionary):
    return dict(
        (dictionary[key],key)
        for key in dictionary
    )

def parse_xml_file(input_file):
    """From Ming, for params"""

    key_value_pairs = {}
    with open(input_file, 'r') as input_file_lines:
        for line in input_file_lines:
            new_line = line.rstrip().replace("<parameter name=\"", "")
            #new_line = new_line.replace("\">", "=")
            new_line = new_line.replace("</parameter>", "")

            splits = new_line.split("\">")
            #print splits
            if(len(splits) != 2):
                continue


            if(splits[0] in key_value_pairs.keys()):
              key_value_pairs[splits[0]].append(splits[1])
            else:
              key_value_pairs[splits[0]] = []
              key_value_pairs[splits[0]].append(splits[1])
    return key_value_pairs

def upload_file_mapping(params = None, params_xml = None):
    try:
        upload_file_mapping_list = params['upload_file_mapping']
    except:
        try:
            upload_file_mapping_list = parse_xml_file(params_xml)['upload_file_mapping']
        except:
            raise ValueError("Neither xml file or parameters provided")
    mapping = {}
    for item in upload_file_mapping_list:
        file_map = item.split("|")
        mapping[os.path.split(file_map[1])[1]] = file_map[0]
    return mapping

def datalist_to_map(datalist):
    with open(datalist, 'r') as datalist_items:
        return dict(
            (ind, os.path.split(filename)[1].replace("\n",""))
            for ind, filename in enumerate(datalist_items)
        )

#
# def id_sqs(database_spectra, spectrum, bucket_size):
#     output = 0
#     bucket = int(float(spectrum.precursor_mz)) - int(float(spectrum.precursor_mz)) % bucket_size
#     lower = []
#     try:
#         lower = database_spectra[bucket - bucket_size]
#     except:
#         lower = []
#     try:
#         db = [
#             database.sqs
#             for database in database_spectra[bucket] + lower
#             if database.scan_number == spectrum.scan_number.replace("out_0_0.","")
#         ]
#         try:
#             output = db[0]
#         except IndexError:
#             pass
#     except KeyError:
#         pass
#     return output
#
#
# def id_peptide(database_spectra, spectrum, bucket_size):
#     peptide = None
#     bucket = int(float(spectrum.precursor_mz)) - int(float(spectrum.precursor_mz)) % bucket_size
#     lower = []
#     try:
#         lower = database_spectra[bucket - bucket_size]
#     except:
#         lower = []
#     try:
#         db = [
#             database.peptide
#             for database in database_spectra[bucket] + lower
#             if database == spectrum
#         ]
#         try:
#             peptide = db[0]
#         except IndexError:
#             pass
#     except KeyError:
#         pass
#     return peptide
#
def purity_from_spectrum_list(spectra):
    only_id_spectra = [
        spectrum
        for spectrum in spectra
        if spectrum.peptide is not None
    ]
    total_spectra_incl_undentified = len(spectra)
    total_spectra_no_undentified = len(only_id_spectra)
    peptide_counts = countby(lambda x: x.peptide, only_id_spectra)
    grouped_spectra = list(countby(lambda x: x.peptide, only_id_spectra).items())
    representative_spectrum = ""
    number = 0.0
    if total_spectra_no_undentified != 0:
        representative_spectrum, number = (sorted(grouped_spectra, key = lambda x: x[1]))[-1]
    return Purity(
        representative_spectrum,
        safe_ratio(number, total_spectra_incl_undentified),
        0 if number == 0 else base.safe_ratio(number, total_spectra_no_undentified),
        peptide_counts
    )

def purity_from_peptide_list(peptides):

    """Given a list of peptides, whats the purity"""

    peptides_to_compare = [
        (strip_peptide(peptide),peptide)
        for peptide in peptides
        if peptide is not None
    ]
    total_spectra_incl_undentified = len(peptides)
    total_spectra_no_undentified = len(peptides_to_compare)
    peptide_counts = countby(lambda x: x, peptides_to_compare)
    grouped_spectra = list(map(lambda x: unpack_peptide(x), countby(lambda x: x[0], peptides_to_compare).items()))
    representative_spectrum = ""
    number = 0
    if total_spectra_no_undentified != 0:
        representative_spectrum, number = (sorted(grouped_spectra, key = lambda x: x[1]))[-1]
    return Purity(
        representative_spectrum,
        safe_ratio(number, total_spectra_incl_undentified),
        0 if number == 0 else safe_ratio(number, total_spectra_no_undentified),
        peptide_counts
    )

def strip_peptide(peptide):

    """Return only the letters in a peptide"""

    return (''.join(list(filter(lambda x: x.isalpha(), peptide))),peptide)

def unpack_peptide(tupled_tuple):
    peptide_search = tupled_tuple[0]
    count = tupled_tuple[1]
    return (peptide_search[1],count)

# def database_by_weight(database_spectra, bucket_size):
#     database = {}
#     for spectrum in database_spectra:
#         bucket = int(float(spectrum.precursor_mz)) - int(float(spectrum.precursor_mz)) % bucket_size
#         try:
#             database[bucket].append(spectrum)
#         except KeyError:
#             database[bucket] = [spectrum]
#     return database

# def read_tsv_database(filename, cluster = False, bucket_size = 1):
#     """Read the database output file -- not the clustering output."""
#     database_spectra = []
#     with open(filename) as tsv:
#         next(tsv)
#         for line in tsv:
#             split_line = line.split("\t")
#             scan = split_line[19]
#             if cluster:
#                 scan = int(split_line[47])
#             database_spectra.append(DatabaseSpectrum(
#                 scan_number=scan,
#                 file_id=split_line[29],
#                 precursor_mz=split_line[36],
#                 charge=split_line[9],
#                 peptide=split_line[0],
#                 is_cluster=cluster
#             ))
#     return database_by_weight(database_spectra, bucket_size)

def median_score(spectra):

    """Until ClusterLike"""

    weights = [
        float(spectrum.precursor_mz)
        for spectrum in spectra
    ]
    return round(statistics.median(weights),3)

def mix_score(spectra):

    """Until ClusterLike"""

    weights = [
        float(spectrum.precursor_mz)
        for spectrum in spectra
    ]
    max_weight = max(weights)
    min_weight = min(weights)
    return round(max_weight - min_weight,3)

def assign_sqs_from_spectra(clusters, cluster_sqs):
    for cluster in clusters:
        try:
            clusters[cluster].sqs = cluster_sqs[clusters[cluster].scan_number].sqs
        except:
            clusters[cluster].sqs = 0
    return clusters

def read_cluster_tsv(cluster_folder, database_spectra, title_prefix,
             database_clusters, data_map = None, bucket_size = 1, only_sqs = False):
    """Read the cluster output file -- not the database output."""
    clusters = {}
    for filename in glob.glob(cluster_folder + "/*.clust"):
        with open(filename) as tsv:
            current_cluster = None
            current_spectra = []
            line_number = 0
            for tsv_line in tsv:
                if tsv_line == '\n':
                    current_cluster.purity = 0 #purity_from_spectrum_list(current_spectra)
                    current_cluster.spectra = current_spectra
                    current_cluster.consensus_peptide = "" #id_peptide(database_clusters, current_cluster, bucket_size)
                    clusters.update({current_cluster.scan_number: current_cluster})
                    current_spectra = []
                else:
                    split_tsv_line = tsv_line.replace('\n', '').split('\t')
                    if len(split_tsv_line) == 4:
                        current_cluster = read_cluster_tsv_line(split_tsv_line, title_prefix, filename)
                    elif len(split_tsv_line) == 7 or len(split_tsv_line) == 8:
                        current_spectra.append(read_spectrum_tsv_line(split_tsv_line, database_spectra, data_map, only_sqs = only_sqs))
                    else:
                        raise ValueError(
                            'Line %d is neither a cluster nor a spectrum\n%s' % (line_number,
                                                                                 tsv_line)
                        )
                line_number += 1
    return clusters

def read_cluster_tsv_line(tsv_line_array, title_prefix, filename):
    return Cluster(
        file_id=filename.replace(".clust",".mgf").replace("clust/",""),
        scan_number=tsv_line_array[0].replace(title_prefix,""),
        precursor_mz=tsv_line_array[2],
        charge=tsv_line_array[3],
        spectra=[],
        purity=None,
        consensus_peptide=None
    )

def read_spectrum_tsv_line(tsv_line_array, database_spectra, data_map, only_sqs, bucket_size = 1):
    spectrum = Spectrum(
        scan_number=tsv_line_array[2],
        file_id=0 if only_sqs else data_map[int(tsv_line_array[1])],
        precursor_mz=tsv_line_array[3],
        charge=tsv_line_array[6],
        similarity=tsv_line_array[4],
        p_value=tsv_line_array[5],
        sqs=0,
        peptide=None
    )
    sqs = 0
    try:
        sqs = float(tsv_line_array[7])
    except IndexError:
        pass
    spectrum.sqs = sqs
    spectrum.peptide = None
    return spectrum

def as_spectra(cluster_folder, database_spectra, title_prefix,
         database_clusters, bucket_size = 1):

    clusters = read_cluster_tsv(cluster_folder, database_spectra, title_prefix,
             database_clusters, bucket_size, only_sqs = True)



    flat_clusters = [
        clusters[key]
        for key in clusters
    ]

    cluster_spectra = dict(
        (spectrum.scan_number, spectrum)
        for cluster in flat_clusters
        for spectrum in cluster.spectra
    )

    return cluster_spectra

# def init_from_pickle(cluster_pickle):
#     clusters = []
#     with open(cluster_pickle,"rb") as f:
#         clusters = pickle.load(f)
#     print("Successfully loaded " + str(len(clusters)) + " clusters from " + cluster_pickle)
#     return cluster_session.ClusterSession(clusters)

def initialize_unidentified_clusters(cluster_folder, clustered_clusters_folder, title_prefix, data_map):
    clusters_as_spectra = as_spectra(clustered_clusters_folder, {}, title_prefix, {})
    clusters = read_cluster_tsv(cluster_folder, {}, title_prefix, {}, data_map)
    return assign_sqs_from_spectra(clusters, clusters_as_spectra)

def find_cluster(clusters, cluster_id):
    found_cluster = None
    for cluster in clusters:
        if int(cluster.scan_number) == cluster_id:
            found_cluster = cluster
            break
    return found_cluster

def initialize_from_proteosafe(clusters_file, spectra_file):
    clusters = {}
    with open(clusters_file, 'r') as cluster_lines:
        cluster_lines.readline()
        for cluster_line in cluster_lines:
            clust_split = cluster_line.split('\t')
            scan_number = clust_split[2]
            new_cluster = Cluster(
                scan_number = scan_number,
                file_id = clust_split[0],
                precursor_mz = clust_split[3],
                charge = int(clust_split[6].replace('+','')),
                spectra = [],
                purity = None,
                consensus_peptide = clust_split[7],
                ms_cluster_charge = clust_split[21],
                spec_prob = clust_split[11],
                p_value = clust_split[12],
                fdr = clust_split[13],
                pep_fdr = clust_split[14]
            )
            clusters[scan_number] = new_cluster
    with open(spectra_file, 'r') as spectrum_lines:
        spectrum_lines.readline()
        for spectrum_line in spectrum_lines:
            spec_split = spectrum_line.split('\t')
            cluster_number = spec_split[15]
            new_spectrum = Spectrum(
                scan_number = spec_split[2],
                file_id = spec_split[0],
                precursor_mz = spec_split[4],
                charge = int(spec_split[6]),
                similarity = None,
                p_value = spec_split[12],
                sqs = spec_split[16],
                peptide = spec_split[7],
                cluster_number = cluster_number,
                spec_prob = spec_split[11],
                fdr = spec_split[13],
                pep_fdr = spec_split[14]
            )
            clusters[cluster_number].spectra.append(new_spectrum)
    return clusters
#
# def load_to_pickle(cluster_folder, clusters_clusters_folder, identified_clusters,
#                    identified_spectra, output_pickle, title_prefix):
#     print("Loading clusters...this might take a moment.")
#     database_spectra = database.read_tsv(identified_spectra)
#     database_clusters = database.read_tsv(identified_clusters, cluster = True)
#     clusters_as_spectra = as_spectra(clusters_clusters_folder, database_spectra, title_prefix, database_clusters)
#     clusters = read_tsv(cluster_folder, database_spectra, title_prefix, database_clusters)
#
#     new_clusters = assign_sqs_from_spectra(clusters, clusters_as_spectra)
#
#     with open(output_pickle,"wb") as f:
#         pickle.dump(new_clusters, f)
#
#     print("Saved " + str(len(clusters)) + " clusters to " + output_pickle)
