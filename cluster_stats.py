from . import cluster_base
import matplotlib
from tabulate import tabulate
from collections import OrderedDict
from IPython import display
from IPython.display import HTML

def str_percent(num):
    return "{0:.2f}".format(num*100) + "%"

class TableSession(object):

    def __init__(self, clusters, x, y = None, z = None, f = None, notes = "", plotter = None, proteosafe_id = ""):
        self.clusters = clusters
        self.x = x
        self.y = y
        self.z = z
        self.f = f
        self.notes = notes
        self.plotter = plotter
        self.proteosafe_id = proteosafe_id
    def display(self):
        t = table(self.clusters, self.x, self.y, self.f)
        return HTML(t)
    def cdisplay(self, func, item_set):
        table = table_data(self.clusters, self.x, self.y, lambda x: x)
        num_total = len(item_set)
        new_total = 0
        new_table = [[] for r in range(0,len(table))]
        for r,row in enumerate(table):
            new_table[r] = [[] for c in range(0,len(row))]
            for c,column in enumerate(row):
                if c == 0:
                    new_table[r][c] = column
                else:
                    column_set = func(column)
                    column_total = column_set & item_set
                    new_total += len(column_total)
                    new_table[r][c] = str_percent(new_total/num_total)
                    item_set = item_set - column_set
        return HTML(tabulate(new_table,self.x.table_header(),tablefmt="html"))
    def inspect(self, x_insp, y_insp, n = 5, all_values = False):
        t = table_data(self.clusters, self.x, self.y, lambda x: x)
        values = [
            cluster.scan_number
            for cluster in t[y_insp][x_insp+1]
        ]
        if not all_values:
            return HTML('&nbsp;&nbsp;&nbsp;'.join([cluster_number_url(value,self.proteosafe_id) for value in values[:n]]))
        else:
            return values
    def plot(self):
        x_plots = [self.x.func(cluster) for cluster in self.clusters]
        if self.y:
            y_plots = [self.y.func(cluster) for cluster in self.clusters]
            return self.plotter.plot(x_plots,y_plots)
        else:
            return self.plotter.plot(x_plots)
    def data(self):
        return table_data(self.clusters, self.x, self.y, lambda x: x)

def inspect(table, x, y):
    to_inspect = table[y][x+1]
    return [
        cluster.scan_number
        for cluster in to_inspect
    ]

def cluster_number_url(cluster_number, proteosafe_id):
    cluster_url1 = "http://ccms-dev1.ucsd.edu/ProteoSAFe/result.jsp?task="
    cluster_url2 = "&view=group_by_spectrum_old_clustered#%7B%22table_sort_history%22%3A%22SpecProb_dsc%3BSpecProb_asc%3BSpecProb_dsc%22%2C%22Scan%23_lowerinput%22%3A%22"
    cluster_url3 = "%22%2C%22Scan%23_upperinput%22%3A%22"
    cluster_url4 = "%22%7D"
    return "<a target=\"_blank\" href=\"" + cluster_url1 + proteosafe_id + cluster_url2 + cluster_number + cluster_url3 + cluster_number + cluster_url4 + "\">"+cluster_number+"</a>"

class ClusterSession(object):

    def __init__(self, clusters, plotter, proteosafe_id = "", default_func = len):
        self.clusters = clusters
        self.default_func = default_func
        self.plotter = plotter
        self.proteosafe_id = proteosafe_id
    # def table(self, x, y = None, z = None, f = None):
    #     func = f if f else self.default_func
    #     t = table(self.clusters, x, y, func)
    #     print(t)
    # def table_data(self, x, y = None, z = None, f = None):
    #     func = lambda x: x
    #     return table_data(self.clusters, x, y, func)
    def table(self, x, y = None, z = None, f = None):
        func = f if f else self.default_func
        return TableSession(self.clusters, x, y, z, func, plotter = self.plotter, proteosafe_id = self.proteosafe_id)


class Header(object):

    def __init__(self, title, func, header_items, continuous = False):
        self.title = title
        self.func = func
        self.header_items = header_items
        self.continuous = continuous
    def to_empty_dict(self):
        return OrderedDict([
                            (str(header_item),[])
                            for header_item in sorted(self.header_items, key = lambda x: x.lower + x.upper)
                            ])
    def table_header(self):
        return [""] + [
                       str(header_item)
                       for i, header_item in enumerate(sorted(self.header_items, key = lambda x: x.lower + x.upper))
                       ]


class HeaderItem(object):

    def __init__(self, lower = -1e309, upper = 1e309, display = None, incl_upper = False):
        self.display = display
        self.lower = lower
        self.upper = upper
        self.incl_upper = incl_upper
    def __str__(self):
        if self.display is not None:
            return self.display
        else:
            if self.lower == self.upper:
                return str(self.lower)
            elif self.lower == -1e309:
                return "- " + str(self.upper)
            elif self.upper == 1e309:
                return str(self.lower) + " +"
            else:
                return str(self.lower) + " - " + str(self.upper)

def group_in_header(items, header, reduce_func = lambda x: x):
    group = header.to_empty_dict()
    for item in items:
        transformed_item = header.func(item)
        for header_item in header.header_items:
            #Logic to decipher between continuous and not
            if transformed_item is not None and (
             transformed_item >= header_item.lower and
             (transformed_item < header_item.upper or
              (transformed_item == header_item.upper and (
                                                          not header.continuous or header_item.incl_upper)
               ))):
               group[str(header_item)].append(item)
               break
    for key in group:
        group[key] = reduce_func(group[key])
    return group

def table_1d_data(to_sort, header_x, func):
    return group_in_header(to_sort, header_x, func).items()

def table_2d_data(to_sort, header_x, header_y, func):
    group_y = group_in_header(to_sort, header_y)
    # for group, stuff in group_y.items():
    #     group_x = group_in_header(stuff, header_x)
    #     print([group] + list(map(lambda x: str(x), group_x.values())))
    return [
            [str(group[0])] +
            list(map(func, group_in_header(group[1], header_x).values()))
            for i, group in enumerate(group_y.items())
            ]

def table_1d(to_sort, header_x, func):
    return tabulate(table_1d_data(to_sort, header_x, func),[header_x.title,""],tablefmt="html")

# def table_2d(to_sort, header_x, header_y, func):
#     group_y = group_in_header(to_sort, header_y)
#     # for group, stuff in group_y.items():
#     #     group_x = group_in_header(stuff, header_x)
#     #     print([group] + list(map(lambda x: str(x), group_x.values())))
#     return "X: " + header_x.title + "\nY: " + header_y.title + "\n" + tabulate([
#         [group[0]] +
#         list(map(func, group_in_header(group[1], header_x).values()))
#         for group in group_y.items()
#     ],header_x.table_header(),tablefmt="grid")

def table_2d(to_sort, header_x, header_y, func):
    return "X: " + header_x.title + "\nY: " + header_y.title + "\n" + tabulate(table_2d_data(to_sort, header_x, header_y, func),header_x.table_header(),tablefmt="html")

def table(to_sort, header_x, header_y = None, func = lambda x: x):
    if header_y:
        return table_2d(to_sort, header_x, header_y, func)
    else:
        return table_1d(to_sort, header_x, func)
#print(group_in_header(group_in_header(to_sort,header_x),header_y))

def table_data(to_sort, header_x, header_y = None, func = lambda x: x):
    if header_y:
        return table_2d_data(to_sort, header_x, header_y, func)
    else:
        return table_1d_data(to_sort, header_x, func)
#print(group_in_header(group_in_header(to_sort,header_x),header_y))


def inspect(table, x, y):
    to_inspect = table[y][x+1]
    return [
            cluster.scan_number
            for cluster in to_inspect
            ]

def pretty_dict(dict_to_pretty, func = lambda x: x):
    return "\n".join([
                      key.replace("_"," ").title() + ": " + str(func(dict_to_pretty[key]))
                      for key in dict_to_pretty
                      ])


def average_sqs_single(cluster):
    output = 0.0
    all_sqs = [
               spectrum.sqs
               for spectrum in cluster.spectra
               ]
    try:
       output = sum(all_sqs)/float(len(all_sqs))
    except:
       pass
    return output




def spectra(clusters):
    return len([
                spectrum
                for cluster in clusters
                for spectrum in cluster.spectra
                ])

def clustered_identified(clusters):
    return len({
               cluster.consensus_peptide
               for cluster in clusters
               if cluster.consensus_peptide is not None
               })

def unclustered_identified(clusters):
    return len({
               cluster.purity.representative_spectrum
               for cluster in clusters
               if cluster.purity.representative_spectrum is not None
               })


def net_pep(clusters):
    return clustered_identified(clusters) - unclustered_identified(clusters)

one_to_fifty = [
    HeaderItem(1,1),
    HeaderItem(2,2),
    HeaderItem(3,3),
    HeaderItem(4,4),
    HeaderItem(5,10),
    HeaderItem(10,20),
    HeaderItem(20,50),
    HeaderItem(50)
    ]

cluster_size = Header(
                      "Cluster Size",
                      lambda x: len(x.spectra),
                      [
                       HeaderItem(1,1),
                       HeaderItem(2,2),
                       HeaderItem(3,3),
                       HeaderItem(4,4),
                       HeaderItem(5,5),
                       HeaderItem(6,6),
                       HeaderItem(7,7),
                       HeaderItem(8,8),
                       HeaderItem(9,9),
                       HeaderItem(10,10),
                       HeaderItem(11,15),
                       HeaderItem(16,20),
                       HeaderItem(21,30),
                       HeaderItem(31,40),
                       HeaderItem(41,50),
                       HeaderItem(51,100),
                       HeaderItem(101,200),
                       HeaderItem(201,500),
                       HeaderItem(501)
                       ]
                      )

coarse_0_1 = [
              HeaderItem(0.0,0.3,"Low (0.0-0.3)"),
              HeaderItem(0.3,0.7,"Medium (0.3-0.7)"),
              HeaderItem(0.7,1.0,"High (0.7-1.0)", incl_upper = True)
              ]

fine_0_1 = [
            HeaderItem(0.0,0.1),
            HeaderItem(0.1,0.2),
            HeaderItem(0.2,0.3),
            HeaderItem(0.3,0.4),
            HeaderItem(0.4,0.5),
            HeaderItem(0.5,0.6),
            HeaderItem(0.6,0.7),
            HeaderItem(0.7,0.8),
            HeaderItem(0.8,0.9),
            HeaderItem(0.9,1.0),
            HeaderItem(1.0,1.0, incl_upper = True)
            ]

identified = [
              HeaderItem(0,0,"No"),
              HeaderItem(1,display="Yes")
              ]

consensus = Header(
                   "Cluster Consensous",
                   lambda x: 1 if x.consensus_peptide else 0,
                   identified
                   )

purity_no_undentified_fine = Header(
                                    "Purity (No unidentified)",
                                    lambda x: x.purity.purity_no_undentified,
                                    fine_0_1,
                                    continuous = True
                                    )

purity_no_undentified_coarse = Header(
                                      "Purity (No unidentified)",
                                      lambda x: x.purity.purity_no_undentified,
                                      coarse_0_1,
                                      continuous = True
                                      )


purity_including_undentified_fine = Header(
                                           "Purity (Incl unidentified)",
                                           lambda x: x.purity.purity_incl_undentified,
                                           fine_0_1,
                                           continuous = True
                                           )

purity_including_undentified_coarse = Header(
                                             "Purity (Incl unidentified)",
                                             lambda x: x.purity.purity_incl_undentified,
                                             coarse_0_1,
                                             continuous = True
                                             )

purity = purity_no_undentified_fine
purity_coarse = purity_no_undentified_coarse
purity_fine = purity_no_undentified_fine

cluster_sqs_coarse = Header(
                            "Cluster Quality",
                            lambda x: x.sqs,
                            coarse_0_1,
                            continuous = True
                            )

cluster_sqs_fine = Header(
                          "Cluster Quality",
                          lambda x: x.sqs,
                          fine_0_1,
                          continuous = True
                          )

spectrum_sqs_coarse = Header(
                             "Average Spectrum Quality",
                             average_sqs_single,
                             coarse_0_1,
                             continuous = True
                             )

spectrum_sqs_fine = Header(
                           "Average Spectrum Quality",
                           average_sqs_single,
                           fine_0_1,
                           continuous = True
                           )

split = Header(
        "Spectral split",
        lambda x: x.of_split,
        one_to_fifty,
        continuous = False
)

cluster_sqs = cluster_sqs_fine
spectrum_sqs = spectrum_sqs_fine
