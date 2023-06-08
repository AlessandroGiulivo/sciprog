###                            Alessandro Giulivo
###                                 10940033
###                 Scientific Programming - Prof Rosario M. piro
###                        Python Project (Project 2)
### Package for disease gene prediction based on an protein-protein interaction network

# This python program provides a "genePred" function that takes:
# - a Protein-Protein Interaction network from an individual species as downloaded from the STRING database,
# – a set of known disease-related “seed” (or reference) genes,
# - (optionally) a set of “candidate” disease genes
# and returns as output a ranked list of the candidate genes, including the scores used for ranking them.
# results are saved into multiple "results*" files.

import numpy as np
import networkx as nx
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

def load_network(filepath, w):
    '''
    This function is used by 'genePred' function for loading the
    PPI Network file into a 'networkX object'
    '''
    file_object = open(filepath, mode="r")
    network_object = nx.DiGraph()
    for line in file_object:
        if line.startswith("#") or line == "\n" or line.upper().startswith("P"):
            continue
        line = line.upper().split()
        if int(line[2]) >= w:
            network_object.add_edge(line[0].split('.')[1], line[1].split('.')[1], weight= int(line[2]))
    network_object = network_object.to_undirected()
    return network_object

def genePred(networkFile, seedFile, candidateFile = None, threshold = 0, barplot=True, drawNetwork = True):
    '''
    The function ranks candidate disease genes according to their vicinity to
    the seed genes on a PPI network of a single species.

    :param networkFile: str; a weighted protein-protein interaction network
                        from the STRING database of a single species
                        (e.g. Homo Sapiens) as a .txt file
    :param seedFile:    str; a set of known disease-related “seed” (or reference) genes in a .txt file
    :param candidateFile: str; (optional) a set of “candidate” disease genes in a .txt file
                        (if not specified, all non-seed genes will be taken as candidates)
    :param threshold:   int; default = 0; minimum cutoff for the interactions to be considered
    :param barplot:     bool; default = True; if True, graphical output file which compares
                        the scores of the individual candidate genes will be produced
    :param drawNetwork: bool; default = True; if True, graphical output file showing how the
                        interesting genes are connected will be produced.
                        Not recommended if parameter 'candidateFile' is not specified.
                        
    returns
            - a results.txt file containing the ranked list of the candidate genes,
                including the scores of their vicinty to the seed genes;
            - a resultsDetailed.txt file with more detailed information on the results.
            - (optionally) a resultBarplot.pdf file that graphically compares
                the scores of the individual candidate genes
            - (optionally) a resultNetwork.pdf file that shows how distant or close
                seed genes and candidate genes are within the PPI network.


    example
            from proj import *
            scores = genePred('hsNetwork.txt', 'seedGenes.txt', 'candidateGenes.txt', threshold = 300)

    notes
            Gene Ids used in 'networkFile', 'seedFile', 'candidateFile' should be of the same type.
    '''

    network = load_network(networkFile, threshold)                      # loading network file into a networkx object
    targets = np.loadtxt(seedFile, dtype=str)                           # loading seed genes
    if candidateFile != None:                                           # loading target genes
        candidates = np.loadtxt(candidateFile, dtype=str)
    else:
        candidates = np.array([gene for gene in list(network.nodes)])

    results = {}
    with open('resultsDetailed.txt', 'w') as df:                            # start writing a 'resultsDetailed.txt' file
        df.write(f'Scientific Programming Project 2\n' +
                 f'Package for disease gene prediction based on an protein-protein interaction network\n' +
                 f'Job datetime: {dt_string}')
        df.write(f'\nParameters used:\n' +
                 f'\tNetwork file: {networkFile}\n' +
                 f'\tTarget genes file: {seedFile}\n' +
                 f'\tCandidate genes file: {candidateFile}\n' +
                 f'\tcombined_score = {threshold}\n\n')
        nodesList = [] #list of nodes to draw later in 'resultNetwork.pdf' file
        for c in candidates:                                                  # take each candidate gene
            df.write(f'\n\nCandidate gene\t{c}')
            results[c] = [0, 0] # candidateGene : [total_distance_to_seed_genes,number_of_seed_genes_connected]
            print(f'Analyzing candidate gene {c}')
            for t in targets:                                                 # and find shortest path to each seed gene
                print(f'\tComputing distance to seed gene {t}:', end=' ')
                try:
                    sp = nx.shortest_path(network, source=c, target=t)
                    d = len(sp)
                    nodesList.extend(sp)
                except:
                    d = None
                    nodesList.extend([c, t])
                if type(d) == int:
                    results[c][0] += d
                    results[c][1] += 1
                df.write(f'\n\tDistance to seed gene\t{t}\t{d}')
                if d == None:
                    print('not found')
                else:
                    print(f'{d}')
            df.write(f'\nTotal distance from targets:\t{results[c][0]}\n' +
                     f'Number of targets connected to candidate:\t{results[c][1]}')
            try:
                df.write(f'\nScore:\t{results[c][0]/(results[c][1]**2)}')
            except:
                continue

    result = {}
    for k, v in results.items():                # compute score for each candidate gene such that it is
        try:                                    # directly proportional to its distance to the seed genes and
            result[k] = v[0]/(v[1]**2)          # inversely proportional to the number of seed genes to which
        except:                                 # it was found to be connected
            continue                            # the lower the score, the higher the gene will be in the ranking

    result = {k: v for k, v in sorted(result.items(), key=lambda item: item[1])} # sort the scores in ascending order
    if result == {}:
        return 'No candidate gene was found to have connections to seed genes'

    with open('results.txt', 'w') as rf:                                            # save the results in "results.txt"
        for k, v in result.items():
            rf.write(f'{k}\t{round(v, 3)}\n')

    if barplot:
        plt.figure()                                                                    # generating a barplot
        sns.barplot(y=list(result.keys()), x=list(result.values()), orient='h')         # comparing scores of each
        plt.title('Candidate Genes Scores')                                             # candidate gene
        plt.ylabel('Genes')
        plt.xlabel('Score')
        plt.tick_params(axis='y', which='major', labelsize=25/len(candidates)**(1/2))
        plt.savefig('resultBarplot.pdf', bbox_inches='tight')                           # saving in "resultBarplot.pdf"

    if drawNetwork:
        genes = set(candidates).union(set(targets))
        net = network.subgraph(set(nodesList))
        nodes = [g for g in genes if g in list(net.nodes)] # removing non-connected genes
        net = network.subgraph(set(nodesList))

        pos = nx.spring_layout(net)
        labels = {}
        for node in net.nodes():
            if node in nodes:
                labels[node] = node

        plt.figure()
        nx.draw_networkx(net, pos=pos, nodelist=nodes,                                      # drawing subnetwork between
                         with_labels=False, node_size=2,                                    # containing shortest paths
                         font_size=0.2, width=0.002)                                        # between seed genes
        nx.draw_networkx_labels(net, pos=pos, labels=labels, font_size=1, font_color='r')   # candidate genes
        plt.title('Connections between Seed Genes and Candidate Genes')
        plt.savefig('resultNetwork.pdf')                                            # saving in "resultNetwork.pdf"
    return result