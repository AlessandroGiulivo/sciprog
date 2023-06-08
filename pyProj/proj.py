import numpy as np
import networkx as nx
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

def load_network(filepath, w):
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

def plots(result, net, nodes, draw_network):
    plt.figure()
    sns.barplot(y=list(result.keys()), x=list(result.values()), orient = 'h')
    plt.title('Candidate Genes Scores')
    plt.ylabel('Genes')
    plt.xlabel('Score')
    plt.savefig('resultBarplot.pdf', bbox_inches='tight')

    if draw_network == True:
        plt.figure()
        nx.draw_networkx(net, nodelist = nodes, with_labels=True, node_size = 5, font_size=3, width=0.2)
        plt.title('Connections between Seed Genes and Candidate Genes')
        plt.savefig('resultNetwork.pdf')

def genePred(networkFile, seedFile, candidateFile = None, threshold = 0, barplot=True, drawNetwork = True):
    '''
    The function ranks candidate disease genes according to their vicinity to
    the seed genes on a PPI network of a single species.

    :param networkFile: a weighted protein-protein interaction network
                        from the STRING database of a single species
                        (e.g. Homo Sapiens) as a .txt file
    :param seedFile: a set of known disease-related “seed” (or reference) genes
    :param candidateFile: optional - a set of “candidate” disease genes
                        (if not specified, all non-seed genes will be taken as candidates)
    :param threshold: int; default = 0; minimum cutoff for the interactions to be considered
    :param barplot: bool; default = True; if True, graphical output file which compares
                        the scores of the individual candidate genes will be produced
    :param drawNetwork: bool; default = True; if True, graphical output file showing how the
                        interesting genes are connected will be produced.
                        Not recommended if parameter 'candidateFile' is not specified.
                        
    :return a results.txt file containing the ranked list of the candidate genes,
                        including the scores of their vicinty to the seed genes;
                        a details.txt file with more detailed information on the results.

    :example from proj import *
             genePred('hsNetwork.txt', 'seedGenes.txt', 'candidateGenes.txt', threshold = 300)


    :note Gene Identifiers used in 'networkFile', 'seedFile', 'candidateFile' should be of the same type.
    '''

    net = np.loadtxt(networkFile, delimiter=' ', skiprows=1, dtype=str)
    network = load_network(networkFile, threshold)
    targets = np.loadtxt(seedFile, dtype=str)
    if candidateFile != None:
        candidates = np.loadtxt(candidateFile, dtype=str)
    else:
        candidates = np.array([gene for gene in list(network.nodes)])

    results = {}
    with open('resultsDetailed.txt', 'w') as df:
        df.write(f'Scientific Programming Project 2\nPackage for disease gene prediction based on an protein-protein interaction network\nJob datetime: {dt_string}')
        df.write(f'\nParameters used:\n\tNetwork file: {networkFile}\n\tTarget genes file: {seedFile}\n\tCandidate genes file: {candidateFile}\n\tcombined_score = {threshold}\n\n')
        nodesList = []
        for c in candidates:
            df.write(f'\n\nCandidate gene\t{c}')
            results[c] = [0, 0]
            print(f'Analyzing candidate gene {c}')
            for t in targets:
                print(f'Computing distance to seed gene {t}:', end=' ')
                try:
                    sp = nx.shortest_path(network, source=c, target=t)
                    d = len(sp)
                    nodesList.extend(sp)
                except:
                    d = None
                if type(d) == int:
                    results[c][0] += d
                    results[c][1] += 1
                df.write(f'\n\tDistance to seed gene\t{t}\t{d}')
                if d == None:
                    print('not found')
                else:
                    print(f'{d}')
            df.write(f'\nTotal distance from targets:\t{results[c][0]}\nNumber of targets connected to candidate:\t{results[c][1]}')
            try:
                df.write(f'\nScore:\t{results[c][0]/(results[c][1]**2)}')
            except:
                continue

    result = {}
    for k, v in results.items():
        try:
            result[k] = v[0]/(v[1]**2)
        except:
            continue
    result = {k: v for k, v in sorted(result.items(), key=lambda item: item[1])}
    if result == {}:
        raise Exception('No candidate gene was found to have connections to seed genes')

    with open('results.txt', 'w') as rf:
        for k, v in result.items():
            rf.write(f'{k}\t{round(v, 3)}\n')

    if barplot == True:
        nodes = list(candidates)
        nodes.extend(list(targets))
        net = network.subgraph(set(nodesList))
        plots(result, net, nodes, drawNetwork)
