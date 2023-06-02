import numpy as np
import networkx as nx
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

'''
if len(sys.argv) < 2:
    raise Exception('network information missing')
elif len(sys.argv) < 3:
    raise Exception('disease genes missing')
networkFile = sys.argv[1]
network = np.loadtxt(networkFile, delimiter=' ', skiprows=1, dtype=str)

targetsFile = sys.argv[2]
targets = np.loadtxt(targetsFile, dtype=str)

if len(sys.argv) > 3:
    candidatesFile = sys.argv[3]
    candidates = np.loadtxt(candidatesFile, dtype=str)
else:
    candidatesFile = None
    candidates = np.array(list(set([gene[0].split('.')[1] for gene in network])))
if len(sys.argv) > 4:
    threshold: float = float(sys.argv[4])
else:
    threshold = 0
'''

def shortest_path(network, source, target, threshold):
    path_list = [[source]]
    path_index = 0
    previous_nodes = {source}
    if source == target:
        return 0

    while path_index < len(path_list):
        current_path = path_list[path_index]
        last_node = current_path[-1]
        next_nodes = [network[i][1].split('.')[1]
                            for i in range(len(network))
                            if network[i][0].split('.')[1] == last_node
                            and int(network[i][2]) > threshold]
        if target in next_nodes:
            current_path.append(target)
            return len(current_path)-1

        for next_node in next_nodes:
            if not next_node in previous_nodes:
                new_path = current_path[:]
                new_path.append(next_node)
                path_list.append(new_path)
                previous_nodes.add(next_node)
        path_index += 1
    return None

def load_network(filepath):
    file_object = open(filepath, mode="r")
    network_object = nx.DiGraph()
    for line in file_object:
        if line.startswith("#") or line == "\n" or line.upper().startswith("P"):
            continue
        line = line.upper().split()
        network_object.add_edge(line[0].split('.')[1], line[1].split('.')[1])
    return network_object

def plots(result, net):
    plt.figure()
    #sns.set()
    sns.barplot(x=list(result.keys()), y=list(result.values()))
    plt.xticks(rotation=80)
    plt.title('Candidate Genes Scores')
    plt.ylabel('Score')
    plt.xlabel('Genes')
    plt.savefig('results_image.pdf')

    plt.figure()
    net_undirected = net.to_undirected()
    nx.draw(net_undirected.subgraph(max(nx.connected_components(net_undirected))), with_labels=True, node_size = 0.5, font_size=0.5, width=0.2)
    plt.savefig('network.pdf')

def genePred(networkFile, seedFile, candidateFile = None, threshold = 0, graphical = True):
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
    :param graphical: bool; default = True; if True, graphical output which compares
                        the scores of the individual candidate genes will be produced
    :return: a results.txt file containing the ranked list of the candidate genes,
                        including the scores of their vicinty to the seed genes;
                        a details.txt file with more detailed information on the results.

    :example: genePred('hsNetwork.txt', 'seedGenes.txt', 'candidateGenes.txt', threshold = 300)
    '''

    network = np.loadtxt(networkFile, delimiter=' ', skiprows=1, dtype=str)
    targets = np.loadtxt(seedFile, dtype=str)
    if candidateFile != None:
        candidates = np.loadtxt(candidateFile, dtype=str)
    else:
        candidates = np.array(list(set([gene[0].split('.')[1] for gene in network])))

    results = {}
    with open('details.txt', 'w') as df:
        df.write(f'Scientific Programming Project 2\nPackage for disease gene prediction based on an protein-protein interaction network\nJob datetime: {dt_string}')
        df.write(f'\nParameters used:\n\tNetwork file: {networkFile}\n\tTarget genes file: {seedFile}\n\tCandidate genes file: {candidateFile}\n\tcombined_score = {threshold}\n\n')
        for c in candidates:
            df.write(f'\n\nCandidate gene\t{c}')
            results[c] = [0, 0]
            print(f'Analyzing candidate gene {c}')
            for t in targets:
                print(f'Computing distance to seed gene {t}:', end=' ')
                d = shortest_path(network, c, t, threshold)
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

    if graphical == True:
        net = load_network(networkFile)
        plots(result, net)
