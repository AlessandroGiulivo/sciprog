import sys
import numpy as np
import networkx as nx
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")


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

def shortest_path(source, target, threshold=0):
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


results = {}
with open('details.txt', 'w') as df:
    df.write(f'Scientific Programming Project 2\nPackage for disease gene prediction based on an protein-protein interaction network\nJob datetime: {dt_string}')
    df.write(f'\nParameters used:\n\tNetwork file: {networkFile}\n\tTarget genes file: {targetsFile}\n\tCandidate genes file: {candidatesFile}\n\tcombined_score = {threshold}\n\n')
    for c in candidates:
        df.write(f'\n\nCandidate gene\t{c}')
        results[c] = [0, 0]
        for t in targets:
            print(f'candidate {c}\ttarget {t}')
            d = shortest_path(c, t, threshold)
            if type(d) == int:
                results[c][0] += d
                results[c][1] += 1
            df.write(f'\n\tDistance to target gene\t{t}\t{d}')
        df.write(f'\nTotal distance from targets:\t{results[c][0]}\nNumber of targets connected to candidate:\t{results[c][1]}')
        try:
            df.write(f'\nScore:\t{results[c][0]/results[c][1]}')
        except:
            continue


results = {k: v[0]/v[1] for k, v in sorted(results.items(), key=lambda item: item[1][0]/item[1][1], reverse=True) if v[1] != 0}


with open('results.txt', 'w') as rf:
    for k, v in results.items():
        rf.write(f'{k}\t{round(v, 3)}\n')

plt.figure()
#sns.set()
sns.barplot(x=list(results.keys()), y=list(results.values()))
plt.savefig('results_image.png')

def load_network(filepath):
    file_object = open(filepath, mode="r")
    network_object = nx.DiGraph()

    for line in file_object:
        if line.startswith("#") or line == "\n" or line.startswith("protein"):
            continue
        line = line.upper().split()
        network_object.add_edge(line[0], line[1])
    return network_object


net = load_network(networkFile)

plt.figure()
nx.draw_networkx(net, node_size = 20, font_size=12)
plt.savefig('network.png')





