from pathlib import Path

import csv
import scipy
import itertools
import collections

import newick

import scipy.stats

c = collections.defaultdict(lambda: scipy.zeros(6, int))

with Path("../beastling/indexes.log").open() as indices:
    for line in csv.DictReader(indices, delimiter="\t"):
        for id, index in line.items():
            if id == "Sample":
                continue
            _, concept = id.rsplit(":", 1)
            c[concept][int(index)] += 1

def distance(d1, d2):
    return (scipy.stats.chisquare(d1+1, d2+1).statistic)

nodes = {newick.Node(n): d for n, d in c.items()}
old_nodes = {}
distances = {(n1, n2): distance(d1, d2)
             for (n1, d1), (n2, d2) in itertools.combinations(nodes.items(), 2)}

while len(nodes) > 1:
    argmin = min(distances, key=distances.get)
    d0 = nodes.pop(argmin[0])
    d1 = nodes.pop(argmin[1])
    old_nodes[argmin[0]] = d0
    old_nodes[argmin[1]] = d1
    d = distances.pop(argmin)
    argmin[0].length = d/2
    argmin[1].length = d/2
    n = newick.Node(None)
    n.add_descendant(argmin[0])
    n.add_descendant(argmin[1])
    for i in list(distances):
        if argmin[0] == i[0] or argmin[1] == i[0]:
            distances.setdefault((n, i[1]), 0)
            distances[n, i[1]] += distances.pop(i) / 2
        if argmin[0] == i[1] or argmin[1] == i[1]:
            distances.setdefault((n, i[0]), 0)
            distances[n, i[0]] += distances.pop(i) / 2
    nodes[n] = d0 + d1

print(list(nodes)[0].newick, end=";\n", file=Path("../beastling/clusters.newick").open("w"))

with Path("../beastling/clusters.csv").open("w") as ocsv:
    for i in list(nodes)[0].get_leaves():
        print(i.name, *old_nodes[i], sep="\t", file=ocsv)
