import sys
import itertools

import ete3

taxa = [
    "A-Atimelang",
    "A-Petleng",
    "A-Takalelang",
    "A-Ulaga",
    "Kafoa",
    "Kamang",
    "Kiraman",
    "Klon",
]

def ancestors(node, tree):
    a = {node}

    r = tree.get_leaves_by_name(node.name)
    for reticulation in r:
        if reticulation in a:
            continue
        a |= ancestors(reticulation, tree)
    if node.up:
        return a | ancestors(node.up, tree)
    else:
        return a

aggregated_distances = {}
for line in sys.stdin:
    t = ete3.Tree(line, format=1)
    for l1, l2 in itertools.combinations(taxa, 2):
        l1, = t.get_leaves_by_name(l1)
        l2, = t.get_leaves_by_name(l2)
        shared_ancestors = sorted(
            ancestors(l1, t) & ancestors(l2, t),
            key=t.get_distance)
        mrca = shared_ancestors[-1]
        aggregated_distances.setdefault((l1.name, l2.name), []).append(
            1 - 2 * t.get_distance(mrca) / (t.get_distance(l1) + t.get_distance(l2)))

from matplotlib import pyplot as plt
ds = sum(aggregated_distances.values(), [])
min_d = min(ds)
max_d = max(ds)

for i in range(len(taxa)):
    for j in range(i):
        plt.subplot(len(taxa) - 1, len(taxa), j * len(taxa) + i)
        try:
            dist = aggregated_distances[taxa[i], taxa[j]]
        except KeyError:
            dist = aggregated_distances[taxa[j], taxa[i]]
        plt.hist(dist, range=(min_d, max_d), bins=20)
plt.show()
