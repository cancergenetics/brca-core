import pandas as pd
import numpy as np
import random
import math
import networkx as nx
import operator
import sys
import pickle
import argparse
from collections import defaultdict
from copy import deepcopy

parser = argparse.ArgumentParser(description='This integrates protein expression with PPI data')
parser.add_argument('-f', action='store', dest='dataset_file', 
    help='The protein expression dataset filename')
parser.add_argument('-n', action='store', dest='dataset_name', 
    help='A shortname for the dataset')
parser.add_argument('-p', action='store', dest='protein_interactions_file', 
    help='The protein interaction dataset name', default="HINT_cocomplex.txt")
parser.add_argument('-d', action='store', dest='degree',type=int, 
    help='Threshold for discarding hubs from the protein interaction network', default=300)
parser.add_argument('-r', action='store', dest='randomisations',type=int, 
    help='Number of randomisations to perform to estimate false discovery rate', default=10)


results = parser.parse_args()
if not results.dataset_file or not results.dataset_name:
    parser.error("Need a protein expression filename and a shortname")

dataset = results.dataset_file
dataset_name = results.dataset_name
ppi_dataset = results.protein_interactions_file
PERMUTATIONS = results.randomisations

print "Using expression dataset", dataset, "with short name", dataset_name
print "Using Protein interaction dataset", ppi_dataset

expression_data = pd.read_csv(dataset, index_col=0, delimiter="\t")
gene_names = defaultdict(str)
with open("entrez_to_symbol.txt","rU") as f:
    for i in f :
        parts = i.strip().split()
        gene_names[int(parts[0])] = parts[1]
expression_data.columns = [int(x) for x in expression_data.columns]
g = nx.Graph()
with open(ppi_dataset,"rU") as f :
    for i in f :
        parts = i.strip().split("\t")
        gene1 = int(parts[0])
        gene2 = int(parts[1])
        llr = float(parts[2])
        g.add_edge(gene1,gene2,weight = llr)
        
print "PPI graph has", len(g.edges()), "edges and", len(g.nodes()), "nodes"

overlap = set([x for x in g.nodes() if x in expression_data])

g = nx.subgraph(g, overlap)

print "PPI graph has", len(g.edges()), "edges and", len(overlap),\
    "nodes when only considering genes in %s" % dataset_name
    
for e in g.selfloop_edges():
    g.remove_edge(e[0], e[1])
    
print "PPI graph has", len(g.edges()), "edges and", len(g.nodes()),\
    "nodes after removing self loops"
    

with open("PPI_graph_filtered.txt","w") as f :
    for edge in g.edges() :
        f.write("%s\t%s\n" % edge)
with open("Background_gene_set.txt","w") as f:
    for n in g.nodes() :
        f.write("%s\n" % n)

overlap = set([x for x in g.nodes() if x in expression_data])
missing_llr = math.log(1./300)

def get_ordered_tuple(item_a, item_b) :
    """
    Takes two items as input & returns an ordered tuple as the result.
    
    Useful for interactions so that (A,B) & (B,A) are not treated seperately
    
    """
    if cmp(item_a,item_b) <= 0:
        return (item_a, item_b)
    else :
        return (item_b, item_a)

def score_ppis(group1, group2):
    total = 0
    for i in group1:
        for j in group2:
            if g.has_edge(i, j):
                total += g[i][j]['weight']
            else:
                total += missing_llr
    return total


def score_correlations(group1, group2):
    corrs = []
    total = 0.0
    for i in group1:
        corrs.extend(expression_data[i][group2])
    return sum(corrs)


def get_integrated_score(group1, group2):
    return score_ppis(group1, group2) + score_correlations(group1, group2)


def score_cluster(cluster):
    group = list(cluster)
    total = 0
    for g1 in range(len(group)):
        for g2 in range(g1+1, len(group)):
            total+=get_integrated_score([group[g1]], [group[g2]])
    return total


def can_connect(group1, group2):
    subgraph = nx.subgraph(g, group1.union(group2))
    if len(subgraph) == 0 :
        print "Null graph", group1, group2
        return False
    if nx.is_connected(subgraph):
        if score_ppis(group1, group2) > 0:
            return score_correlations(group1, group2) > 0
        else:
            return False
    else:
        return False


def simple_can_connect(group1,group2):
    return nx.is_connected(nx.subgraph(g, group1.union(group2)))

def get_triplets(graph) :
    ''' Returns all 3-cliques in the graph
    '''
    triplets = set()
    done = set()
    for n1 in graph :
        neighbours = set(graph[n1])
        for n2 in neighbours :
            for n3 in neighbours.intersection(graph[n2]) :
                triplets.add(tuple(sorted([n1,n2,n3])))
    return triplets
    

for randomisation in range(PERMUTATIONS+1):
    if randomisation !=0:
        # Randomise correlation matrix ids
        print "Performing randomisation %s" % randomisation
        id_index = list(expression_data.index)
        random.shuffle(id_index)
        expression_data.index = id_index
        expression_data.columns = id_index
        sys.stdout.flush()
    
    temp_graph = nx.Graph()
    temp_score = {}
    for i in g.edges():
        if score_correlations([i[0]],[i[1]]) > 0: 
            temp_graph.add_edge(i[0],i[1])
    
    for i in temp_graph.edges() :
        if score_correlations([i[0]],[i[1]]) > 0:
            temp_score[get_ordered_tuple(i[0], i[1])] = get_integrated_score([i[0]], [i[1]])
    
    for i in get_triplets(temp_graph) :
        temp_score[i] = score_cluster(i)
    
    sorted_scores = sorted(temp_score.iteritems(),key=lambda(x,y):(y,x),reverse=True)
    index = max(overlap) + 1
    clusters = defaultdict(set)
    assigned = set()
    for i in sorted_scores :
        if not assigned.intersection(i[0]) :
            clusters[index] = set(i[0])
            assigned = assigned.union(set(i[0]))
            index += 1
    for i in overlap.difference(assigned):
        clusters[i].add(i)
    cluster_distance = {}
    # initially we only consider connected nodes
    for c1 in clusters.keys():
        for c2 in clusters.keys():
            if c1 < c2 and can_connect(clusters[c1], clusters[c2]):
                cluster_distance[get_ordered_tuple(c1, c2)] = get_integrated_score(clusters[c1], clusters[c2])
    threshold = 0
    max_pair, max_sim = max(cluster_distance.iteritems(), key=operator.itemgetter(1))
    move = "merge"
    removal_scores = {}
    switch_scores = {}
    all_distances = cluster_distance.copy()
    
    while max_sim > threshold and len(clusters) > 1:
        to_update = []
        if move == "merge":
            clusters[index] = clusters[max_pair[0]].union(clusters[max_pair[1]])
            del clusters[max_pair[0]]
            del clusters[max_pair[1]]
            to_update.append(index)
        elif move == "remove":
            clusters[max_pair[1]].remove(max_pair[0])
            clusters[max_pair[0]] = set([max_pair[0]])
            to_update.append(max_pair[0])
            to_update.append(max_pair[1])
        elif move == "switch":
            for c in clusters.keys():
                if max_pair[0] in clusters[c]:
                    clusters[c].remove(max_pair[0])
                    to_update.append(c)
            clusters[max_pair[1]].add(max_pair[0])
            to_update.append(max_pair[1])
        
        touched = set(to_update).union(max_pair)
        for c in cluster_distance.keys():
            if c[0] in touched or c[1] in touched:
                del cluster_distance[c]
        for c in switch_scores.keys():
            if c[0] in touched or c[1] in touched:
                del switch_scores[c]
        for c in removal_scores.keys():
            if c[0] in touched or c[1] in touched:
                del removal_scores[c]
        
        for update in to_update:
            for c in clusters.keys():
                if c!=update and can_connect(clusters[c], clusters[update]):
                    cluster_distance[get_ordered_tuple(c, update)] = get_integrated_score(clusters[c], clusters[update])
        
        multinode_clusters = [x for x in clusters if len(clusters[x])>1]
        # Get score for removing element from cluster
        for c in multinode_clusters:
            for gene in clusters[c]:
                if (gene,c) not in removal_scores:
                    if simple_can_connect(clusters[c].difference([gene]), set()) or len(clusters[c].difference(set([gene]))) < 2:
                        score = -1 * get_integrated_score([gene], clusters[c].difference([gene]))
                        removal_scores[(gene,c)] = score
                    else:
                        removal_scores[(gene,c)] = -1000000
        
        # Get score for switching element from cluster
        for c in multinode_clusters:
            for gene in clusters[c]:
                for other_clust in clusters:
                    if c!=other_clust:
                        if (gene,other_clust) not in switch_scores:
                            if can_connect(set([gene]), clusters[other_clust]):
                                score = get_integrated_score(set([gene]), clusters[other_clust])
                                switch_scores[(gene,other_clust)] = score + removal_scores[(gene,c)]
                            else:
                                switch_scores[(gene,other_clust)] = -1000000
        index+=1
        max_sim = 0
        move = "merge"
        
        if len(cluster_distance) > 0:
            max_pair, max_sim = max(cluster_distance.iteritems(), key=operator.itemgetter(1))
        if len(switch_scores) > 0:
            max_pair_switch, max_sim_switch = max(switch_scores.iteritems(), key=operator.itemgetter(1))
            if max_sim_switch > max_sim:
                max_sim = max_sim_switch
                max_pair = max_pair_switch
                move = "switch"
                #print "SWITCH", max_sim, max_pair
        if len(removal_scores) > 0:
            max_pair_remove, max_sim_remove = max(removal_scores.iteritems(), key=operator.itemgetter(1))
            if max_sim_remove > max_sim:
                max_sim = max_sim_remove
                max_pair = max_pair_remove
                move = "remove"
                #print "REMOVE", max_sim, max_pair
    
    cluster_scores = {}
    for x in [c for c in clusters if len(clusters[c])>1]:
        cluster_scores[x] = score_cluster(clusters[x])
    if randomisation == 0:
        real_clusters = deepcopy(clusters)
        real_cluster_scores = deepcopy(cluster_scores)
    else:
        with open("results/Randomised_clusters_%s_%s_%s.txt" % (dataset_name, ppi_dataset.split('.')[0], randomisation), "w") as f:
            for x in cluster_scores:
                f.write("%s\t%s\n" % (len(clusters[x]), cluster_scores[x]))

sorted_clusters = sorted(real_cluster_scores.items(), key=operator.itemgetter(1))
with open("results/%s_%s_clusters_common.txt" % (dataset_name, ppi_dataset.split('.')[0]), "w") as f:
    for s in sorted_clusters:
        f.write("%s\t%s\t%s\n" % (s[0], s[1], ",".join([gene_names[gene] for gene in real_clusters[s[0]]])))
with open("results/%s_%s_%s_clusters_entrez.txt" % (dataset_name, ppi_dataset.split('.')[0]), "w") as f:
    for s in sorted_clusters:
        f.write("%s\t%s\t%s\n" % (s[0], s[1], ",".join([str(gene) for gene in real_clusters[s[0]]])))
