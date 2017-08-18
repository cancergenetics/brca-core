import pandas as pd
import numpy as np
import argparse
import os
import csv

parser = argparse.ArgumentParser(description='This compares the output of real vs randomised expression/PPI integration')
parser.add_argument('-n', action='store', dest='dataset_name', 
    help='A shortname for the dataset')
parser.add_argument('-p', action='store', dest='ppi_dataset', type=str,
    help='The PPI dataset', default="weighted_ppis.txt")
parser.add_argument('-f', action='store', dest='fdr', type=float, 
    help='The False Discovery Rate to be used', default=0.1)
results = parser.parse_args()
if not results.dataset_name:
    parser.error("Need a shortname")

# Edit these for each different data resource
dataset_name = results.dataset_name
ppi_dataset = results.ppi_dataset.split('.')[0]
fdr_threshold = results.fdr

rand_output_name = "Randomised_clusters_%s_%s" % (dataset_name, ppi_dataset)

rand_runs = []

for i in os.listdir("./") :
    if i.startswith(rand_output_name) :
        with open(i, "r") as f :
            llrs = []
            for x in f :
                parts = x.split()
                llrs.append((int(parts[0]),float(parts[1])))
            rand_runs.append(llrs)
            print i, len(llrs)
print len(rand_runs)

real_output_name = "%s_%s_clusters_entrez.txt" % (dataset_name, ppi_dataset)
 
print "Comparing ", len(rand_runs), "random runs against", real_output_name

cluster_scores = {}
cluster_members = {}
with open(real_output_name,"r") as f: 
    reader = csv.reader(f,delimiter="\t")
    for line in reader :
        cluster_scores[line[0]] = float(line[1])
        cluster_members[line[0]] = set(line[2].split(","))

observed_scores = sorted(cluster_scores.values(),reverse=True) 

threshold = 0
for i in observed_scores :
    real_count = sum([len(cluster_members[x]) for x in cluster_scores if cluster_scores[x] >= i])
    random_counts = []
    for run in rand_runs :
        random_counts.append(sum([x[0] for x in run if x[1] >= i]))
    fdr = np.median(random_counts) / float(real_count)
    if fdr < fdr_threshold :
        threshold = i
        print "At a threshold of %s we observe %s real vs %s expected genes in clusters" % (i, real_count, np.median(random_counts))          
        print "We estimate the FDR at ", fdr

output_name = "%s_%s_clusters_entrez_fdr%s.txt" % (dataset_name, ppi_dataset,fdr_threshold)
with open(output_name,"w") as f :
    for i in cluster_scores :
        if cluster_scores[i] >= threshold :
            f.write("%s\t%s\t%s\n" % (i,cluster_scores[i],",".join(cluster_members[i])))