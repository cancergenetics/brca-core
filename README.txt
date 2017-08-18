The input files required for this script can (temporarily) be obtained from https://www.dropbox.com/s/pbuc87swsx6eb22/input_files.zip?dl=0

The clustering.py script integrates protein expression with PPI data to identify potential protein complexes. It requires a weighted protein protein interaction network and a weighted protein co-expression matrix. Both are expected to contain log-likelihood ratios. 

Run "python clustering.py -h" to see the below output

usage: clustering.py [-h] [-f DATASET_FILE] [-n DATASET_NAME]
                     [-p PROTEIN_INTERACTIONS_FILE] [-d DEGREE]
                     [-r RANDOMISATIONS]

optional arguments:
  -h, --help            show this help message and exit
  -f DATASET_FILE       The protein expression dataset filename
  -n DATASET_NAME       A shortname for the dataset
  -p PROTEIN_INTERACTIONS_FILE
                        The weighted protein interaction dataset name
  -d DEGREE             Threshold for discarding hubs from the protein
                        interaction network (not used currently)
  -r RANDOMISATIONS     Number of randomisations to perform to estimate false
                        discovery rate (default is 10)

Using the sample data you can run clustering.py as follows :

python clustering.py -f brca_llrs.txt -n TCGA_BRCA -p weighted_ppis.txt

This will create a putative clustering in the results folder (TCGA_BRCA_weighted_ppis_entrez.txt) along with clusterings generated using random permutations (e.g. Randomised_clusters_TCGA_BRCA_weighted_ppis_1.txt). This will also output a filtered protein interaction network (containing only proteins present in the expression set) and a file containing a list of proteins common to the expression set and interaction network (which can be used as background for gene ontology enrichment).

The process_results.py script takes the output of the clustering.py script and estimates the false discovery rate (FDR) for the identified complexes. By default it will output those complexes with an estimated FDR < 10%

Run "python process_results.py -h" to see the below output

usage: process_results.py [-h] [-n DATASET_NAME] [-p PPI_DATASET] [-f FDR]

This compares the output of real vs randomised expression/PPI integration

optional arguments:
  -h, --help       show this help message and exit
  -n DATASET_NAME  A shortname for the dataset
  -p PPI_DATASET   The PPI dataset
  -f FDR           The False Discovery Rate to be used (defaults to 10%)
  
Using the sample data you can run clustering.py as follows :

python process_results.py -n TCGA_BRCA -p weighted_ppis.txt

This will create two new files in the results folder (TCGA_BRCA_weighted_ppis_entrez_fdr0.1.txt) and (TCGA_BRCA_weighted_ppis_common_fdr0.1.txt) which contain the complexes identified at the chosen FDR. The _entrez_ file contains the complexes annotated using their ENTREZ IDs while _common_ uses common gene names from HGNC 
