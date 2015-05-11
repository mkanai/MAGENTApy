import numpy as np
from utils import *


def multi_list2cell_GeneSetDB(filename):
    output_cellarray = np.zeros(BUFLINENUM, dtype=('a{},a{},{}i'.format(MAXIDLEN, MAXIDLEN, MAXGENESETNUM)))
    num_genes_per_gene_set = np.zeros(BUFLINENUM)
    geneset = np.tile(-1, MAXGENESETNUM)
    counter = 0
    with open(filename, "r") as reader:
        for line in reader:
            tokens = line.strip('\n').split('\t')
            ntokens = len(tokens)
            num_genes_per_gene_set[counter] = ntokens - 2
            if num_genes_per_gene_set[counter] > MAXGENESETNUM:
                print num_genes_per_gene_set[counter]
            geneset[:num_genes_per_gene_set[counter]] = map(lambda x: int(x), tokens[2:])
            output_cellarray[counter] = (tokens[0], tokens[1], geneset)
            counter += 1
            geneset.fill(-1)
            if (counter % BUFLINENUM) + 1 == BUFLINENUM:
                output_cellarray = np.concatenate((output_cellarray, np.zeros(BUFLINENUM, dtype=output_cellarray.dtype)))
                num_genes_per_gene_set = np.concatenate((num_genes_per_gene_set, np.zeros(BUFLINENUM)))
    return (output_cellarray, num_genes_per_gene_set)
