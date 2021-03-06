import numpy as np
from utils import *


def multi_list2cell_GeneSetDB(filename):
    output_cellarray = np.zeros(BUFLINENUM, dtype=('a{0},a{1},{2}i'.format(MAXIDLEN, MAXIDLEN, MAXGENENUM_PER_GENESET)))
    num_genes_per_gene_set = np.zeros(BUFLINENUM)
    geneset = np.tile(-1, MAXGENENUM_PER_GENESET)
    counter = 0
    with open(filename, "r") as reader:
        for line in reader:
            tokens = np.array(line.strip('\n').split('\t'))
            genes = np.unique(tokens[2:])
            num_genes_per_gene_set[counter] = len(genes) 
            if num_genes_per_gene_set[counter] > MAXGENENUM_PER_GENESET:
                raise ValueError('A number of genes in ({0}, {1}): {2} exceeds MAXGENENUM_PER_GENESET. Try increase it.'.format(tokens[0], tokens[1], num_genes_per_gene_set[counter]))
            geneset[:num_genes_per_gene_set[counter]] = map(lambda x: int(x), genes)
            output_cellarray[counter] = (tokens[0], tokens[1], geneset)
            counter += 1
            geneset.fill(-1)
            if (counter % BUFLINENUM) + 1 == BUFLINENUM:
                output_cellarray = np.concatenate((output_cellarray, np.zeros(BUFLINENUM, dtype=output_cellarray.dtype)))
                num_genes_per_gene_set = np.concatenate((num_genes_per_gene_set, np.zeros(BUFLINENUM)))
    return (output_cellarray, num_genes_per_gene_set)
