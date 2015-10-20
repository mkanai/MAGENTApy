import numpy as np

from .utils import *


def ExtractGeneScoreBestSNP_PvalZscore_NumSNPsPerGene(GeneSubsetChrPos, All_SNP_scores_pos, interval_up, interval_down, strand, best_pval_or_z, SNP_rs):
    num_genes, c = GeneSubsetChrPos.shape
    Scores_X = np.zeros((num_genes, 4))
    Best_SNP_rs = np.zeros(num_genes, dtype="a{0}".format(MAXIDLEN))
    num_SNPs_per_gene = np.zeros(num_genes, dtype="i")
    count_num_Genes_NaN_scores = 0

    for gene in range(num_genes):
        if strand[gene] == 1:
            find_SNPs_near_gene = logical_and(np.equal(All_SNP_scores_pos[:, 0], GeneSubsetChrPos[gene, 0]),
                                              np.greater_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 1]-interval_up)),
                                              np.less_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 2]+interval_down)))
        else:
            find_SNPs_near_gene = logical_and(np.equal(All_SNP_scores_pos[:, 0], GeneSubsetChrPos[gene, 0]),
                                              np.greater_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 1]-interval_down)),
                                              np.less_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 2]+interval_up)))

        if np.any(find_SNPs_near_gene):
            num_stat_not_NaN = np.sum(np.logical_not(np.isnan(All_SNP_scores_pos[find_SNPs_near_gene][:, 2])))
            if num_stat_not_NaN == 0:
                count_num_Genes_NaN_scores = count_num_Genes_NaN_scores + 1
                Scores_X[gene, :] = np.tile(np.nan, 4)
                Best_SNP_rs[gene] = ''
                num_SNPs_per_gene[gene] = np.nan
            else:
                num_SNPs_per_gene[gene] = num_stat_not_NaN

                if best_pval_or_z == 0:
                    max_pos_val = np.max(All_SNP_scores_pos[find_SNPs_near_gene][:, 2])
                    min_neg_val = np.min(All_SNP_scores_pos[find_SNPs_near_gene][:, 2])
                    if abs(max_pos_val) >= abs(min_neg_val):
                        best_val = max_pos_val
                    else:
                        best_val = min_neg_val

                    best_val = max_pos_val

                    if strand[gene] == 1:
                        find_max_abs_Zscore = logical_and(np.equal(All_SNP_scores_pos[:, 2], best_val),
                                                          np.equal(All_SNP_scores_pos[:, 0], GeneSubsetChrPos[gene, 0]),
                                                          np.greater_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 1]-interval_up)),
                                                          np.less_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 2]+interval_down)))
                    else:
                        find_max_abs_Zscore = logical_and(np.equal(All_SNP_scores_pos[:, 2], best_val),
                                                          np.equal(All_SNP_scores_pos[:, 0], GeneSubsetChrPos[gene, 0]),
                                                          np.greater_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 1]-interval_down)),
                                                          np.less_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 2]+interval_up)))
                else:
                    min_pval = np.min(All_SNP_scores_pos[find_SNPs_near_gene][:, 3])

                    if strand[gene] == 1:
                        find_max_abs_Zscore = logical_and(np.equal(All_SNP_scores_pos[:, 3], min_pval),
                                                          np.equal(All_SNP_scores_pos[:, 0], GeneSubsetChrPos[gene, 0]),
                                                          np.greater_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 1]-interval_up)),
                                                          np.less_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 2]+interval_down)))
                    else:
                        find_max_abs_Zscore = logical_and(np.equal(All_SNP_scores_pos[:, 3], min_pval),
                                                          np.equal(All_SNP_scores_pos[:, 0], GeneSubsetChrPos[gene, 0]),
                                                          np.greater_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 1]-interval_down)),
                                                          np.less_equal(All_SNP_scores_pos[:, 1], (GeneSubsetChrPos[gene, 2]+interval_up)))

                Scores_X[gene, :] = All_SNP_scores_pos[find_max_abs_Zscore][0, :4]

                if len(SNP_rs) > 1:
                    Best_SNP_rs[gene] = (SNP_rs[find_max_abs_Zscore])[0]
                else:
                    Best_SNP_rs[gene] = ''
        else:
            Scores_X[gene, :] = np.tile(np.nan, 4)
            Best_SNP_rs[gene] = ''

    print('There are {0:d} genes with SNPs in their target region that were not assigned an association score (NaN).\n'.format(count_num_Genes_NaN_scores))
    return (Scores_X, num_SNPs_per_gene, Best_SNP_rs)
