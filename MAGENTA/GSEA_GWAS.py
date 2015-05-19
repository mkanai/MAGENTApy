import numpy as np
from utils import *


def GSEA_GWAS(Uncorr_score, Corr_score, top_percen_cutoffs, num_rounds, find_gene_set, score_signif_direct, choose_unique_genes, num_perm_limit):
    initial_num_genes_subset = len(find_gene_set)

    original_num_rounds = num_rounds
    b2 = find_gene_set

    GeneSetScore = np.abs(Corr_score[b2])
    Uncorr_score_geneset = np.abs(Uncorr_score[b2, :])
    GeneSetScore = GeneSetScore[np.logical_not(np.isnan(GeneSetScore))]

    Uncorr_score_all_sorted_ind = np.argsort(np.abs(Corr_score))
    # Sorted_all_uncorr_SNP_Chrpos = Uncorr_score[Uncorr_score_all_sorted_ind, :2]
    Sorted_all_uncorr_SNP_Chrpos = Uncorr_score[Uncorr_score_all_sorted_ind, 0] * 1e11 + Uncorr_score[Uncorr_score_all_sorted_ind, 1]
    score = Corr_score[Uncorr_score_all_sorted_ind]

    cutoffs = np.percentile(score, top_percen_cutoffs)
    cutoffs_size = len(cutoffs)
    gene_set_size = len(GeneSetScore)

    obs_num_above_cutoff = np.zeros(cutoffs_size)
    obs_fract_above_cutoff = np.zeros(cutoffs_size)
    count_num_permut_above_perc_cutoff = np.zeros(cutoffs_size)

    if score_signif_direct == 1:
        obs_num_above_cutoff = np.sum(np.less_equal(GeneSetScore, cutoffs[:, None]), axis=1)
    elif score_signif_direct == 0:
        obs_num_above_cutoff = np.sum(np.greater_equal(GeneSetScore, cutoffs[:, None]), axis=1)
    obs_fract_above_cutoff = 1.0 * obs_num_above_cutoff / gene_set_size

    fract_permut_above_perc_cutoff = np.ones(cutoffs_size) * (0.99 / num_rounds)
    record_rand_gs_num = np.zeros((original_num_rounds, cutoffs_size))
    record_rand_gs_fract = np.zeros((original_num_rounds, cutoffs_size))

    score_indexes = np.arange(len(score))
    Updated_gene_ind = np.zeros(gene_set_size, dtype="int")

    while np.any(np.less_equal(fract_permut_above_perc_cutoff, 1.0 / num_rounds)):
        for i in range(num_rounds):
            # rand_geneset_ind = np.random.choice(score_indexes, gene_set_size, replace=False)
            rand_geneset_ind = np.unique(np.random.randint(len(score), size=gene_set_size))
            if choose_unique_genes == 1:
                # temp = Sorted_all_uncorr_SNP_Chrpos[rand_geneset_ind, 0] * 1e11 + Sorted_all_uncorr_SNP_Chrpos[rand_geneset_ind, 1]
                temp = Sorted_all_uncorr_SNP_Chrpos[rand_geneset_ind]
                _, Unique_genes_ind = np.unique(temp, return_index=True)

                Saved_gene_ind = rand_geneset_ind[Unique_genes_ind]

                saved_gene_size = len(Saved_gene_ind)
                while saved_gene_size < gene_set_size:
                    remain_num = gene_set_size - saved_gene_size
                    # setdiff
                    # remain_ind = np.where(np.logical_not(np.in1d(np.arange(len(score)), Saved_gene_ind)))[0]
                    remain_ind = np.setdiff1d(score_indexes, Saved_gene_ind, assume_unique=True)

                    # rand_geneset_ind_add = np.random.choice(remain_ind, remain_num, replace=False)
                    rand_geneset_ind_add = remain_ind[np.unique(np.random.randint(len(remain_ind), size=remain_num))]
                    # Updated_gene_ind = np.unique(np.concatenate((Saved_gene_ind, rand_geneset_ind_add)))
                    # Updated_gene_ind = np.union1d(Saved_gene_ind, rand_geneset_ind_add)
                    # Updated_gene_ind = np.concatenate((Saved_gene_ind, rand_geneset_ind_add))
                    Updated_gene_ind[:saved_gene_size] = Saved_gene_ind
                    Updated_gene_ind[saved_gene_size:] = rand_geneset_ind_add
                    # temp = Sorted_all_uncorr_SNP_Chrpos[Updated_gene_ind, 0] * 1e11 + Sorted_all_uncorr_SNP_Chrpos[Updated_gene_ind, 1]
                    temp = Sorted_all_uncorr_SNP_Chrpos[Updated_gene_ind]
                    _, Unique_genes_ind = np.unique(temp, return_index=True)
                    Saved_gene_ind = Updated_gene_ind[Unique_genes_ind]
                    saved_gene_size = len(Saved_gene_ind)
            else:
                Saved_gene_ind = rand_geneset_ind

            rand_geneset_score = score[Saved_gene_ind]
            rand_geneset_score_size = len(rand_geneset_score)

            if score_signif_direct == 1:
                num_above_cutoffs = np.sum(np.less_equal(rand_geneset_score, cutoffs[:, None]), axis=1)
                count_num_permut_above_perc_cutoff += np.greater_equal(num_above_cutoffs, obs_num_above_cutoff)
            elif score_signif_direct == 0:
                num_above_cutoffs = np.sum(np.greater_equal(rand_geneset_score, cutoffs[:, None]), axis=1)
                count_num_permut_above_perc_cutoff += np.greater_equal(num_above_cutoffs, obs_num_above_cutoff)

            if num_rounds <= 10000:
                record_rand_gs_num[i, :] = num_above_cutoffs
                record_rand_gs_fract[i, :] = 1.0 * num_above_cutoffs / rand_geneset_score_size

        fract_permut_above_perc_cutoff = 1.0 * count_num_permut_above_perc_cutoff / num_rounds

        num_rounds *= 100
        if num_rounds > num_perm_limit:
            fract_permut_above_perc_cutoff[np.less_equal(fract_permut_above_perc_cutoff, 1.0 / num_rounds)] = 999

    fract_permut_above_perc_cutoff[np.equal(fract_permut_above_perc_cutoff, 999)] = 0.99 / (num_rounds / 100)
    EnrichPval_multipleCutoffs = fract_permut_above_perc_cutoff
    record_rand_gs_num_output = record_rand_gs_num[:original_num_rounds, :]
    record_rand_gs_fract_output = record_rand_gs_fract[:original_num_rounds, :]

    return (EnrichPval_multipleCutoffs, obs_num_above_cutoff, obs_fract_above_cutoff, record_rand_gs_num_output, record_rand_gs_fract_output)


def GSEA_RankSumStat(Corr_score, find_rnd_gene_set, score_signif_direct):

    find_outside_gene_set = np.where(np.logical_not(np.in1d(np.arange(len(Corr_score)), find_rnd_gene_set)))[0]

    RankSum_stat, RankSum_pval = sp.stats.ranksums(np.abs(Corr_score[find_rnd_gene_set]), np.abs(Corr_score[find_outside_gene_set]))
    GSEA_RankSum_z_pval = np.array((-RankSum_stat, RankSum_pval / 2.0 if score_signif_direct == 1 else -RankSum_pval / 2.0))

    return GSEA_RankSum_z_pval


def GSEA_GWAS_RankSum(Uncorr_score, Corr_score, top_percen_cutoffs, num_rounds, find_gene_set, score_signif_direct, choose_unique_genes):

    initial_num_genes_subset = len(find_gene_set)

    b2 = find_gene_set

    GeneSetScore = np.abs(Corr_score[b2])
    Uncorr_score_geneset = np.abs(Uncorr_score[b2, :])

    Uncorr_score_geneset = Uncorr_score_geneset[np.logical_not(np.isnan(GeneSetScore))]
    GeneSetScore = GeneSetScore[np.logical_not(np.isnan(GeneSetScore))]

    Uncorr_score_all_sorted_ind = np.argsort(np.abs(Corr_score))
    Sorted_all_uncorr_SNP_Chrpos = Uncorr_score[Uncorr_score_all_sorted_ind, :2]
    score = Corr_score[Uncorr_score_all_sorted_ind]

    cutoffs = np.percentile(score, top_percen_cutoffs)
    gene_set_size = len(GeneSetScore)

    rand_geneset_score = np.tile(np.nan, (gene_set_size, num_rounds))
    rand_geneset_find = np.tile(np.nan, (gene_set_size, num_rounds))

    for i in range(num_rounds):
        rand_geneset_ind = np.random.choice(np.arange(len(score)), gene_set_size)
        Unique_genes = np.ascontiguousarray(Sorted_all_uncorr_SNP_Chrpos[rand_geneset_ind, :2]).view(np.dtype((np.void, Sorted_all_uncorr_SNP_Chrpos[rand_geneset_ind, :2].dtype.itemsize * Sorted_all_uncorr_SNP_Chrpos[rand_geneset_ind, :2].shape[1])))
        _, Unique_genes_ind = np.unique(Unique_genes, return_index=True)

        while len(Saved_gene_ind) < gene_set_size:
            remain_num = gene_set_size- len(Saved_gene_ind)
            remain_ind = np.where(np.logical_not(np.in1d(np.arange(len(score)), Saved_gene_ind)))[0]

            rand_geneset_ind_add = np.random.choice(remain_ind, remain_num, replace=False)
            # Updated_gene_ind = np.unique(np.concatenate((Saved_gene_ind, rand_geneset_ind_add)))
            Updated_gene_ind = np.union1d(Saved_gene_ind, rand_geneset_ind_add)
            Unique_genes = np.ascontiguousarray(Sorted_all_uncorr_SNP_Chrpos[Updated_gene_ind, :2]).view(np.dtype((np.void, Sorted_all_uncorr_SNP_Chrpos[Updated_gene_ind, :2].dtype.itemsize * Sorted_all_uncorr_SNP_Chrpos[Updated_gene_ind, :2].shape[1])))
            _, Unique_genes_ind = np.unique(Unique_genes, return_index=True)
            Saved_gene_ind = Updated_gene_ind[Unique_genes_ind]

        rand_geneset_score[:, i] = score[Saved_gene_ind]
        rand_geneset_find[:, i] = Saved_gene_ind.T

    Obs_GS_RankSum_z_pval = GSEA_RankSumStat(Corr_score, find_gene_set, score_signif_direct)

    for p in range(num_rounds):
        Rand_GS_RankSum_z_pval[p, :2] = GSEA_RankSumStat(Corr_score,rand_geneset_find[:, p], score_signif_direct)

    GSEA_p[0] = 1.0 * np.sum(np.greater_equal(Rand_GS_RankSum_z_pval[:, 0], Obs_GS_RankSum_z_pval[0])) / num_rounds
    GSEA_p[1] = 1.0 * np.sum(np.less(Rand_GS_RankSum_z_pval[:, 0], Obs_GS_RankSum_z_pval[0])) / num_rounds

    if GSEA_p[0] == 0:
        GSEA_p[0] = 0.99 / num_rounds

    if GSEA_p[1] == 0:
        GSEA_p[1] = 0.99 / num_rounds

    if Obs_GS_RankSum_z_pval[0] >= 0:
        GSEA_p[2] = 1
    else:
        GSEA_p[2] = 0

    return GSEA_p
