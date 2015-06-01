# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np
import pandas as pd
import click
import datetime
import time
import os
import yaml

from MAGENTA.utils import *
from MAGENTA.Converting_GWAS_TwoTailedPval_to_Zscores import *
from MAGENTA.ExtractGeneScoreBestSNP_PvalZscore_NumSNPsPerGene import *
from MAGENTA.Calculate_NumIndepSNPsPerGene import *
from MAGENTA.Calculate_NumHotspots_per_gene import *
from MAGENTA.GeneScoreRegressCorr import *
from MAGENTA.multi_list2cell_GeneSetDB import *
from MAGENTA.GSEA_GWAS import *


@click.command()
@click.option('--config', default="config.yml", help='Configuration file.')
def magenta(config):
    with open(config, 'r') as f:
        config = yaml.load(f)

    GWAS_SNP_score_file_name = config['GWAS_SNP_score_file_name']
    GeneSet_db_file_name = config['GeneSet_db_file_name']
    Flag_gene_set_file_name = config['Flag_gene_set_file_name']
    match_genes_GeneID = config['match_genes_GeneID']
    GWAS_OR_filename = config['GWAS_OR_filename']
    SNP_rs_num_file_name = config['SNP_rs_num_file_name']
    exp_label = config['exp_label']
    code_name = config['code_name']
    Genome_build = config['Genome_build']
    Flag_gs = config['Flag_gs']
    Remove_gs = config['Remove_gs']
    Remove_HLA = config['Remove_HLA']
    HLA_region = config['HLA_region']
    num_gs_simul = config['num_gs_simul']
    num_perm_limit = config['num_perm_limit']
    Gene_boundr_upstr = config['Gene_boundr_upstr']
    Gene_boundr_downstr = config['Gene_boundr_downstr']
    run_GSEA = config['run_GSEA']
    GSEA_method = config['GSEA_method']
    top_percen_cutoffs = config['top_percen_cutoffs']
    min_gs_size = config['min_gs_size']
    max_gs_size = config['max_gs_size']
    choose_unique_genes = config['choose_unique_genes']
    print_trait_gene_names_in_gs = config['print_trait_gene_names_in_gs']
    print_gene_scores = config['print_gene_scores']
    print_rs_num = config['print_rs_num']
    print_best_SNP_OR = config['print_best_SNP_OR']
    calculate_FDR = config['calculate_FDR']
    seed = config['seed']
    if seed == "time" or seed is None:
        seed = time.time()
    seed = int(seed)
    np.random.seed(seed)

    analysis_start_time = time.time()
    if Genome_build == 'NCBI36':
        HumanGeneChrPos = pd.read_csv('AllHumanGeneChrPosStrand_18434Genes_RefFlat_111909', header=None, delim_whitespace=True)
        hotspot_boundaries = pd.read_csv('hotspot_boundaries_b36', header=None, delim_whitespace=True)
        PrunedSNPsChrPos = pd.read_csv('CEUPruned_SNP_ChrNumPos_hg18_030811', header=None, delim_whitespace=True).values
        AllRefFlatGene = pd.read_csv('RefFlatGeneSymbolGeneID_18434Genes_111909', header=None, delim_whitespace=True).values
        AllRefFlatGeneID = AllRefFlatGene[:, 0]
        AllRefFlatGeneNames = AllRefFlatGene[:, 1]
    else:
        HumanGeneChrPos = pd.read_csv("AllHumanGeneChrPosStrand_RefSeq_hg19_072111", header=None, delim_whitespace=True).values
        hotspot_boundaries = pd.read_csv("hotspot_boundaries_b37_hg19_072111", header=None, delim_whitespace=True).values
        PrunedSNPsChrPos = pd.read_csv("CEU_HapMap_pruned_SNPs_ChrNumPos_hg19_072111", header=None, delim_whitespace=True).values
        AllRefFlatGene = pd.read_csv('AllHumanGeneNames_RefSeq_hg19_072111', header=None, delim_whitespace=True).values
        AllRefFlatGeneID = AllRefFlatGene[:, 0]
        AllRefFlatGeneNames = AllRefFlatGene[:, 1]

    calculate_Scores_confounders = 1
    GeneScoreMetric = 1
    score_signif_direct = 1

    if Remove_HLA == 1:
        find_not_HLA_reg = np.logical_not(inHLAregion(HumanGeneChrPos, HLA_region['st'], HLA_region['en']))
        AllGeneChrPos = HumanGeneChrPos[find_not_HLA_reg]
    else:
        AllGeneChrPos = HumanGeneChrPos


    todays_date = datetime.datetime.today().strftime("%b%d_%y")

    # % (0) PRINT INTO LOG FILE:

    # % make directory
    Output_dir = get_valid_filename('Output_MAGENTA_{0}_{1}perm_{2}'.format(exp_label, num_gs_simul, todays_date))
    if not os.path.exists(Output_dir):
        os.mkdir(Output_dir)

    print('All output data is saved in the following directory: {0}\n'.format(Output_dir))
    GSEA_log_file_name = '{0}/MAGENTA_run_{1}_{2}perm_{3}.log'.format(Output_dir, exp_label, num_gs_simul, todays_date)
    logger = Logger(GSEA_log_file_name)
    logger.log('MAGENTA written by Ayellet Segre, Altshuler and Daly Labs, Date: {0}\n\n', todays_date)
    logger.log('Summary file for running MAGENTA (Gene Set Enrichment Analysis (GSEA) on Genome-wide association study (GWAS) variant results).\n')
    logger.log('Program run: {0}\n', code_name)
    logger.log('GWAS used: {0}\n', exp_label)
    logger.log('Random seed: {0:d}\n', seed)
    logger.log('Number of randomly sampled gene sets for GSEA-GWAS p-value estimation is: {0}\n', num_gs_simul)
    logger.log('Gene boundaries used for mapping SNPs onto genes are: {0}kb upstream to most extreme gene transcript start position, and {1}kb downstream to most extreme gene transcript end position, taking gene orientation into account.\n', Gene_boundr_upstr / 1000, Gene_boundr_downstr / 1000)

    # (1) Open GSEA-GWAS results file

    if run_GSEA == 1:
        GSEA_results_file_name = get_valid_filename('MAGENTA_pval_GeneSetEnrichAnalysis_{0}_{1}kb_upstr_{2}kb_downstr_{3}perm_{4}.results'.format(exp_label, Gene_boundr_upstr / 1000, Gene_boundr_downstr / 1000, num_gs_simul, todays_date), Output_dir)
        results = Writer(GSEA_results_file_name, bufsize=1)

        if GSEA_method == 1:
            if calculate_FDR == 1:
                results.write('DB\tGS\tORIG_GS_SIZE\tEFF_GS_SIZE\t#_GENES_ABS_LIST\t#_GENES_WO_SCORE\t#_GENES_REM_PROXIM\tMED_GENE_SIZE_KB\tMEAN_GENE_SIZE_KB\tNOMINAL_GSEA_PVAL_95PERC_CUTOFF\tFDR_95PERC_CUTOFF\tEXP_#_GENES_ABOVE_95PERC_CUTOFF\tOBS_#_GENES_ABOVE_95PERC_CUTOFF\tNOMINAL_GSEA_PVAL_90PERC_CUTOFF\tFDR_90PERC_CUTOFF\tEXP_#_GENES_ABOVE_90PERC_CUTOFF\tOBS_#_GENES_ABOVE_90PERC_CUTOFF\t#_GENES_FLAGGED\tFLAGGED_GENE_NAMES\n')
            else:
                results.write('DB\tGS\tORIG_GS_SIZE\tEFF_GS_SIZE\t#_GENES_ABS_LIST\t#_GENES_WO_SCORE\t#_GENES_REM_PROXIM\tMED_GENE_SIZE_KB\tMEAN_GENE_SIZE_KB\tNOMINAL_GSEA_PVAL_95PERC_CUTOFF\tEXP_#_GENES_ABOVE_95PERC_CUTOFF\tOBS_#_GENES_ABOVE_95PERC_CUTOFF\tNOMINAL_GSEA_PVAL_90PERC_CUTOFF\tEXP_#_GENES_ABOVE_90PERC_CUTOFF\tOBS_#_GENES_ABOVE_90PERC_CUTOFF\t#_GENES_FLAGGED\tFLAGGED_GENE_NAMES\n')
        elif GSEA_method == 4:
            results.write('DB\tGS\tORIG_GS_SIZE\tEFF_GS_SIZE\t#_GENES_ABS_LIST\t#_GENES_WO_SCORE\t#_GENES_REM_PROXIM\tMED_GENE_SIZE_KB\tMEAN_GENE_SIZE_KB\tRANKSUM_PVAL_HIGHER_TAIL\t#_GENES_FLAGGED\tFLAGGED_GENE_NAMES\n')

    if match_genes_GeneID == 1:
        AllGene_Names = AllRefFlatGeneID
    else:
        AllGene_Names = AllRefFlatGeneNames

    if Remove_HLA == 1:
        AllGene_Names = AllGene_Names[find_not_HLA_reg]
        AllRefFlatGeneID = AllRefFlatGeneID[find_not_HLA_reg]
        AllRefFlatGeneNames = AllRefFlatGeneNames[find_not_HLA_reg]

    # % (3.1) Read in GWA SNP z-scores and p-values after genomic control
    GWAS_SNP_ChrNumPosZscoresPval = pd.read_csv(GWAS_SNP_score_file_name, header=None, delim_whitespace=True).values
    n_r, n_c = GWAS_SNP_ChrNumPosZscoresPval.shape

    if n_c < 4 and np.max(GWAS_SNP_ChrNumPosZscoresPval[:, 2]) <= 1:
        Input_file_SNPpval_only = 1
    else:
        Input_file_SNPpval_only = 0

    if print_best_SNP_OR == 1:
        if GWAS_OR_filename:
            GWAS_SNP_OR_L95_U95 = pd.read_csv(GWAS_OR_filename, header=None, delim_whitespace=True).values

    if Input_file_SNPpval_only == 1:
        GWAS_SNPChrNumPos_ZscoresCalculatedFromPval = Converting_GWAS_TwoTailedPval_to_Zscores(GWAS_SNP_ChrNumPosZscoresPval)
        GWAS_SNP_ChrNumPosZscoresPval = GWAS_SNPChrNumPos_ZscoresCalculatedFromPval


    #  (3.2) Read in all SNP rs numbers in same order as SNP positions and z-scores/p-values in input table one column of all SNP rs numbers

    if SNP_rs_num_file_name and print_rs_num == 1:
        AllSNP_rs_num = pd.read_csv(SNP_rs_num_file_name, header=None).values
    else:
        AllSNP_rs_num = np.zeros(len(GWAS_SNP_ChrNumPosZscoresPval), dtype="str")

    # (4) Calculate 'gene association score' from local SNP association z-scores Best SNP per gene scoring metric

    StrandOrient = AllGeneChrPos[:, 5]
    if GeneScoreMetric == 1:
        best_pval_or_z = 1
        BestSNPPerGeneChrNumPosZScores, NumSNPs_per_gene, Best_SNP_rs = ExtractGeneScoreBestSNP_PvalZscore_NumSNPsPerGene(AllGeneChrPos[:, :3], GWAS_SNP_ChrNumPosZscoresPval, Gene_boundr_upstr, Gene_boundr_downstr, StrandOrient, best_pval_or_z, AllSNP_rs_num)

        Uncorr_score = np.abs(BestSNPPerGeneChrNumPosZScores)

        # % (5) Correct gene score according to confounders using step-wise multivariate linear regression analysis
        interval = Gene_boundr_upstr+Gene_boundr_downstr
        AllGeneSizes = AllGeneChrPos[:, 2] - AllGeneChrPos[:, 1]
        GeneSize_plus_interval_in_kb = (AllGeneSizes + interval)/1000.0

        start_time = time.time()
        Num_Indep_SNPs_per_gene = Calculate_NumIndepSNPsPerGene(AllGeneChrPos, PrunedSNPsChrPos, Gene_boundr_upstr, Gene_boundr_downstr)
        end_time = time.time()

        logger.log('It took {0:.4f} seconds or {1:.4f} minutes to calculate number of independent SNPs per gene.\n', end_time-start_time, (end_time-start_time)/60)

        start_time = time.time()
        NumHotspots_per_gene = Calculate_NumHotspots_per_gene(AllGeneChrPos[:, :3], hotspot_boundaries, Gene_boundr_upstr, Gene_boundr_downstr, StrandOrient)
        end_time = time.time()
        logger.log('It took {0:.4f} seconds or {1:.4f} minutes to calculate number of hotspots per gene.\n', end_time-start_time, (end_time-start_time)/60)

        confounders = np.vstack((GeneSize_plus_interval_in_kb,
                                 NumSNPs_per_gene/GeneSize_plus_interval_in_kb,
                                 Num_Indep_SNPs_per_gene/GeneSize_plus_interval_in_kb,
                                 NumHotspots_per_gene/GeneSize_plus_interval_in_kb)).T

        start_time = time.time()
        RegressCorrGeneScores_pval, Residuals = GeneScoreRegressCorr(confounders, Uncorr_score[:, 2])
        Corr_score = RegressCorrGeneScores_pval
        end_time = time.time()
        logger.log('It took {0:.4f} seconds or {1:.4f} minutes to correct gene scores with stepwise regression.\n', end_time-start_time, (end_time-start_time)/60)
        np.savetxt(get_valid_filename('RegressCorrGeneScores_pval_{0}kb_upstr_{1}kb_downstr_{2}'.format(Gene_boundr_upstr/1000, Gene_boundr_downstr/1000, todays_date), Output_dir), RegressCorrGeneScores_pval)

    elif GeneScoreMetric == 2:  # TO BE ADDED
        raise NotImplementedError("TO BE IMPLEMENTED.")
        # SetScore = [];
        # Corr_score = SetScore(find_corr_score_isan);

    # % (6) Remove genes without gene scores (NaN), and if requested - remove input subset of genes from gene list and from analyses (e.g. known disease genes)
    # % Adjust all relevant gene matrices (e.g. Corrected and uncorrected gene scores, Gene chromosome positions, Gene IDs/Names accordingly

    find_corr_score_isan = np.logical_and(np.logical_not(np.isnan(Corr_score)), np.logical_not(np.isnan(Uncorr_score[:, 2])))
    AllGene_Names_score_isan = AllGene_Names[find_corr_score_isan]

    if Flag_gene_set_file_name:
        if Remove_gs == 1:
            if match_genes_GeneID == 1:
                Flag_gs_GeneNames = pd.read_csv(Flag_gene_set_file_name, header=None, delim_whitespace=True).values
            else:
                Flag_gs_GeneNames = pd.read_csv(Flag_gene_set_file_name, header=None, delim_whitespace=True).values  # well, there should be some difference, though

            Flag_gs_GeneNames = np.unique(Flag_gs_GeneNames)

            AllGenes_wo_flag_gs_ind = np.where(np.logical_not(np.in1d(AllGene_Names, Flag_gs_GeneNames)))[0]
            find_score_isan_not_flag_gs_unsorted = np.intersect1d(AllGenes_wo_flag_gs_ind, np.where(find_corr_score_isan))
            find_score_isan_not_flag_gs = np.sort(find_score_isan_not_flag_gs_unsorted)
    else:
        find_score_isan_not_flag_gs = find_corr_score_isan

    Corr_score_isan = Corr_score[find_score_isan_not_flag_gs]
    Uncorr_score_isan = np.abs(BestSNPPerGeneChrNumPosZScores[find_score_isan_not_flag_gs, :])

    AllGene_Names_isan = AllGene_Names[find_score_isan_not_flag_gs]
    NumSNPs_per_gene_isan = NumSNPs_per_gene[find_score_isan_not_flag_gs]
    AllGeneSizes_isan = AllGeneSizes[find_score_isan_not_flag_gs]



    # (7) Read in gene sets
    gs_GeneIDs_cellarray, num_genes_gs = multi_list2cell_GeneSetDB(GeneSet_db_file_name)
    find_gs_min_num_genes = np.logical_and(np.greater_equal(num_genes_gs, min_gs_size), np.less_equal(num_genes_gs, max_gs_size))
    num_gene_sets_min_num_genes = np.sum(find_gs_min_num_genes)
    find_gs_min_num_genes[len(gs_GeneIDs_cellarray)-1] = True
    gs_GeneIDs_cellarray_min_gene_num = gs_GeneIDs_cellarray[find_gs_min_num_genes]
    gs_GeneIDs_cellarray_min_gene_num[num_gene_sets_min_num_genes][0] = "TAIL"
    num_genes_in_gs_min_gene_num = num_genes_gs[find_gs_min_num_genes]
    unique_resources = gs_GeneIDs_cellarray_min_gene_num[0][0]

    # (8) PRINT INTO LOG FILE:
    logger.log('The input file name of gene sets analyzed is: {0}\n', GeneSet_db_file_name)
    logger.log('Number of gene sets with at least {0} genes is {1}.\n', min_gs_size, num_gene_sets_min_num_genes)
    logger.log('Gene sets were taken from the following resources:\n')
    logger.log('{0}\t', unique_resources)
    if print_gene_scores == 1:
        logger.log('Gene scores for the inputed set of gene-sets were printed into files with the suffix .geneset.\n')
    logger.log('\n')

    if run_GSEA == 1:
        logger.log('The GSEA results file and other output files are saved in the directory: {0}\n', Output_dir)
        logger.log('The GSEA results file name is: {0}\n', GSEA_results_file_name)
        logger.log('Here is a key of the header subtitles for each columns:\n')
        logger.log('DB=databse.\n')
        logger.log('GS=gene set.\n')
        logger.log('ORIG_GS_SIZE=Original number of genes per gene set.\n')
        logger.log('EFF_GS_SIZE=Effective number of genes per gene set analyzed by GSEA, after removing genes that were not assigned a gene score (e.g. no SNPs in their region), or\n')
        logger.log('after adjusting for physcial clustering of genes in a given gene set (removing all but one gene from a subset of genes assigned the same best SNP,\n')
        logger.log('keeping the gene with the most significant gene score).\n')
        logger.log('#_GENES_ABS_LIST = number of genes absent from the full human genes list used for the analysis.\n')
        logger.log('#_GENES_WO_SCORE = number of genes that did not get assigned a score since they had no SNPs with association statistics in their extended gene region boundaries.\n')
        logger.log('#_GENES_REM_PROXIM = Number of genes removed due to physical proximity adjustment (See above).\n')
        logger.log('MED_GENE_SIZE_KB = median gene size of all genes in gene set in kilobase units.\n')
        logger.log('MEAN_GENE_SIZE_KB = mean gene size of all genes in gene set in kilobase units.\n')
        logger.log('NOMINAL_GSEA_PVAL_95PERC_CUTOFF = GSEA p-value using 95 percentile of all gene scores for the enrichment cutoff.\n')
        logger.log('FDR_95PERC_CUTOFF = Estimated false discovery rate (q-value) using 95 percentile cutoff.\n')
        logger.log('EXP_#_GENES_ABOVE_95PERC_CUTOFF = Expected number of genes with a corrected gene p-value above the 95 percentile enrichment cutoff.\n')
        logger.log('OBS_#_GENES_ABOVE_95PERC_CUTOFF = Observed number of genes with a corrected gene p-value above the 95 percentile enrichment cutoff.\n')
        logger.log('NOMINAL_GSEA_PVAL_90PERC_CUTOFF = GSEA p-value using 90 percentile of all gene scores for the enrichment cutoff.\n')
        logger.log('FDR_90PERC_CUTOFF =  Estimated false discovery rate (q-value) using 90 percentile cutoff.\n')
        logger.log('EXP_#_GENES_ABOVE_90PERC_CUTOFF = Expected number of genes with a corrected gene p-value above the 90 percentile enrichment cutoff.\n')
        logger.log('OBS_#_GENES_ABOVE_90PERC_CUTOFF = Observed number of genes with a corrected gene p-value above the 90 percentile enrichment cutoff.\n')
        logger.log('#_GENES_FLAGGED = Number of genes in gene set that belong to the subset of genes flagged by the user (e.g. genes near validated disease/trait SNPs).\n')
        logger.log('FLAGGED_GENE_NAMES = Gene symbols of genes in list of flagged genes list that belong to each gene set.\n')
    else:
        logger.log('GSEA was not run on the input gene sets, as requested by the user.\n')

    # (9) LOOP OVER ALL GENE SETS

    EnrichPval_multipleCutoffs_all_genesets = np.zeros((0, len(top_percen_cutoffs)))
    Presence_flag_gs_genes = np.zeros(num_gene_sets_min_num_genes)
    Eff_num_genes_gs = np.zeros(num_gene_sets_min_num_genes)


    Record_results_all_genesets = np.zeros((num_gene_sets_min_num_genes, len(top_percen_cutoffs) + 10))
    Rand_gs_fract_p95 = np.zeros((num_gs_simul, num_gene_sets_min_num_genes))
    Rand_gs_fract_p75 = np.zeros((num_gs_simul, num_gene_sets_min_num_genes))


    first_db = gs_GeneIDs_cellarray_min_gene_num[0][0]
    gs_ind = 0

    db_name = np.zeros(num_gene_sets_min_num_genes, dtype='a{0}'.format(MAXIDLEN))
    GeneSetLabel_array = np.zeros(num_gene_sets_min_num_genes, dtype="a{0}".format(MAXIDLEN))
    GeneSetGeneNames_cell_of_arrays = np.zeros(num_gene_sets_min_num_genes, dtype=",".join(["a{0}".format(MAXIDLEN)] * MAXGENENUM_PER_GENESET))
    Original_gs_size = np.zeros(num_gene_sets_min_num_genes)
    Num_genes_no_match_list = np.zeros(num_gene_sets_min_num_genes)
    num_genes_gs_wo_score = np.zeros(num_gene_sets_min_num_genes)
    num_genes_removed_overlap = np.zeros(num_gene_sets_min_num_genes)
    Med_genesize_gs = np.zeros(num_gene_sets_min_num_genes)
    Mean_genesize_gs = np.zeros(num_gene_sets_min_num_genes)
    Num_LargeGenes_gs = np.zeros(num_gene_sets_min_num_genes)
    Num_SmallGenes_gs = np.zeros(num_gene_sets_min_num_genes)

    obs_fract_above_cutoff_p95 = np.zeros(num_gene_sets_min_num_genes)
    obs_fract_above_cutoff_p75 = np.zeros(num_gene_sets_min_num_genes)
    obs_num_gs_above_cutoff_p95 = np.zeros(num_gene_sets_min_num_genes)
    obs_num_gs_above_cutoff_p75 = np.zeros(num_gene_sets_min_num_genes)

    for gs in range(num_gene_sets_min_num_genes):
        db_name[gs_ind] = gs_GeneIDs_cellarray_min_gene_num[gs][0]
        GeneSetLabel = gs_GeneIDs_cellarray_min_gene_num[gs][1]
        GeneSetLabel_array[gs_ind] = GeneSetLabel

        GeneSetGeneNames = gs_GeneIDs_cellarray_min_gene_num[gs][2]
        GeneSetGeneNames = GeneSetGeneNames[np.greater_equal(GeneSetGeneNames, 0)]

        # % (9.1) Make sure gene set ID/name vector is unique
        GeneSetGeneNames = np.unique(GeneSetGeneNames)
        GeneSetGeneNames_cell_of_arrays[gs_ind] = GeneSetGeneNames
        Original_gs_size[gs_ind] = len(GeneSetGeneNames)

        # (9.2) Find genes that have a match with our whole human genes list
        GS_GeneNames_match = np.intersect1d(AllGene_Names, GeneSetGeneNames)
        Num_genes_no_match_list[gs_ind] = Original_gs_size[gs_ind] - len(GS_GeneNames_match)

        # (9.3) Remove genes without a corrected score and without the flagged gene subset if requested
        find_gs_genes_isan_sorted = np.in1d(AllGene_Names_isan, GS_GeneNames_match)
        GS_GeneNames_match_isan = AllGene_Names_isan[find_gs_genes_isan_sorted]
        find_gs_genes_w_score_sorted = np.in1d(AllGene_Names_score_isan, GS_GeneNames_match)
        GS_GeneNames_match_w_score = AllGene_Names_score_isan[find_gs_genes_w_score_sorted]
        num_genes_gs_wo_score[gs_ind] = len(GS_GeneNames_match) - len(GS_GeneNames_match_w_score)

        find_gs_genes_isan = np.sort(np.where(find_gs_genes_isan_sorted)[0])
        Uncorr_score_gs_isan = Uncorr_score_isan[find_gs_genes_isan, :]
        Corr_score_gs_isan = Corr_score_isan[find_gs_genes_isan]
        GS_GeneNames_match_isan = AllGene_Names_isan[find_gs_genes_isan]


        # (9.4) For subsets of genes that share the same best-SNP-per-gene score, choose the gene with the best score, remove the rest

        if choose_unique_genes == 1:
            Uncorr_score_geneset_pval_sorted_ind = np.argsort(np.abs(Corr_score_gs_isan))
            Uncorr_score_SNP_Chrpos_geneset_sorted = Uncorr_score_gs_isan[Uncorr_score_geneset_pval_sorted_ind, :2]
            # Unique_genes = np.ascontiguousarray(Uncorr_score_SNP_Chrpos_geneset_sorted).view(np.dtype((np.void, Uncorr_score_SNP_Chrpos_geneset_sorted.dtype.itemsize * Uncorr_score_SNP_Chrpos_geneset_sorted.shape[1])))
            temp = Uncorr_score_SNP_Chrpos_geneset_sorted[:, 0] * 1e11 + Uncorr_score_SNP_Chrpos_geneset_sorted[:, 1]
            _, Unique_genes_ind = np.unique(temp, return_index=True)
            num_genes_removed_overlap[gs_ind] = len(GS_GeneNames_match_isan) - len(Unique_genes_ind)
            Eff_num_genes_gs[gs_ind] = len(Unique_genes_ind)

            s1 = GS_GeneNames_match_isan[Uncorr_score_geneset_pval_sorted_ind]
            GS_GeneNames_match_isan_no_overlap = s1[Unique_genes_ind]
            find_gene_set_ind_sorted = np.where(np.in1d(AllGene_Names_isan, GS_GeneNames_match_isan_no_overlap))[0]

            find_gene_set_ind = np.sort(find_gene_set_ind_sorted)
            GS_GeneNames_final = GS_GeneNames_match_isan_no_overlap
        else:
            # num_genes_removed_overlap = 0
            Eff_num_genes_gs[gs_ind] = len(find_gs_genes_isan)

            GS_GeneNames_final = GS_GeneNames_match_isan
            find_gene_set_ind = find_gs_genes_isan

        GS_genesize_vect = AllGeneSizes_isan[find_gene_set_ind]

        Med_genesize_gs[gs_ind] = np.median(GS_genesize_vect[np.logical_not(np.isnan(GS_genesize_vect))]) / 1000.0
        Mean_genesize_gs[gs_ind] = np.mean(GS_genesize_vect[np.logical_not(np.isnan(GS_genesize_vect))]) / 1000.0

        Num_LargeGenes_gs[gs_ind] = len(np.greater_equal(GS_genesize_vect, 100000))
        Num_SmallGenes_gs[gs_ind] = len(np.less_equal(GS_genesize_vect, 10000))

        # (9.5) Record presence of genes from input subset of genes in given gene set (e.g. known disease genes)

        if Flag_gene_set_file_name and Flag_gs == 1:
            if match_genes_GeneID == 1:
                Flag_gs_GeneNames = pd.read_csv(Flag_gene_set_file_name, header=None, delim_whitespace=True).values
            q = np.intersect1d(GS_GeneNames_match, np.unique(Flag_gs_GeneNames))
            if len(q) > 0:
                Presence_flag_gs_genes[gs_ind] = len(q)
        else:
            Presence_flag_gs_genes[gs_ind] = 0


        # (9.6) Print gene scores for each gene set

        if print_gene_scores == 1:
            output_genescore_file_name = get_valid_filename('GeneAssocScores_{0}_{1}_{2}kb_upstr_{3}kb_downstr_{4}.geneset'.format(exp_label, GeneSetLabel, Gene_boundr_upstr/1000, Gene_boundr_downstr/1000, todays_date), Output_dir)
            geneset = Writer(output_genescore_file_name)

            geneset.write('Database\tGene_Set\tGene_Symbol\tEntrez_ID\tGene_p-value\tGene_Chr_Num\tGene_Start_Pos\tGene_End_Pos\tGene_Size_kb\tNum_SNPs_per_Gene\tNum_Indep_HapMap_SNPs_per_Gene\tNum_RecombHotspots_per_Gene\t')

            if SNP_rs_num_file_name and print_rs_num == 1:
                geneset.write('Best_SNP_rs\t')
                if num_col_OR == 1:
                    geneset.write('Best_SNP_rs\t')
            geneset.write('Best_SNP_Chr_Num\tBest_SNP_Chr_Pos\tBest_SNP_Z\tBest_SNP_pval\t')

            if print_best_SNP_OR == 1:
                num_row_OR, num_col_OR = GWAS_SNP_OR_L95_U95.shape
                if num_col_OR == 1:
                    geneset.write('Best_SNP_ODDS_RATIO_OR_EFFECT\t')
                elif num_col_OR > 1:
                    geneset.write('Best_SNP_ODDS_RATIO_OR_EFFECT\tBest_SNP_L95_CI\tBest_SNP_U95_CI\t')

            geneset.write('1=Flagged_gene\n')

            for i in range(len(GeneSetGeneNames)):
                find_gene_a = np.in1d(AllGene_Names, GeneSetGeneNames[i])
                MatchedGene = AllGene_Names[find_gene_a]
                if np.sum(find_gene_a) == 1:
                    find_gene_a = np.where(find_gene_a)[0][0]
                    geneset.write('\t'.join((db_name[gs_ind],
                                             GeneSetLabel_array[gs_ind],
                                             AllRefFlatGeneNames[find_gene_a],
                                             "{0}".format(AllRefFlatGeneID[find_gene_a]),
                                             "NaN" if np.isnan(Corr_score[find_gene_a]) else "{0:E}".format(Corr_score[find_gene_a]),
                                             '\t'.join(map(lambda x: str(x), AllGeneChrPos[find_gene_a, :3])),
                                             "{0:1.0f}".format(round(AllGeneSizes[find_gene_a] / 1000.0, 0)),
                                             "{0}".format(NumSNPs_per_gene[find_gene_a]),
                                             "{0}".format(Num_Indep_SNPs_per_gene[find_gene_a]),
                                             "{0}".format(NumHotspots_per_gene[find_gene_a]))) + '\t')

                    if SNP_rs_num_file_name and print_rs_num == 1 and calculate_Scores_confounders == 1:
                        geneset.write('{0}\t'.format(Best_SNP_rs[find_gene_a]))

                    geneset.write(format_NaN('{0:1.0f}\t', Uncorr_score[find_gene_a, 0]))
                    geneset.write(format_NaN('{0:1.0f}\t', Uncorr_score[find_gene_a, 1]))
                    geneset.write(format_NaN('{0:5.4f}\t', Uncorr_score[find_gene_a, 2]))
                    geneset.write(format_NaN('{0:E}\t', Uncorr_score[find_gene_a, 3]))

                    if print_best_SNP_OR == 1:
                        num_row_OR, num_col_OR = GWAS_SNP_OR_L95_U95.shape
                        if np.logical_not(np.isnan(Best_SNP_rs[find_gene_a])):
                            find_best_SNP_rs = np.where(np.equal(AllSNP_rs_num, Best_SNP_rs[find_gene_a]))[0]
                            for odds in range(num_col_OR):
                                geneset.write(format_NaN('{0:5.4f}\t', GWAS_SNP_OR_L95_U95[find_best_SNP_rs, odds]))
                        else:
                            geneset.write('\t'.join(("NaN") * num_col_OR) + '\t')

                    if Flag_gene_set_file_name and Flag_gs == 1:
                        if len(np.intersect1d(GeneSetGeneNames[i], Flag_gs_GeneNames)):
                            geneset.write('1\n')
                        else:
                            geneset.write('0\n')
                    else:
                        geneset.write('\n')

        # (9.7) Apply GSEA-GWAS to given gene set



        if run_GSEA == 1:
            if GSEA_method == 1:
                # start_time = time.time()
                EnrichPval_multipleCutoffs, obs_num_gs_above_cutoff, obs_fract_above_cutoff, record_rand_gs_num, record_rand_gs_fract = GSEA_GWAS(Uncorr_score_isan, Corr_score_isan, top_percen_cutoffs, num_gs_simul, find_gene_set_ind, score_signif_direct, choose_unique_genes, num_perm_limit)
                # end_time = time.time()

                # logger.log('It took {:.4f} seconds or {:.4f} minutes to run GSEA-GWAS for {}.\n', end_time-start_time, (end_time-start_time)/60, GeneSetLabel)

                final_num_simul, cut = record_rand_gs_fract.shape
                Rand_gs_fract_p95[:final_num_simul, gs_ind] = record_rand_gs_fract[:, 0]
                Rand_gs_fract_p75[:final_num_simul, gs_ind] = record_rand_gs_fract[:, 1]
                obs_fract_above_cutoff_p95[gs_ind] = obs_fract_above_cutoff[0]
                obs_fract_above_cutoff_p75[gs_ind] = obs_fract_above_cutoff[1]
                obs_num_gs_above_cutoff_p95[gs_ind] = obs_num_gs_above_cutoff[0]
                obs_num_gs_above_cutoff_p75[gs_ind] = obs_num_gs_above_cutoff[1]

            elif GSEA_method == 4:
                EnrichPval_multipleCutoffs = GSEA_GWAS_RankSum(Uncorr_score_isan, Corr_score_isan, top_percen_cutoffs, num_gs_simul, find_gene_set_ind, score_signif_direct, choose_unique_genes)
                top_percen_cutoffs = EnrichPval_multipleCutoffs

            EnrichPval_multipleCutoffs_all_genesets = np.vstack((EnrichPval_multipleCutoffs_all_genesets, EnrichPval_multipleCutoffs))

            Record_results_all_genesets[gs, :10] = np.array((Original_gs_size[gs_ind],
                                                             Eff_num_genes_gs[gs_ind],
                                                             Num_genes_no_match_list[gs_ind],
                                                             num_genes_gs_wo_score[gs_ind],
                                                             num_genes_removed_overlap[gs_ind],
                                                             Presence_flag_gs_genes[gs_ind],
                                                             Med_genesize_gs[gs_ind],
                                                             Mean_genesize_gs[gs_ind],
                                                             Num_LargeGenes_gs[gs_ind],
                                                             Num_SmallGenes_gs[gs_ind]))

            if len(EnrichPval_multipleCutoffs):
                if GSEA_method == 1:
                    Record_results_all_genesets[gs, 10:(len(top_percen_cutoffs) + 10)] = EnrichPval_multipleCutoffs
                elif GSEA_method == 4:
                    Record_results_all_genesets[gs, 10] = EnrichPval_multipleCutoffs[0]

            if calculate_FDR == 1:
                if gs_GeneIDs_cellarray_min_gene_num[gs][0] != gs_GeneIDs_cellarray_min_gene_num[gs + 1][0] or gs == (num_gene_sets_min_num_genes - 1):

                    gs_ind += 1
                    score_cutoffs = np.percentile(Corr_score_isan, top_percen_cutoffs)

                    Norm_rand_gs_fract_p95 = np.zeros((num_gs_simul, gs_ind))
                    Norm_rand_gs_fract_p75 = np.zeros((num_gs_simul, gs_ind))
                    Norm_obs_fract_above_cutoff_p95 = np.zeros(gs_ind)
                    Norm_obs_fract_above_cutoff_p75 = np.zeros(gs_ind)

                    FDR_gs_p95 = np.zeros(gs_ind)
                    FDR_gs_p75 = np.zeros(gs_ind)

                    for each_gs in range(gs_ind):
                        # sp.stats.mstats.zscore
                        Norm_rand_gs_fract_p95[:, each_gs] = (Rand_gs_fract_p95[:, each_gs] - np.mean(Rand_gs_fract_p95[:, each_gs])) / np.std(Rand_gs_fract_p95[:, each_gs])
                        Norm_obs_fract_above_cutoff_p95[each_gs] = (obs_fract_above_cutoff_p95[each_gs] - np.mean(Rand_gs_fract_p95[:, each_gs])) / np.std(Rand_gs_fract_p95[:, each_gs])

                        Norm_rand_gs_fract_p75[:, each_gs] = (Rand_gs_fract_p75[:, each_gs] - np.mean(Rand_gs_fract_p75[:, each_gs])) / np.std(Rand_gs_fract_p75[:, each_gs])
                        Norm_obs_fract_above_cutoff_p75[each_gs] = (obs_fract_above_cutoff_p75[each_gs] - np.mean(Rand_gs_fract_p75[:, each_gs])) / np.std(Rand_gs_fract_p75[:, each_gs])

                    n_simul, n_gs = Norm_rand_gs_fract_p95.shape
                    num_rand_gs = n_simul * n_gs

                    for each_gs in range(gs_ind):
                        try:
                            FDR_gs_p95[each_gs] = (1.0 * np.sum(np.greater_equal(Norm_rand_gs_fract_p95, Norm_obs_fract_above_cutoff_p95[each_gs])) / num_rand_gs) / (1.0 * np.sum(np.greater_equal(Norm_obs_fract_above_cutoff_p95, Norm_obs_fract_above_cutoff_p95[each_gs])) / len(Norm_obs_fract_above_cutoff_p95))
                        except ZeroDivisionError:
                            FDR_gs_p95[each_gs] = np.nan
                        try:
                            FDR_gs_p75[each_gs] = (1.0 * np.sum(np.greater_equal(Norm_rand_gs_fract_p75, Norm_obs_fract_above_cutoff_p75[each_gs])) / num_rand_gs) / (1.0 * np.sum(np.greater_equal(Norm_obs_fract_above_cutoff_p75, Norm_obs_fract_above_cutoff_p75[each_gs])) / len(Norm_obs_fract_above_cutoff_p75))
                        except ZeroDivisionError:
                            FDR_gs_p75[each_gs] = np.nan

                        if FDR_gs_p95[each_gs] > 1:
                            FDR_gs_p95[each_gs] = 1
                        if FDR_gs_p75[each_gs] > 1:
                            FDR_gs_p75[each_gs] = 1

                        results.write('\t'.join((db_name[each_gs],
                                                 GeneSetLabel_array[each_gs],
                                                 "{0:1.0f}".format(Original_gs_size[each_gs]),
                                                 "{0:1.0f}".format(Eff_num_genes_gs[each_gs]),
                                                 "{0:1.0f}".format(Num_genes_no_match_list[each_gs]),
                                                 "{0:1.0f}".format(num_genes_gs_wo_score[each_gs]),
                                                 "{0:1.0f}".format(num_genes_removed_overlap[each_gs]),
                                                 "{0:1.0f}".format(Med_genesize_gs[each_gs]),
                                                 "{0:1.0f}".format(Mean_genesize_gs[each_gs]))))

                        if len(EnrichPval_multipleCutoffs):
                            if GSEA_method == 1:
                                results.write('\t' +
                                              '\t'.join(('{0:E}'.format(EnrichPval_multipleCutoffs_all_genesets[each_gs, 0]),
                                                         format_NaN('{0:E}', FDR_gs_p95[each_gs], sep=''),
                                                         '{0:1.0f}'.format(round(top_percen_cutoffs[0] / 100.0 * Eff_num_genes_gs[each_gs])),
                                                         '{0:1.0f}'.format(obs_num_gs_above_cutoff_p95[each_gs]),
                                                         '{0:E}'.format(EnrichPval_multipleCutoffs_all_genesets[each_gs, 1]),
                                                         format_NaN('{0:E}', FDR_gs_p75[each_gs], sep=''),
                                                         '{0:1.0f}'.format(round(top_percen_cutoffs[0] / 100.0 * Eff_num_genes_gs[each_gs])),
                                                         '{0:1.0f}'.format(obs_num_gs_above_cutoff_p75[each_gs]))))
                            elif GSEA_method == 4:
                                results.write('\t{0:E}', EnrichPval_multipleCutoffs_all_genesets[each_gs, 0])

                        results.write('\t{0:1.0f}', Presence_flag_gs_genes[each_gs])

                        if Flag_gene_set_file_name and Flag_gs == 1 and print_trait_gene_names_in_gs == 1:
                            Flagged_GeneNames_in_GeneSet = np.interest1d(GeneSetGeneNames_cell_of_arrays[each_gs], Flag_gs_GeneNames)
                            if len(Flagged_GeneNames_in_GeneSet):
                                results.write('\t')
                                for u in range(len(Flagged_GeneNames_in_GeneSet) - 1):
                                    find_flag_gene_name = np.where(np.in1d(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet[u]))[0]
                                    if np.sum(find_flag_gene_name):
                                        Flag_gene_name = AllRefFlatGeneNames[find_flag_gene_name[0]]
                                        results.write("{0},", Flag_gene_name)

                                find_flag_gene_name = np.where(np.in1d(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet[len(Flagged_GeneNames_in_GeneSet)]))[0]
                                if np.sum(find_flag_gene_name):
                                    Flag_gene_name = AllRefFlatGeneNames[find_flag_gene_name[0]]
                                    results.write(Flag_gene_name)
                        results.write('\n')

                    gs_ind = -1
                    EnrichPval_multipleCutoffs_all_genesets = np.zeros((0, len(top_percen_cutoffs)))
            else:
                # (9.8) Print Gene Set Enrichment results for each gene set into output file
                results.write('\t'.join((db_name[each_gs],
                                         GeneSetLabel_array[each_gs],
                                         "{0:1.0f}".format(Original_gs_size[each_gs]),
                                         "{0:1.0f}".format(Eff_num_genes_gs[each_gs]),
                                         "{0:1.0f}".format(Num_genes_no_match_list[each_gs]),
                                         "{0:1.0f}".format(num_genes_gs_wo_score[each_gs]),
                                         "{0:1.0f}".format(num_genes_removed_overlap[each_gs]),
                                         "{0:1.0f}".format(Med_genesize_gs[each_gs]),
                                         "{0:1.0f}".format(Mean_genesize_gs[each_gs]))))

                if len(EnrichPval_multipleCutoffs):
                    if GSEA_method == 1:
                        for y in range(len(EnrichPval_multipleCutoffs)):
                            results.write('\t' +
                                          '\t'.join(('{0:E}'.format(EnrichPval_multipleCutoffs[y]),
                                                     '{0:1.0f}'.format(round(top_percen_cutoffs[y] / 100.0 * Eff_num_genes_gs[gs_ind])),
                                                     '{0:1.0f}'.format(obs_num_gs_above_cutoff[y]))))
                    elif GSEA_method == 4:
                        results.write('\t{0:E}', EnrichPval_multipleCutoffs[0])

                results.write('\t{0:1.0f}', Presence_flag_gs_genes[gs_ind])

                if Flag_gene_set_file_name and Flag_gs == 1 and print_trait_gene_names_in_gs == 1:
                    Flagged_GeneNames_in_GeneSet = np.interest1d(GeneSetGeneNames_cell_of_arrays[gs_ind], Flag_gs_GeneNames)
                    if len(Flagged_GeneNames_in_GeneSet):
                        results.write('\t')
                        for u in range(len(Flagged_GeneNames_in_GeneSet) - 1):
                            find_flag_gene_name = np.where(np.in1d(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet[u]))[0]
                            if np.sum(find_flag_gene_name):
                                Flag_gene_name = AllRefFlatGeneNames[find_flag_gene_name[0]]
                                results.write("{0},", Flag_gene_name)

                        find_flag_gene_name = np.where(np.in1d(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet[len(Flagged_GeneNames_in_GeneSet)]))[0]
                        if np.sum(find_flag_gene_name):
                            Flag_gene_name = AllRefFlatGeneNames[find_flag_gene_name[0]]
                            results.write(Flag_gene_name)
                results.write('\n')
        gs_ind += 1

    if run_GSEA == 1:
        # save_com = ['save ', num2str(Output_dir) ,'/GeneSetEnrichPval_AllGeneSets_', num2str(exp_label) , '_',  num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(num_gs_simul), 'perm_', num2str(todays_date) , ' EnrichPval_multipleCutoffs_all_genesets'];
        pass

    run_time = time.time() - analysis_start_time
    logger.log("Time it took the program to run: {0} minutes = {1} hours = {2} days\n", run_time / 60, run_time / 3600, run_time / (3600 * 24))

if __name__ == '__main__':
    magenta()
