import numpy as np
import scipy as sp
import scipy.stats
from stepwisefit import *


def GeneScoreRegressCorr(confounders, GeneScores):
    num_genes, num_confounders = confounders.shape
    cutoff = 0.05
    find_score_isan = np.logical_and(np.logical_not(np.isnan(GeneScores)), np.isfinite(GeneScores))
    beta, se, pval, inmodel, stats, nextstep, history = stepwisefit(confounders[find_score_isan], np.abs(GeneScores[find_score_isan]))
    Residuals = np.abs(GeneScores) - stats.intercept
    for j in range(num_confounders):
        if pval[j] <= cutoff:
            Residuals = Residuals - beta[j]*confounders[:, j]

    find_residuals_isan = np.logical_and(np.logical_not(np.isnan(Residuals)), np.isfinite(Residuals))
    RegressCorrGeneScores_pval = 1 - sp.stats.norm.cdf(Residuals, np.mean(Residuals[find_residuals_isan]), np.std(Residuals[find_residuals_isan]))
    return (RegressCorrGeneScores_pval, Residuals)
