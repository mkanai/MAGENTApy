import numpy as np
import scipy as sp
import scipy.stats


def Converting_GWAS_TwoTailedPval_to_Zscores(GWAS_SNPChrNumPos_Pval):
    pval = GWAS_SNPChrNumPos_Pval[:, 2]
    zscores = np.sqrt(sp.stats.chi2.ppf(1 - pval, 1))
    GWAS_SNPChrNumPos_ZscorePval = np.concatenate((GWAS_SNPChrNumPos_Pval[:, :2], zscores[:, None], GWAS_SNPChrNumPos_Pval[:, 2:None]), axis=1)
    return GWAS_SNPChrNumPos_ZscorePval
