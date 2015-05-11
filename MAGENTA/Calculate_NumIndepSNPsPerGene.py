import numpy as np
from utils import *


def Calculate_NumIndepSNPsPerGene(AllHumanGeneChrPos, PrunedSNPsChrPos, boundr_upstr, boundr_downstr):
    n_r, n_c = AllHumanGeneChrPos.shape
    Num_Indep_SNPs_per_gene = np.zeros(n_r, dtype="i")
    strand = np.all(np.equal(AllHumanGeneChrPos[:, 5], 1))

    for i in range(n_r):
        if strand:
            Num_Indep_SNPs_per_gene[i] = np.sum(logical_and(np.equal(PrunedSNPsChrPos[:, 0], AllHumanGeneChrPos[i, 0]),
                                                            np.greater_equal(PrunedSNPsChrPos[:, 1], (AllHumanGeneChrPos[i, 1] - boundr_upstr)),
                                                            np.less_equal(PrunedSNPsChrPos[:, 1], (AllHumanGeneChrPos[i, 2] + boundr_downstr))))
        else:
            Num_Indep_SNPs_per_gene[i] = np.sum(logical_and(np.equal(PrunedSNPsChrPos[:, 0], AllHumanGeneChrPos[i, 0]),
                                                            np.greater_equal(PrunedSNPsChrPos[:, 1], (AllHumanGeneChrPos[i, 1] - boundr_downstr)),
                                                            np.less_equal(PrunedSNPsChrPos[:, 1], (AllHumanGeneChrPos[i, 2] + boundr_upstr))))
    return Num_Indep_SNPs_per_gene
