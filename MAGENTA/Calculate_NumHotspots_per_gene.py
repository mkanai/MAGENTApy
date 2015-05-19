import numpy as np
from utils import *


def Calculate_NumHotspots_per_gene(GeneChrPos, hotspotPos, Gene_boundr_upstr, Gene_boundr_downstr, StrandOrient):
    num_genes, x = GeneChrPos.shape
    NumHotspots_per_gene = np.zeros(num_genes, dtype="i")
    for i in range(num_genes):
        if StrandOrient[i] == 1:
            find_num_hotspots = np.logical_and(np.equal(hotspotPos[:, 0], GeneChrPos[i, 0]),
                                               np.logical_or(np.logical_and(np.greater_equal(hotspotPos[:, 1], (GeneChrPos[i, 1] - Gene_boundr_upstr)),
                                                                            np.less_equal(hotspotPos[:, 1], (GeneChrPos[i, 2] + Gene_boundr_downstr))),
                                                             np.logical_and(np.greater_equal(hotspotPos[:, 2], (GeneChrPos[i, 1] - Gene_boundr_upstr)),
                                                                            np.less_equal(hotspotPos[:, 2], (GeneChrPos[i, 2] + Gene_boundr_downstr)))))
        else:
            find_num_hotspots = np.logical_and(np.equal(hotspotPos[:, 0], GeneChrPos[i, 0]),
                                               np.logical_or(np.logical_and(np.greater_equal(hotspotPos[:, 1], (GeneChrPos[i, 1] - Gene_boundr_downstr)),
                                                                            np.less_equal(hotspotPos[:, 1], (GeneChrPos[i, 2] + Gene_boundr_upstr))),
                                                             np.logical_and(np.greater_equal(hotspotPos[:, 2], (GeneChrPos[i, 1] - Gene_boundr_downstr)),
                                                                            np.less_equal(hotspotPos[:, 2], (GeneChrPos[i, 2] + Gene_boundr_upstr)))))
        if np.any(find_num_hotspots):
            NumHotspots_per_gene[i] = np.sum(find_num_hotspots)
    return NumHotspots_per_gene
