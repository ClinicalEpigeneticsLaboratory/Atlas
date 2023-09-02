from os import makedirs
from os.path import join
import pandas as pd
from tqdm import tqdm
import scipy.stats as sts


def build_genome_map(manifest: pd.DataFrame, methylation: pd.DataFrame) -> dict:
    cmaps = {}
    cpgs = set(methylation.index)
    manifest = manifest.loc[list(cpgs.intersection(set(manifest.index)))]

    for chr_ in tqdm(manifest.CHR.unique()):
        temp = manifest[manifest.CHR == chr_]
        temp_distance_matrix = temp["MAPINFO"]
        temp_distance_matrix = temp_distance_matrix.sort_values(ascending=True)
        cmaps[chr_] = temp_distance_matrix

    return cmaps


def estimate_associations_within_blocks(cmaps: dict, methylation: pd.DataFrame, distance: int = 10000):
    stats_per_chr = {}
    already_assigned_cpgs = []

    for chr_, cmap in cmaps.items():
        stats = []

        for cpg in tqdm(cmap.index, desc=f"CHR{chr_}"):
            if cpg in already_assigned_cpgs:
                continue

            cpg_met_levels = methylation[cpg]

            consecutive_cpgs_distance = cmap - cmap.loc[cpg]
            consecutive_cpgs_distance = consecutive_cpgs_distance[(consecutive_cpgs_distance > 0) &
                                                                  (consecutive_cpgs_distance < distance)]

            for consecutive_cpg in consecutive_cpgs_distance.index:
                consecutive_cpg_met_levels = methylation[consecutive_cpg]
                corr = sts.pearsonr(cpg_met_levels, consecutive_cpg_met_levels)
                coef, pval = corr.statistic, corr.pvalue

                already_assigned_cpgs.append(consecutive_cpg)
                stats.append({
                    "CpG 1": cpg,
                    "CpG 2": consecutive_cpg,
                    "CpG 1 POS": cmap.loc[cpg],
                    "CpG 2 POS": cmap.loc[consecutive_cpg],
                    "Distance": consecutive_cpgs_distance.loc[consecutive_cpg],
                    "CHR": chr_,
                    "corr-coef": coef,
                    "p-value": pval
                })

        stats_per_chr[chr_] = pd.DataFrame(stats)

    return stats_per_chr


def export_blocks_stats(corr_matrix: dict, path: str):
    output = join(path, "corr")
    makedirs(output, exist_ok=True)
    for chr_, matrix in corr_matrix.items():
        matrix.to_parquet(join(output, f"{chr_}.parquet"))
