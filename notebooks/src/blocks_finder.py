import pandas as pd
import scipy.stats as sts
from tqdm import tqdm
from .block import Block

def load(data: pd.DataFrame | str) -> pd.DataFrame:
    if isinstance(data, str):
        if data.endswith(".csv"):
            return pd.read_csv(data, index_col=0, dtype_backend="pyarrow")
        return pd.read_parquet(data, dtype_backend="pyarrow")
    return data


class BlocksFinder:
    def __init__(self, methylation: pd.DataFrame | str, manifest: pd.DataFrame | str):
        self.methylation = load(methylation).T
        self.manifest = load(manifest)
        self.manifest = self.manifest.loc[list(set.intersection(set(self.methylation.columns),
                                                                set(self.manifest.index)))]
        self.genome_map = None

    def build_genome_map(self) -> dict:
        genome_map = {}
        for chr_ in tqdm(self.manifest.CHR_hg38.unique()):
            temp = self.manifest[self.manifest.CHR_hg38 == chr_][["Start_hg38", "End_hg38"]]
            temp = temp.sort_values("Start_hg38", ascending=True)
            genome_map[chr_] = temp

        self.genome_map = genome_map

    @staticmethod
    def find_blocks(chr_: str, genome_map: dict, distance_threshold: int = 10000, corr_threshold: float = 0.7, n_threshold: int = 3) -> list[Block]:
        assigned_to_block_cpgs = set()
        chromosome_map = genome_map[chr_]
        identified_blocks = []

        for cpg in tqdm(chromosome_map.index, desc=chr_):
            if cpg not in assigned_to_block_cpgs:
                cpg_position = chromosome_map.loc[cpg, "Start_hg38"]

                distance = chromosome_map["End_hg38"] - cpg_position
                selected = distance[(distance > 0) & (distance <= distance_threshold)]

                if len(selected) < n_threshold:
                    continue

                methylation = self.methylation[selected.index]
                for consecutive_cpg in selected.index:
                    corr, _ = sts.pearsonr(methylation[cpg],
                                           methylation[consecutive_cpg])
                    if abs(corr) < corr_threshold:
                        selected = selected.drop(consecutive_cpg)

                if len(selected) < n_threshold:
                    continue

                block_end_position = cpg_position + selected.max()
                cpgs_in_block = selected.index.tolist()

                single_block = Block(cpg, chr_, cpg_position, block_end_position, cpgs_in_block)
                identified_blocks.append(single_block)
                assigned_to_block_cpgs.add(cpg)

        return identified_blocks

