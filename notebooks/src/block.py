from dataclasses import dataclass
from dataclasses import field

import pandas as pd


@dataclass(frozen=True)
class Block:
    cpg: str
    block_chr: str
    block_start: int
    block_end: int
    cpgs_in_block: list

    @property
    def block_id(self) -> str:
        return f"{self.block_chr}:{self.block_start}-{self.block_end}"

    @property
    def block_size(self) -> int:
        return self.block_end - self.block_start

    @property
    def n_cpgs(self) -> int:
        return len(self.cpgs_in_block)


@dataclass(frozen=True)
class AnnotatedBlock:
    block: Block
    correlated_genes: list
    stats_table: pd.DataFrame

    @property
    def block_chr(self) -> str:
        return self.block.block_chr

    @property
    def block_id(self) -> str:
        return self.block.block_id



@dataclass
class BlockCollector:
    blocks: list[Block] | list[AnnotatedBlock] = field(default_factory=list)

    def extract_block(self, chromosome: str) -> list[Block]:
        selected = []
        for block in self.blocks:
            if block.block_chr == chromosome:
                selected.append(block)

        return selected

    def add_single(self, block: Block | AnnotatedBlock) -> None:
        self.blocks.append(block)

    def add_many(self, blocks: list[Block] | list[AnnotatedBlock]) -> None:
        self.blocks.extend(blocks)




