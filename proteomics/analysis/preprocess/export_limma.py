from pathlib import Path
from typing import NamedTuple

import pandas as pd

from proteomics.analysis.io.load_metadata import Contrast

__all__ = ["make_limma_contrasts", "LimmaInputs"]


def subset_genes(
    gene_list_file: Path | None,
    gene_list_col: str | None,
    genes_sheet: str | None,
) -> bool:
    return not (
        gene_list_file is None or gene_list_col is None or genes_sheet is None
    )


class LimmaInputs(NamedTuple):
    """
    Inputs into a function that runs limma.
    """
    counts_file: Path
    metadata_file: Path
    contrast_name: str
    contrast_1: str
    contrast_2: str


def make_limma_contrasts(
    df: pd.DataFrame,
    *,
    output_dir: Path,
    contrast_list: list[Contrast],
    gene_list_file: Path | None,
) -> list[LimmaInputs]:
    if not output_dir.exists():
        output_dir.mkdir()
        print(f"Created directory: {output_dir}")

    limma_inputs = []
    for contrast in contrast_list:
        g1, g2 = contrast["contrast"]

        counts_df = df[contrast["group_0_cols"] + contrast["group_1_cols"]]

        if subset_genes(
            gene_list_file=gene_list_file,
            gene_list_col=contrast["gene_list_col"],
            genes_sheet=contrast["genes_sheet"],
        ):
            print("Subsetting genes")
            genes_df = pd.read_excel(
                gene_list_file, sheet_name=contrast["genes_sheet"]
            )[contrast["gene_list_col"]]
            print(
                f"Counts df genes: {counts_df.shape[0]}, Gene list: {genes_df.shape[0]}"
            )
            counts_df = counts_df.loc[genes_df]
            print(f"Counts df genes after subsetting: {counts_df.shape[0]}")

        counts_df.to_csv(output_dir / f"{g1}_vs_{g2}_counts.csv")

        metadata = pd.DataFrame(
            {
                "Sample": counts_df.columns,
                "Group": [
                    g1 if col in contrast["group_0_cols"] else g2
                    for col in counts_df.columns
                ],
            }
        )

        metadata.to_csv(output_dir / f"{g1}_vs_{g2}_metadata.csv", index=False)

        limma_inputs.append(
            LimmaInputs(
                counts_file=output_dir / f"{g1}_vs_{g2}_counts.csv",
                metadata_file=output_dir / f"{g1}_vs_{g2}_metadata.csv",
                contrast_name=f"{g1}_vs_{g2}",
                contrast_1=g1,
                contrast_2=g2,
            )
        )

    return limma_inputs
