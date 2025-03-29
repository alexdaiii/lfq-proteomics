from pathlib import Path
from typing import Literal

import pandas as pd
from pydantic import Field, BaseModel, FilePath
from pydantic_core import from_json

from proteomics.analysis.deg_analysis.base_args import (
    RDegPlotArgs,
    RConfig,
    RunRDegAnalysis,
)


class EnrichmentArgs(RDegPlotArgs):
    organism: Literal["mmu", "hsa"] = Field(
        ...,
        description="The organism to use for the enrichment analysis. Must be one of 'mmu' or 'hsa'",
    )


# type EnrichResult = {
#     file_path: string;
#     ont: "BP" | "CC" | "MF" | "KEGG" | null;
#     result_type: "enrichResult" | "compareResult" | "plot" | "geneIds";
# }
class EnrichResult(BaseModel):
    file_path: FilePath
    ont: Literal["BP", "CC", "MF", "KEGG", "None"]
    result_type: Literal["enrichResult", "compareResult", "plot", "geneIds"]


def parse_json_str(outline: str) -> EnrichResult | None:
    try:
        return EnrichResult.model_validate(from_json(outline))
    except Exception:
        return None


class EnrichmentAnalysis(EnrichmentArgs, RunRDegAnalysis[list[EnrichResult]]):
    def get_r_script(self) -> Path:
        return Path(__file__).parent / "run_enrichment.R"

    def process_stdout(self, stdout: str) -> list[EnrichResult]:
        results = []

        for line in stdout:
            for outline in line.split("\n"):
                enrich_result = parse_json_str(outline)
                if enrich_result:
                    results.append(enrich_result)

        return results


def fix_result_for_excel(result: EnrichResult) -> None:
    """
    Excel messes up the formatting of the gene ratio and bg ratio columns and
    assumes they are dates. This function fixes that by adding a ' to the front
    of the values so that it is treated as a string.
    """
    df = pd.read_csv(result.file_path)

    cols_to_fix = [
        "GeneRatio",
        "BgRatio",
    ]

    for col in cols_to_fix:
        print(f"Fixing column {col}")
        df[col] = df[col].apply(lambda x: f'"{x}"')

    print(f"Writing fixed file to {result.file_path}")
    df.to_csv(result.file_path, index=False)


def run_enrichment_r(
    *,
    r_config: RConfig,
    output_dir: Path,
    deg_results: Path,
    experiment: str,
    enrichment_args: EnrichmentArgs,
) -> list[EnrichResult]:
    enrichment_analysis = EnrichmentAnalysis(
        output_dir=output_dir,
        input_file=deg_results,
        rscript_bin=r_config.rscript_bin,
        experiment=experiment,
        **enrichment_args.model_dump(),
    )

    results: list[EnrichResult] = enrichment_analysis.run_analysis()

    for result in results:
        if (result.result_type == "enrichResult") or (
            result.result_type == "compareResult"
        ):
            fix_result_for_excel(result)

    return results


def main():
    r_config = RConfig(
        rscript_bin="/Users/ac4294/.pyenv/versions/miniforge3-24.11.3-0/envs/liver_prot_s2808d_2025/bin/Rscript"
    )
    output_dir = Path("/Users/ac4294/Downloads")
    deg_results = Path("/Users/ac4294/dev/2025/yang_liver_proteomics_s2808d/"
                       ".metaflow/ProteomicsAnalysis/1743014843985429/"
                       "_results/S2808D_vs_WT/limma_outputs/S2808D_vs_WT_deg_limma.csv")
    experiment = "DEG"

    enrichment_args = EnrichmentArgs(
        organism="mmu",
        gene_column="X",
        lfc_column="logFC",
        pval_column="adj.P.Val",
        fc_threshold=1.5,
        pval_threshold=0.05,
        width=6,
        height=8,
    )

    enrichment_results = run_enrichment_r(
        r_config=r_config,
        output_dir=output_dir,
        deg_results=deg_results,
        experiment=experiment,
        enrichment_args=enrichment_args,
    )

    print("The enrichment results")
    print(enrichment_results)


if __name__ == "__main__":
    main()
