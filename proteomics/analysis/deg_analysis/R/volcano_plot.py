import re
from pathlib import Path

from proteomics.analysis.deg_analysis.base_args import (
    RDegPlotArgs,
    RConfig,
    RunRDegAnalysis,
)


class VolcanoArgs(RDegPlotArgs):
    upregulated_color: str = "red"
    downregulated_color: str = "royalblue"
    not_sig_color: str = "black"


class VolcanoPlot(VolcanoArgs, RunRDegAnalysis[Path]):
    def get_r_script(self) -> Path:
        return Path(__file__).parent / "volcano-plot.R"

    def process_stdout(self, stdout: str) -> Path | None:
        try:
            cleaned_string = re.sub(r"^\[\d+]\s*", "", stdout[-1].strip())
            volcano_png = Path(cleaned_string.strip('"'))
            # make sure the file exists
            if not volcano_png.exists() or not volcano_png.suffix == ".png":
                raise FileNotFoundError(
                    f"File {volcano_png} is not a png file"
                )
            # make sure its a png file
        except FileNotFoundError as e:
            print(f"Error: {e}")
            volcano_png = None

        return volcano_png


def run_volcano_plot_r(
    *,
    r_config: RConfig,
    output_dir: Path,
    deg_results: Path,
    experiment: str,
    volcano_args: VolcanoArgs,
) -> Path | None:
    volcano_plot = VolcanoPlot(
        output_dir=output_dir,
        input_file=deg_results,
        rscript_bin=r_config.rscript_bin,
        experiment=experiment,
        **volcano_args.model_dump(),
    )

    return volcano_plot.run_analysis()


def main():
    r_config = RConfig(
        rscript_bin="/Users/ac4294/.pyenv/versions/miniforge3-24.11.3-0/envs/liver_prot_s2808d_2025/bin/Rscript"
    )
    output_dir = Path("/Users/ac4294/Downloads")
    deg_results = Path(
        "/Users/ac4294/dev/2025/yang_liver_proteomics_s2808d/data/example/_results-1742436624597126/S2808D_HFD_vs_WT_HFD/limma_outputs/S2808D_HFD_vs_WT_HFD_deg_limma.csv"
    )
    experiment = "S2808D_HFD_vs_WT_HFD"
    volcano_args = VolcanoArgs(
        gene_column="X",
        lfc_column="logFC",
        pval_column="adj.P.Val",
        fc_threshold=1.5,
        pval_threshold=0.05,
        width=6,
        height=8,
    )

    enrichment_results = run_volcano_plot_r(
        r_config=r_config,
        output_dir=output_dir,
        deg_results=deg_results,
        experiment=experiment,
        volcano_args=volcano_args,
    )

    print("The volcano plot results")
    print(enrichment_results)


if __name__ == "__main__":
    main()
