from pathlib import Path

from pydantic import DirectoryPath, Field, FilePath

from proteomics.analysis.deg_analysis.base_args import (
    DegAnalysisArgs,
    RunRMixin,
    RConfig,
)


class LimmaArgs(DegAnalysisArgs):
    """
    Creates the arguments for running the limma analysis in R.

    - make_option(c("--counts"), type = "character", help = "Path to the counts file"),
    - make_option(c("--metadata"), type = "character", help = "Path to the metadata file"),
    - make_option(c("--output_dir"), type = "character", help = "Output directory for results"),
    - make_option(c("--sig_output_dir"), type = "character", help = "Output directory for significant results"),
    - make_option(c("--contrast_name"), type = "character", help = "Name of the contrast"),
    - make_option(c("--contrast_1"), type = "character", help = "First group for contrast"),
    - make_option(c("--contrast_2"), type = "character", help = "Second group for contrast"),

    """

    counts: FilePath = Field(..., description="The count file")
    metadata: FilePath = Field(..., description="The metadata file")
    output_dir: DirectoryPath = Field(..., description="The output directory")
    sig_output_dir: DirectoryPath = Field(
        ..., description="The output directory for significant genes"
    )
    contrast_name: str = Field(..., description="Name of the contrast")
    contrast_1: str = Field(..., description="First group for contrast")
    contrast_2: str = Field(..., description="Second group for contrast")


class RunLimma(LimmaArgs, RunRMixin[Path]):
    def get_r_script(self) -> Path:
        return Path(__file__).parent / "limma.R"

    def process_stdout(self, stdout: list[str]) -> Path | None:
        try:
            result_path = Path(stdout[-1].strip())
            # make sure the file exists
            assert result_path.exists()
        except Exception as e:
            print(f"Error: {e}")
            result_path = None

        return result_path


def run_limma_r(
    *,
    r_config: RConfig,
    output_dir: Path,
    counts: Path,
    metadata: Path,
    sig_output_dir: Path,
    contrast_name: str,
    contrast_1: str,
    contrast_2: str,
) -> Path | None:
    limma = RunLimma(
        rscript_bin=r_config.rscript_bin,
        counts=counts,
        metadata=metadata,
        output_dir=output_dir,
        sig_output_dir=sig_output_dir,
        contrast_name=contrast_name,
        contrast_1=contrast_1,
        contrast_2=contrast_2,
    )

    return limma.run_analysis()


def main():
    r_config = RConfig(
        rscript_bin="/Users/ac4294/.pyenv/versions/miniforge3-24.11.3-0/envs/liver_prot_s2808d_2025/bin/Rscript"
    )
    output_dir = Path("/Users/ac4294/Downloads")
    counts = Path(
        "/Users/ac4294/dev/2025/yang_liver_proteomics_s2808d/data/example/_results-1742436624597126/limma_inputs/S2808D_vs_WT_counts.csv"
    )
    metadata = Path(
        "/Users/ac4294/dev/2025/yang_liver_proteomics_s2808d/data/example/_results-1742436624597126/limma_inputs/S2808D_vs_WT_metadata.csv"
    )
    sig_output_dir = output_dir
    contrast_name = "S2808D_vs_WT"
    contrast_1 = "S2808D"
    contrast_2 = "WT"

    result = run_limma_r(
        r_config=r_config,
        output_dir=output_dir,
        counts=counts,
        metadata=metadata,
        sig_output_dir=sig_output_dir,
        contrast_name=contrast_name,
        contrast_1=contrast_1,
        contrast_2=contrast_2,
    )

    print(result)


if __name__ == "__main__":
    main()
