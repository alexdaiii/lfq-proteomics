from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from matplotlib.colors import LinearSegmentedColormap
from metaflow import (
    FlowSpec,
    step,
    current,
    card,  # type: ignore
    IncludeFile,
)
from metaflow.cards import Image, Markdown, Table
from PIL import Image as PILImage

from proteomics.analysis.deg_analysis.R.limma import run_limma_r
from proteomics.analysis.deg_analysis.R.run_enrichment import (
    EnrichmentArgs,
    run_enrichment_r,
)
from proteomics.analysis.deg_analysis.R.volcano_plot import (
    VolcanoArgs,
    run_volcano_plot_r,
)
from proteomics.analysis.deg_analysis.base_args import RConfig
from proteomics.analysis.deg_analysis.fix_kegg_ids import fix_kegg_ids
from proteomics.analysis.deg_analysis.heatmap import (
    make_heatmap_sample,
    MakeHeatmapOtherKwargs,
)
from proteomics.analysis.io.get_raw_data import (
    load_imputed_counts,
    ImputedIntensity,
)
from proteomics.analysis.io.load_metadata import MetadataMaps, make_metadata, \
    validate_metadata, create_contrast_from_metadata
from proteomics.analysis.pca.pca import run_pca
from proteomics.analysis.preprocess import (
    PreprocessParams,
    normalize,
)
from proteomics.analysis.preprocess.export_limma import (
    make_limma_contrasts,
    LimmaInputs,
)
from proteomics.utils.base_params import BaseParams
from proteomics.utils.metaflow_util import get_task_output, get_run_output


class ParameterFile(BaseParams):
    imputed_data: ImputedIntensity
    preprocess: PreprocessParams
    r: RConfig
    heatmap: MakeHeatmapOtherKwargs
    volcano: VolcanoArgs
    enrich: EnrichmentArgs


def config_file_parser(config: str) -> ParameterFile:
    config = yaml.safe_load(config)
    print("Loaded:")
    print(config)
    return ParameterFile(**config)


class ProteomicsAnalysis(FlowSpec):
    params_file: str = IncludeFile(
        "params_file",
        # default=(get_project_root() / "config.yaml").__str__(),
    )

    parameters: ParameterFile

    _run_output_dir: Path
    preprocess_config: PreprocessParams
    raw_input_config: ImputedIntensity
    r_config: RConfig

    raw_counts: pd.DataFrame
    metadata_maps: MetadataMaps

    counts_norm: pd.DataFrame

    pca_df: pd.DataFrame

    limma_inputs: list[LimmaInputs]
    limma_input: LimmaInputs

    result_path: Path | None
    kegg_results: list[Path] | None
    gene_ids: pd.DataFrame | None

    def create_output_dir(self, stub: str, parents: bool = False) -> Path:
        output_dir = self._run_output_dir / stub
        output_dir.mkdir(exist_ok=True, parents=parents)
        print(f"Creating directory: {output_dir}")
        return output_dir

    @card
    @step
    def start(self):
        self.parameters = config_file_parser(self.params_file)

        print("Config loaded and validated")

        self.preprocess_config = self.parameters.preprocess
        self.raw_input_config = self.parameters.imputed_data
        self.r_config = self.parameters.r

        self._run_output_dir = get_run_output(self) / "_results"
        print(f"Creating output directory: {self._run_output_dir}")
        self._run_output_dir.mkdir(exist_ok=True)
        print("Starting analysis")

        self.next(self.load_raw_counts_and_metadata)

    @card
    @step
    def load_raw_counts_and_metadata(self):
        print("Loading raw counts and metadata")

        self.raw_counts = load_imputed_counts(
            input_file=self.raw_input_config.input_file,
            sheet_name=self.raw_input_config.sheet_name,
            index_col=self.raw_input_config.index_col,
        )

        counts_info = (
            f"Loaded proteomics data with {self.raw_counts.shape[0]} genes "
            f"and {self.raw_counts.shape[1]} samples"
        )
        print(counts_info)
        current.card.append(Markdown(counts_info))

        self.metadata_maps = make_metadata(
            metadata_file=Path(self.raw_input_config.metadata_file),
            metadata_sheet_name=self.raw_input_config.metadata_sheet_name,
            metadata_index_col=self.raw_input_config.metadata_index_col,
            condition_col=self.raw_input_config.condition_col,
        )

        validate_metadata(
            self.raw_counts,
            self.metadata_maps,
        )

        self.next(self.normalize_data)

    @card
    @step
    def normalize_data(self):
        plot_dir = self.create_output_dir("normalization")

        print("Normalizing data")
        normalize_rt = normalize(
            df=self.raw_counts,
            plot_dir=plot_dir,
            n_rows=self.preprocess_config.n_rows,
            n_cols=self.preprocess_config.n_cols,
        )

        self.counts_norm = normalize_rt.df_norm
        current.card.append(Markdown("## Before Normalization"))
        current.card.append(
            Image.from_matplotlib(normalize_rt.before_norm_fig),
            "Before normalization",
        )
        current.card.append(Markdown("## After Normalization"))
        current.card.append(
            Image.from_matplotlib(normalize_rt.after_norm_fig),
            "After normalization",
        )

        print("Data normalized")
        self.next(self.pca, self.export_limma_contrasts)

    @card
    @step
    def pca(self):
        plot_dir = self.create_output_dir("pca")

        print("Running PCA")

        pca_output = run_pca(
            df_norm=self.counts_norm,
            metadata_maps=self.metadata_maps,
            scree_fig_name=plot_dir / "scree.png",
            pca_fig_name=plot_dir / "pca.png",
        )

        self.pca_df = pca_output.pca_df

        current.card.append(Markdown("## PCA"))
        current.card.append(Image.from_matplotlib(pca_output.pca_fig, "PCA"))

        current.card.append(Markdown("## Scree plot"))
        current.card.append(
            Markdown(
                f"Explained variance of PC 1 and 2: {pca_output.scree_output.explained_variance:.2f}%"
            )
        )
        current.card.append(
            Image.from_matplotlib(pca_output.scree_output.figure, "Scree plot")
        )

        print("PCA finished")

        self.next(self.join_pca_and_limma)

    @card
    @step
    def export_limma_contrasts(self):
        limma_input_dir = self.create_output_dir("limma_inputs")

        print("Getting contrasts")
        contrast_list = create_contrast_from_metadata(
            self.metadata_maps,
            self.preprocess_config.contrasts,
        )

        print("Exporting limma contrasts")
        self.limma_inputs = make_limma_contrasts(
            df=self.counts_norm,
            output_dir=limma_input_dir,
            contrast_list=contrast_list,
            gene_list_file=Path(self.preprocess_config.gene_input_file)
            if self.preprocess_config.gene_input_file
            else None,
        )

        print("Limma contrasts exported")

        self.next(self.run_limma, foreach="limma_inputs")

    @card
    @step
    def run_limma(self):
        print("Running limma")
        self.limma_input = self.input

        self.create_output_dir(f"{self.limma_input.contrast_name}")

        limma_output_dir = self.create_output_dir(
            f"{self.limma_input.contrast_name}/limma_outputs"
        )
        limma_sig_results_dir = self.create_output_dir(
            f"{self.limma_input.contrast_name}/limma_sig_results"
        )

        current.card.append(
            Markdown(f"## Limma for {self.limma_input.contrast_name}")
        )
        current.card.append(
            Markdown(
                f"Formula: {self.limma_input.contrast_1} - {self.limma_input.contrast_2}"
            )
        )

        self.result_path = run_limma_r(
            r_config=self.r_config,
            output_dir=limma_output_dir,
            counts=self.limma_input.counts_file,
            metadata=self.limma_input.metadata_file,
            sig_output_dir=limma_sig_results_dir,
            contrast_name=self.limma_input.contrast_name,
            contrast_1=self.limma_input.contrast_1,
            contrast_2=self.limma_input.contrast_2,
        )

        self.next(self.run_heatmap, self.run_volcano_plot, self.run_enrichment)

    @card
    @step
    def run_volcano_plot(self):
        volcano_output = self.create_output_dir(
            f"{self.limma_input.contrast_name}/volcano"
        )

        if not self.result_path:
            print("No result path found for volcano plot")
            current.card.append(
                Markdown("No DEG result found for volcano plot")
            )
            self.next(self.join_post_deg)

        print("Running volcano plot")

        current.card.append(
            Markdown(f"### Volcano plot for {self.limma_input.contrast_name}")
        )

        volcano_png = run_volcano_plot_r(
            r_config=self.r_config,
            output_dir=volcano_output,
            deg_results=self.result_path,
            experiment=self.limma_input.contrast_name,
            volcano_args=self.parameters.volcano,
        )

        if volcano_png:
            print("Plotting volcano plot: ", volcano_png)
            current.card.append(
                Image.from_pil_image(
                    PILImage.open(volcano_png), "Volcano plot"
                )
            )

        self.next(self.join_post_deg)

    @card
    @step
    def run_heatmap(self):
        heatmap_output = self.create_output_dir(
            f"{self.limma_input.contrast_name}/heatmap"
        )

        if not self.result_path:
            print("No result path found for heatmap plot")
            current.card.append(
                Markdown("No DEG result found for heatmap plot")
            )
            self.next(self.join_post_deg)

        # Define the colors
        colors = ["royalblue", "white", "red"]
        # Create the colormap
        cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

        heatmap = make_heatmap_sample(
            counts_file=self.limma_input.counts_file,
            metadata_file=self.limma_input.metadata_file,
            limma_results_file=self.result_path,
            output_dir=heatmap_output,
            make_heatmap_kwargs=MakeHeatmapOtherKwargs(
                title=self.limma_input.contrast_name,
                cmap=cmap,
                **self.parameters.heatmap.model_dump(
                    exclude={"cmap", "title"}, exclude_unset=True
                ),
            ),
        )

        current.card.append(
            Markdown(f"### Heatmap for {self.limma_input.contrast_name}")
        )
        current.card.append(
            Image.from_pil_image(PILImage.open(heatmap), "Heatmap")
        )

        self.next(self.join_post_deg)

    @card
    @step
    def run_enrichment(self):
        enrich_output = self.create_output_dir(
            f"{self.limma_input.contrast_name}/enrichment"
        )

        if not self.result_path:
            print("No result path found for enrichment analysis")
            current.card.append(
                Markdown("No DEG result found for enrichment analysis")
            )
            self.next(self.join_post_deg)

        print("Running enrichment analysis")

        enrich_results = run_enrichment_r(
            r_config=self.r_config,
            output_dir=enrich_output,
            deg_results=self.result_path,
            experiment=self.limma_input.contrast_name,
            enrichment_args=self.parameters.enrich,
        )
        print(enrich_results)

        self.kegg_results = []
        self.gene_ids = None

        current.card.append(
            Markdown(f"## Enrichment for {self.limma_input.contrast_name}")
        )
        for result in enrich_results:
            current.card.append(Markdown(f"### {result.ont} enrichment"))

            if (
                result.result_type == "enrichResult"
            ):
                current.card.append(
                    Table.from_dataframe(pd.read_csv(result.file_path).head())
                )
                if result.ont == "KEGG":
                    self.kegg_results.append(result.file_path)
            elif result.result_type == "compareResult":
                current.card.append(
                    Table.from_dataframe(pd.read_csv(result.file_path).head()
                                         )
                )
                self.kegg_results.append(result.file_path)
            elif result.result_type == "geneIds":
                print("Reading gene ids")
                self.gene_ids = pd.read_csv(result.file_path)
            elif result.result_type == "plot":
                current.card.append(
                    Image.from_pil_image(PILImage.open(result.file_path))
                )
            else:
                current.card.append(Markdown(f"Unknown result type: {result}"))

        self.next(self.fix_kegg_gene_ids)

    @step
    def fix_kegg_gene_ids(self):
        kegg_fixed_dir = self.create_output_dir(
            f"{self.limma_input.contrast_name}/kegg_fixed"
        )

        if self.gene_ids is None or len(self.gene_ids) == 0:
            print("No gene ids or kegg results found")
            current.card.append(Markdown("No gene ids or kegg results found"))
            self.next(self.join_post_deg)

        for kegg_result in self.kegg_results:
            print("Fixing KEGG ids for ", kegg_result)
            kegg_fixed = fix_kegg_ids(
                gene_ids=self.gene_ids,
                kegg_df=pd.read_csv(kegg_result),
            )
            kegg_fixed.to_csv(kegg_fixed_dir / kegg_result.name)

        print("KEGG ids fixed")

        self.next(self.join_post_deg)

    @step
    def join_post_deg(self, _inputs: Any):
        self.next(self.join_deg_analysis)

    @step
    def join_deg_analysis(self, _inputs: Any):
        self.next(self.join_pca_and_limma)

    @step
    def join_pca_and_limma(self, _inputs: Any):
        self.next(self.end)

    @step
    def end(self):
        print(get_run_output(self))
        print(get_task_output(self))
        pass


if __name__ == "__main__":
    ProteomicsAnalysis()
