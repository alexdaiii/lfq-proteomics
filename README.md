# LFQ Proteomics Analysis

A Metaflow project to analyze label free proteomics output from Spectronaut.

This project does not perform imputation of the data
([yet](./proteomics/analysis/preprocess/impute_data.py)),
so a full matrix with imputed values is required from Perseus or
[MissForest](https://github.com/yuenshingyan/MissForest)

Data Required:

- A Excel file with the following sheets:
    - **counts matrix**
        - Can actually be named anything, but should be a full rank matrix
          where the rows are the features and the columns are the observations.
        - The first column are the the protein/gene names.
        - All subsequent columns are the observations, with the column header
          being the sample name.
    - **metadata**
        - A DeSeq2 style metadata table with the following columns:
            - **sample**
                - The sample name. Must match the column headers in the counts
                  matrix.
                  Copy and paste transpose the count matrix header if necessary.
                - This should be the first column.
            - **condition**
                - The condition of the sample. Can actually have any name,
                  but by default it will look for "condition".

It performs the following steps:

1. **Median Center Normalization and Log2 Transformation:**
    - Normalizes the data to the median and applies log2 transformation.

2. **Principal Component Analysis (PCA):**
    - Performs PCA to visualize the variability in the data.

3. **Differential Expression Analysis Using limma:**
    - Creates a design matrix using the formula
      `CONDITION_vs_CONTROL = CONDITION - CONTROL`.
    - Uses a linear model for full biological replicates:
        - Each sample is a completely independent biological replicate
          (e.g., samples from different patients, mice, etc.).
        - No batch/technical replicates column is included.
        - Model formula: `~0 + Condition`.

4. **Volcano Plot of the Results:**
    - Generates a volcano plot to visualize differentially expressed proteins.

5. **Heatmap of All Differentially Expressed Genes:**
    - Creates a heatmap to visualize the expression levels of
      all differentially expressed genes.

6. **Enrichment Analysis of GO and KEGG Pathways:**
    - Performs enrichment analysis to identify significantly
      enriched Gene Ontology (GO) terms and KEGG pathways.

## How to run the project

### Docker

1. Create a config file. See [config-example.yaml](./config-example.yaml) for an
   example.
2. Run the metaflow script.
    - You will need to mount an output directory
      to `/home/mambauser/code/.metaflow` to save the results.
    - You will also need to mount the config file and input data
      to somewhere inside the container.
    - The --max-workers flag is optional and can limit the number of
      workers metaflow uses if you are running out of memory.
    - Since the workdir in the contaienr is `/home/mambauser` the metaflow
      local storage is at `/home/mambauser/.metaflow`.

```shell
docker run \
-v ./data:/home/mambauser/data \
-v ./config-example.yaml:/home/mambauser/config.yaml \
-v ./output:/home/mambauser/.metaflow \
--pull missing alexdaiii/lfq-proteomics-limma \
run \
--max-workers 2 \
--params_file /home/mambauser/config.yaml
```

### Conda

1. Clone or download the repository
2. Add the project to the PYTHONPATH

```shell
export PYTHONPATH=$(pwd):$PYTHONPATH
```

3. Create a config file. See [config-example.yaml](./config-example.yaml) for an
   example.
4. Run the metaflow script

```shell
python proteomics/run_analysis.py run \
--params_file ./config-example.yaml
```

## How build the project

### Conda

This is only for x86 since some packages are not available for Apple Silicon.

```shell
conda env create -f environment.yaml
```

### Docker

An example of how to build the image for ARM64.

```shell
docker buildx build \
--platform linux/arm64 \
-t alexdaiii/lfq-proteomics-limma:0.0.1 \
-t alexdaiii/lfq-proteomics-limma:latest \
--load \
.
```

## TODO:

- Add imputation of missing values using MissForest.
- Add step to compare standard t-test and limma results.
- Github action to build and push the docker image. (currently only an ARM image
  on [Docker Hub](https://hub.docker.com/r/alexdaiii/lfq-proteomics-limma))