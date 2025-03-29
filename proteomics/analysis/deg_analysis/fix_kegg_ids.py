import pandas as pd

from proteomics.utils.base_params import get_project_root


def fix_kegg_ids(
    *,
    gene_ids: pd.DataFrame,
    kegg_df: pd.DataFrame,
    gene_id_col: str = "geneID",
    symbol_col: str = "SYMBOL",
    entrez_col: str = "ENTREZID",
) -> pd.DataFrame:
    """
    Converts the ENTREZID in the kegg_df to gene symbols using the gene_ids DataFrame.

    Args:
        gene_ids: A DataFrame with columns for gene symbols and entrez ids.
        kegg_df: A DataFrame with a column for gene ids.
        gene_id_col: The column in kegg_df that contains the gene ids.
        symbol_col: The column in gene_ids that contains the gene symbols.
        entrez_col: The column in gene_ids that contains the entrez ids.

    Returns: A DataFrame with the gene ids in kegg_df replaced with gene symbols.
    """
    id_to_gene_map = {
        str(entrez): symbol
        for symbol, entrez in gene_ids[[symbol_col, entrez_col]].itertuples(
            index=False
        )
    }

    kegg_parsed = kegg_df.copy()
    kegg_parsed[gene_id_col] = kegg_parsed[gene_id_col].astype(str).apply(
        lambda x: [id_to_gene_map.get(y) for y in x.split("/")]
    )

    return kegg_parsed


def main():
    gene_ids = pd.read_csv(
        get_project_root()
        / "data/example/_results/S2808D_HFD_vs_S2808D/enrichment/gene_ids.csv"
    )
    kegg_df = pd.read_csv(
        get_project_root()
        / "data/example/_results/S2808D_HFD_vs_S2808D/enrichment/KEGG_enrich.csv"
    )
    kegg_contrast = pd.read_csv(
        get_project_root()
        / "data/example/_results/S2808D_HFD_vs_S2808D/enrichment/KEGG_contrast.csv"
    )

    fixed_kegg_df = fix_kegg_ids(gene_ids=gene_ids, kegg_df=kegg_df)
    fixed_kegg_contrast = fix_kegg_ids(
        gene_ids=gene_ids, kegg_df=kegg_contrast
    )

    print(fixed_kegg_df.head())
    print(fixed_kegg_contrast.head())

    print(fixed_kegg_df["geneID"])


if __name__ == "__main__":
    main()
