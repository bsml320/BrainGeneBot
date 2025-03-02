import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import random
import torch
from Bio import Entrez
def enrichment_barplot(enrichr_results, n=20, scheme='blues'):
    """
    Plots enrichment results in barplot form

    Parameters
    ----------
    enrichr_results: pd.DataFrame
      - result dataframe from enrichrpy.enrichr.get_pathway_enrichment
    n: int
      - plot top N pathways, default=20
    scheme: str
      - altair color scheme to use. schemes listed here https://vega.github.io/vega/docs/schemes/
    """
    source = enrichr_results.copy()
    source['Num hits'] = [len(ls) for ls in source['Overlapping genes']]
    source['-log10(FDR)'] = -np.log10(source['Adjusted p-value'])
    source['Pathway'] = source['Term name'].to_list()

    if n is not None:
        source = source.sort_values('Adjusted p-value').iloc[:n]

    c = alt.Chart(source).mark_bar().encode(
        x=alt.X('-log10(FDR)'),
        y=alt.Y('Pathway', sort={"encoding": "x", "order": "descending"}),
        color=alt.Color('Num hits', scale=alt.Scale(scheme=scheme, domainMin=0))
    )
    xrule = (
        alt.Chart()
        .mark_rule(strokeDash=[8, 6], color="red", strokeWidth=2)
        .encode(x=alt.datum(-np.log10(.05)))
    )

    return c + xrule


def enrichment_dotplot(enrichr_results, n=20, hue='Z-score', scheme='viridis', log=True):
    """
    Plots enrichment results in dotplot form

    Parameters
    ----------
    enrichr_results: pd.DataFrame
      - result dataframe from enrichrpy.enrichr.get_pathway_enrichment
    n: int
      - plot top N pathways, default=20
    hue: str
      - variable to color the dotplot by, default='Combined score'
    scheme: str
      - altair color scheme to use. schemes listed here https://vega.github.io/vega/docs/schemes/
    """
    source = enrichr_results.copy()
    source['Num hits'] = [len(ls) for ls in source['Overlapping genes']]
    source['-log10(FDR)'] = -np.log10(source['Adjusted p-value'])
    source['Pathway'] = source['Term name'].to_list()
    source[f'log({hue})'] = np.log(source[hue])

    if n is not None:
        source = source.sort_values('Adjusted p-value').iloc[:n]

    c = alt.Chart(source).mark_circle().encode(
        x=alt.X('-log10(FDR):Q'),
        y=alt.Y('Pathway', sort={"encoding": "x", "order": "descending"}),
        size=alt.Size('Num hits'),
        color=alt.Color(hue if not log else f'log({hue})', scale=alt.Scale(scheme=scheme, domainMin=0))
    )
    xrule = (
        alt.Chart()
        .mark_rule(strokeDash=[8, 6], color="red", strokeWidth=2)
        .encode(x=alt.datum(-np.log10(.05)))
    )

    return (c + xrule).configure_axis(grid=True)


def seed_it(seed):
    """
    set random seed for reproducibility
    :param seed: seed number
    :type seed: int
    :return:
    :rtype:
    """
    random.seed(seed)
    os.environ["PYTHONSEED"] = str(seed)
    np.random.seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.enabled = True
    torch.manual_seed(seed)


def rsid2gene(rsids):
    names = []
    for snp_id in rsids:
        numeric_id = snp_id.replace('rs', '')
        record = Entrez.read(Entrez.elink(dbfrom="snp",
                                          id=numeric_id,
                                          db="gene"))
        # Check if LinkSetDb and links are present
        if 'LinkSetDb' in record[0] and len(record[0]['LinkSetDb']) > 0:
            links = record[0]['LinkSetDb'][0]['Link']
            for gene_link in links:
                gene_id = gene_link['Id']
                handle = Entrez.esummary(db="gene", id=gene_id)
                uid_record = Entrez.read(handle)
                handle.close()

                uid_summary = uid_record["DocumentSummarySet"]['DocumentSummary'][0]
                names.append(uid_summary['Name'])
        else:
            print(f"No gene links found for {snp_id}")
    return names



def read_csv_as_string_table(csv_filename, max_col_width=5000, delimiter=","):
    """
    Reads a CSV file and converts it into a structured string table for direct API input.
    Handles special cases like missing values, encoding issues, and long text truncation.

    Parameters:
        csv_filename (str): Path to the CSV file.
        max_col_width (int): Max width of each column for display (truncation).
        delimiter (str): Delimiter used in the CSV file.

    Returns:
        str: The CSV content formatted as a structured string table.
    """
    try:
        # Read CSV into a DataFrame with automatic encoding detection
        df = pd.read_csv(csv_filename, delimiter=delimiter, encoding="utf-8", dtype=str, na_filter=False)

        # Truncate long fields for readability
        df = df.applymap(lambda x: (x[:max_col_width] + "...") if len(x) > max_col_width else x)

        # Ensure all missing values are replaced with a placeholder
        df.fillna("N/A", inplace=True)

        # Convert DataFrame to a formatted string table (without index)
        table_str = df.to_string(index=False)

        return table_str

    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return None  # Return None if the file cannot be read
