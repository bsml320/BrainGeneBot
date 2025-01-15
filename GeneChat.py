import csv
import json
import re
import os
import enrichrpy.enrichr as een
from utilis import enrichment_barplot, enrichment_dotplot
import pandas as pd
# import gseapy

def generate_prompt_gene(data_entries, **keywords):
        """
        Generate a detailed prompt for deciding between displaying a protein interaction network or performing gene enrichment analysis based on given genetic data inputs related to a specific disease.

        Args:
        **keywords (dict): Keyword arguments where keys are descriptors like 'Disease', 'Genes', and other relevant factors, and values are their respective details.

        Returns:
        str: A formatted prompt string designed to output a clear decision on the analysis method along with a list of genes.
        """
        # Handle additional notes separately if provided
        additional_notes = keywords.pop('AdditionalNotes', None)
        notes_section = ""
        if additional_notes:
            notes_section = f"\n\nAdditional Notes:\n{additional_notes}"

        # Format the data entries into a JSON-like string

        if 'Disease' not in keywords:
            keywords['Disease'] = "Alzheimer’s Disease"  # Set default disease if not specified
        if 'num_gene' not in keywords:
            keywords['num_gene'] = "100"

        prompt_template = f"""
        You are a skilled genetic data analyst. Given the following {keywords['Disease']} information and related genetic factors, determine which analysis would be most beneficial: protein interaction network display or gene enrichment analysis.

        Disease and Gene Information:
        {{
            {data_entries}
        }}

        Based on the information above, decide which type of analysis to prioritize and list the top 20 genes related to the disease. Please provide your recommendation in the following JSON format:

        {{
            "recommended_function": "<1  (protein_interaction_network_display); 2  (gene_enrichment_analysis)>",
            "gene_list": [
                "Gene1",
                "Gene2",
                "Gene3"
                // Add additional genes as determined necessary
            ]
        }}

        {notes_section}
        Ensure that the JSON output is properly formatted and valid, and includes a brief rationale for the chosen analysis type.
        """
        return prompt_template

def get_input_for_string_api(response_text):
    """
    Processes an LLM's JSON-like response to extract the gene list for the STRING API based on the recommended function number.

    Args:
    response_text (str): JSON-like string response from an LLM, which includes a numeric function identifier.
    function_mapping (dict): Dictionary mapping function numbers to canonical function names, e.g., {"1": "protein interaction network display", "2": "gene enrichment analysis"}.

    Returns:
    list: List of genes for the STRING API if the response recommends a protein interaction network, otherwise an empty list.
    """
    # Extracting the recommended function number and gene list from the response
    function_match = re.search(r'"recommended_function":\s*"(\d+)\s*\(([^)]+)\)"', response_text)
    function_number = function_match.group(1) if function_match else None

    gene_list_match = re.search(r'"gene_list":\s*\[\s*([^]]+)\s*\]', response_text)
    gene_list = re.findall(r'"([^"]+)"', gene_list_match.group(1)) if gene_list_match else []

    # Check if the function number corresponds to protein interaction network display

    return function_number, gene_list  # Return the gene list if it's for protein interaction

def fetch_string_network(genes, save_dir=r'C:\Users\gqu\OneDrive - UTHealth Houston\projects\Genevic\scirpt\saved_file'):
        """
        Fetches a network image from the STRING database for a list of genes.

        Args:
        genes (list of str): List of gene identifiers.

        Returns:
        None: Saves images to files.
        """
        string_api_url = "https://version-12-0.string-db.org/api"
        output_format = "image"
        method = "network"
        ##
        ## Construct URL
        ##

        request_url = "/".join([string_api_url, output_format, method])

        ## For each gene call STRING

        for gene in genes:
            ##
            ## Set parameters
            ##

            params = {

                "identifiers": gene,  # your protein
                "species": 9606,  # species NCBI identifier
                "add_white_nodes": 15,  # add 15 white nodes to my protein
                "network_flavor": "confidence",  # show confidence links
                "caller_identity": "www.awesome_app.org"  # your app name
                # "network_type": "physical"

            }

            ##
            ## Call STRING
            ##

            response = requests.post(request_url, data=params)

            ##
            ## Save the network to file
            ##

            file_name = "%s_network.png" % gene
            print("Saving interaction network to %s" % file_name)
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            with open(os.path.join(save_dir,file_name), 'wb') as fh:
                fh.write(response.content)

if __name__ == "__main__":
    import google.generativeai as genai
    genai.configure(api_key="AIzaSyA9zIygboCR5nrj3DJWFjg1fUBQbluv0Yo")
    model = genai.GenerativeModel("gemini-1.5-flash")

    Gene_list_1= ['APOC1',
 'APOE',
 'ARID1B',
 'CCDC83',
 'CD2AP',
 'CLPTM1',
 'CR1',
 'CSMD1',
 'D2HGDH',
 'DSG4',
 'EPHA5',
 'EPHX2',
 'FAM209B',
 'GPC6',
 'LINC00943',
 'LINC01824',
 'LINC02147',
 'LOC101928219',
 'LOC102724793',
 'LOC105370500',
 'LOC105373605',
 'LOC105374773',
 'LOC105375056',
 'LOC105376942',
 'LOC105377949',
 'MAG',
 'MRPL23',
 'NME8',
 'OR52N4',
 'SLC6A18',
 'TOMM40',
 'TPM4',
 'TRAPPC6A',
 'TREM2',
 'TRIP4',
 'UBE2L3',
 'ZAN',
 'ZYX']
    Gene_list_2 = ['ADAM17',
 'BIN1',
 'CADPS2',
 'CARD11',
 'CARD11-AS1',
 'CCDC8',
 'CEACAM16',
 'CERT1',
 'COL23A1',
 'CSMD1',
 'EPHA1-AS1',
 'EXOC3L2',
 'FAM135B',
 'FAM193B-DT',
 'GPR4',
 'INPP5D',
 'LILRB1',
 'LINC00519',
 'LINC02161',
 'LOC105372334',
 'LOC105375721',
 'LOC107984118',
 'MAL2',
 'MARCHF7',
 'MCTP1',
 'MYPOP',
 'NKPD1',
 'PPP1R13L',
 'RTN2',
 'SIGLEC11',
 'SLC24A4',
 'SLC9A3',
 'TBC1D9',
 'TOMM40']
    prompt = generate_prompt_gene(
        data_entries=' '.join(Gene_list_1),
        Disease="Alzheimer's Disease",
        # AdditionalNotes="Consider the role of these genes in amyloid beta processing and their impact on neurodegeneration. Do gene enrichment analysis",
        AdditionalNotes=" Limited the number of output genes to be 10, keep the most significant ones, and show protein protein interactions",
    )
    print(prompt)
    response = model.generate_content(prompt)
    print(response.text)

    func_number, gene_list = get_input_for_string_api(response.text)
    print(func_number, gene_list)
    func_number = '1'
    import requests
    from time import sleep

    num_gene = 5
    # Example usage of the function
    # my_genes = ["YMR055C", "YFR028C", "YNL161W", "YOR373W", "YFL009W", "YBR202W"]
    if func_number == '1':
        fetch_string_network(gene_list)
    elif func_number == '2':
        df = een.get_pathway_enrichment(gene_list, gene_set_library='GO_Biological_Process_2021')
        df.to_csv("saved_file/pathway_enrichment.csv")
        pd.set_option('display.max_rows', None)  # or set a high number like 1000
        pd.set_option('display.max_columns', None)  # or set a number that covers all your columns

        # Optionally, adjust the maximum width of each column to see all the data in each cell
        pd.set_option('display.max_colwidth', None)  # or use a specific large number
        print(df.head())
        chart1 = enrichment_barplot(df, n=num_gene)
        # chart1.display()
        chart1.save('saved_file/enrichment_barplot.html')
        chart2 = enrichment_dotplot(df, n=num_gene, hue='Z-score', log=True)
        # chart2.display()
        chart2.save('saved_file/enrichment_dotplot.html')


