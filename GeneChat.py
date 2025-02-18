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
        You are a skilled genetic data analyst. Given the provided information about {keywords['Disease']} including genetic factors and associated data. 
        
        The data provided contains detailed disease and gene information:
        {{
            {data_entries}
        }}

        Based on this information, your task is to list the top 20 genes relevant to the disease for protein interaction network display. Please present your recommendation in the following JSON format:

        {{
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

    gene_list_match = re.search(r'"gene_list":\s*\[\s*([^]]+)\s*\]', response_text)
    gene_list = re.findall(r'"([^"]+)"', gene_list_match.group(1)) if gene_list_match else []

    # Check if the function number corresponds to protein interaction network display

    return gene_list  # Return the gene list if it's for protein interaction

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
    Gene_list_2 = [ 'APOE', 'TOMM40', 'APOC1', 'TREM2', 'NECTIN2', 'BCAM',
    'LOC105373605', 'SORL1', 'NCK2', 'TBC1D9', 'LOC105377557',
    'FAM193B-DT', 'SPOCK3', 'CERT1', 'CD2AP', 'MCTP1', 'CNTN2',
    'LINC02161', 'OSMR-DT', 'LINC01484', 'COL23A1', 'SLC9A3',
    'KLF7', 'FAM135B', 'TUSC3', 'LOC105375721', 'PTK2B',
    'CSMD1', 'MAL2-AS1', 'MAL2','LOC102724977', 'SLC24A4',
    'GPR68', 'LINC00519', 'PLA2R1', 'MARCHF7', 'INPP5D',
    'ADAM17', 'BIN1', 'EPHA1-AS1', 'CARD11'
    'CADPS2', 'JAZF1-AS1', 'TNS3', 'BCL3', 'RTN2',
    'ADAMTS10', 'MYPOP', 'CEACAM16-AS1', 'CEACAM16', 'PPP1R37', 'SIGLEC11']
    prompt = generate_prompt_gene(
        data_entries=' '.join(Gene_list_2),
        Disease="Alzheimer's Disease",
        AdditionalNotes="Consider the role of these genes in amyloid beta processing and their impact on neurodegeneration. Do gene enrichment analysis",
        # AdditionalNotes=" Limited the number of output genes to be 10, keep the most significant ones, and show protein protein interactions",
    )
    print(prompt)
    response = model.generate_content(prompt)
    print(response.text)
    gene_list = get_input_for_string_api(response.text)
    print(gene_list)
    func_number = '1'
    import requests
    from time import sleep

    num_gene = 100
    # Example usage of the function
    # my_genes = ["YMR055C", "YFR028C", "YNL161W", "YOR373W", "YFL009W", "YBR202W"]
    if func_number == '1':
        fetch_string_network(gene_list)
    elif func_number == '2':
        df = een.get_pathway_enrichment(Gene_list_1, gene_set_library='GO_Biological_Process_2021')
        df.to_csv("saved_file/pathway_enrichment_mc.csv")
        pd.set_option('display.max_rows', None)  # or set a high number like 1000
        pd.set_option('display.max_columns', None)  # or set a number that covers all your columns

        # Optionally, adjust the maximum width of each column to see all the data in each cell
        pd.set_option('display.max_colwidth', None)  # or use a specific large number
        print(df.head())
        chart1 = enrichment_barplot(df, n=num_gene)
        # chart1.display()
        chart1.save('saved_file/enrichment_barplot_mc.html')
        chart2 = enrichment_dotplot(df, n=num_gene, hue='Z-score', log=True)
        # chart2.display()
        chart2.save('saved_file/enrichment_dotplot_mc.html')


