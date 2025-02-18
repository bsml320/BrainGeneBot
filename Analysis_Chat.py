import csv
import json
import re
from Bio import Entrez
from tabulate import tabulate

import csv
import time
from Bio import Entrez

# Always set your email for NCBI Entrez requests.
Entrez.email = "your.email@example.com"


def generate_prompt_analy(**keywords):
    """
    Generate a detailed prompt for designing a genetic data search strategy using the given keyword inputs.

    Args:
    **keywords (dict): Keyword arguments where keys are the types of data points (e.g., disease, gene, additional notes) and values are their respective names.

    Returns:
    str: A formatted prompt string.
    """
    # Separate 'AdditionalNotes' from other keywords if present
    additional_notes = keywords.pop('AdditionalNotes', None)
    notes_section = ""
    if additional_notes:  # Include additional notes if provided
        notes_section = f"""
        Additional Notes:
        {additional_notes}
        """

    # Create a string of key-value pairs in JSON format for the rest of the keywords
    # data_entries = ", ".join(f'"{key}": "{value}"' for key, value in keywords.items())
    # if not data_entries:  # Use a generic example if no other keywords are provided
    #     data_entries = '"Example": {"Disease": "Alzheimer’s Disease", "Gene": "APOE"}'
    data_entries = ", ".join(
        f'"{key}": {json.dumps(value)}'
        for key, value in keywords.items()
    )

    # Fallback if nothing is provided:
    if not data_entries:
        # Use a generic example with multiple genes
        data_entries = '"Example": {"Disease": "Alzheimer\'s Disease", "Gene": ["APOE", "TP53"]}'

    prompt_template = f"""
    You are a skilled genetic data analyst specializing in translational medicine.

    Make use of your professional knowledge to analyze the given data from ClinVar for diseases and drug responses. 
    Summarize your findings gene by gene with detailed illustrations and supporting evidence.


    Data:
    {{
        {data_entries}
    }}

    Output Format:

    Represent your answer in JSON format, adhering to the following structure:

    {{
        
    "Gene Name": <>,
    ["Variant": <>,
    "Clinical Significance": <>,
    "Condition(s)": <>,
    "Last Evaluated": <>,
    "Submitter Information": <>]
    "Findings": <Provide a comprehensive analysis and structured summary of your findings, with a particular focus on their clinical relevance to Alzheimer’s disease.>
        
    }}

    Ensure that the JSON is properly formatted and valid.
    {notes_section}
    """
    return prompt_template

def generate_prompt_geneenrich(data_entries, **keywords):
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
        You are a skilled genetic data analyst with expertise in gene enrichment analysis. Using the provided gene enrichment results associated to {keywords['Disease']}, 
        apply your professional knowledge to analyze the data and generate a structured summary.

        Given Information:
        {{
            {data_entries}
        }}

        Please provide your recommendation in the following JSON format:

        {{
            - **Key Functional Insights:** [Describe any observed patterns, clusters, or biological significance.]  
            - **Disease Relevance:** [Explain how the identified genes contribute to {keywords['Disease']}.]  
            - **Recommended Next Steps:** [Suggest further analysis or validation techniques.]  
        }}

        {notes_section}
        Ensure that the JSON output is properly formatted and valid.
        """
    return prompt_template

def query_clinvar_for_pair(gene, rsid):
    """
    Queries ClinVar with a gene and rsID using ESearch and returns a list of ClinVar IDs.
    """
    query = f"{gene}[gene] AND {rsid}"
    with Entrez.esearch(db="clinvar", term=query, retmode="xml") as handle:
        search_results = Entrez.read(handle, validate=False)
    id_list = search_results.get("IdList", [])
    return id_list

def get_clinvar_summary(clinvar_id):
    """
    Retrieves a summary for a given ClinVar record ID using ESummary.
    """
    with Entrez.esummary(db="clinvar", id=clinvar_id, retmode="xml") as handle:
        summary_data = Entrez.read(handle, validate=False)
    return summary_data

def extract_clinvar_info(clinvar_record):
    """
    Extracts key information from the ClinVar ESummary record.
    Returns a list of cleaned records.
    """
    if not clinvar_record or "DocumentSummarySet" not in clinvar_record:
        return []

    dss = clinvar_record["DocumentSummarySet"]
    db_build = dss.get("DbBuild", "N/A")
    summaries = dss.get("DocumentSummary", [])
    if not summaries:
        return []

    extracted_records = []

    for ds in summaries:
        # Basic info
        accession = ds.get("accession", "N/A")
        accession_version = ds.get("accession_version", "N/A")
        title = ds.get("title", "N/A")
        obj_type = ds.get("obj_type", "N/A")

        # Variation details (take first entry if available)
        vs_list = ds.get("variation_set", [])
        if vs_list:
            vs = vs_list[0]
            variation_details = {
                "Measure ID": vs.get("measure_id", "N/A"),
                "Variation Name": vs.get("variation_name", "N/A"),
                "cDNA Change": vs.get("cdna_change", "N/A"),
                "Aliases": "; ".join(vs.get("aliases", [])) if vs.get("aliases") else "N/A"
            }
        else:
            variation_details = "N/A"

        # Supporting submissions
        submissions = ds.get("supporting_submissions", {})
        supporting_submissions = {
            "SCV": submissions.get("scv", "N/A") if submissions.get("scv") else "N/A",
            "RCV": submissions.get("rcv", "N/A") if submissions.get("rcv") else "N/A"
        }

        # Germline classification details
        germline_class = ds.get("germline_classification", {})
        traits = germline_class.get("trait_set", [])
        trait_names = [trait.get("trait_name", "N/A") for trait in traits if trait.get("trait_name")]
        germline_classification = {
            "Description": germline_class.get("description", "N/A"),
            "Last Evaluated": germline_class.get("last_evaluated", "N/A"),
            "Review Status": germline_class.get("review_status", "N/A"),
            "Traits": trait_names if trait_names else "N/A"
        }

        # Gene information
        genes = ds.get("genes", [])
        gene_info = []
        if genes:
            for gene_dict in genes:
                gene_info.append({
                    "Symbol": gene_dict.get("symbol", "N/A"),
                    "GeneID": gene_dict.get("GeneID", "N/A"),
                    "Strand": gene_dict.get("strand", "N/A"),
                    "Source": gene_dict.get("source", "N/A")
                })
        else:
            gene_info = "N/A"

        # Build a cleaned record
        record_data = {
            "Accession": accession,
            "Accession Version": accession_version,
            "Title": title,
            "Object Type": obj_type,
            "Variation Details": variation_details,
            "Supporting Submissions": supporting_submissions,
            "Germline Classification": germline_classification,
            "Genes": gene_info,
            "DbBuild": db_build
        }

        extracted_records.append(record_data)

    return extracted_records

def flatten_extracted_record(record, gene, rsid):
    """
    Flattens the nested ClinVar record into a single dictionary suitable for CSV output.
    """
    flat = {}
    flat["Gene"] = gene
    flat["rsID"] = rsid
    flat["Accession"] = record.get("Accession", "N/A")
    flat["Accession Version"] = record.get("Accession Version", "N/A")
    flat["Title"] = record.get("Title", "N/A")
    flat["Object Type"] = record.get("Object Type", "N/A")

    # Variation Details
    var_details = record.get("Variation Details", {})
    if isinstance(var_details, dict):
        flat["Measure ID"] = var_details.get("Measure ID", "N/A")
        flat["Variation Name"] = var_details.get("Variation Name", "N/A")
        flat["cDNA Change"] = var_details.get("cDNA Change", "N/A")
        flat["Aliases"] = var_details.get("Aliases", "N/A")
    else:
        flat["Measure ID"] = "N/A"
        flat["Variation Name"] = "N/A"
        flat["cDNA Change"] = "N/A"
        flat["Aliases"] = "N/A"

    # Supporting Submissions
    submissions = record.get("Supporting Submissions", {})
    scv = submissions.get("SCV", "N/A") if isinstance(submissions, dict) else "N/A"
    rcv = submissions.get("RCV", "N/A") if isinstance(submissions, dict) else "N/A"
    if isinstance(scv, list):
        flat["SCV"] = "; ".join(scv)
    else:
        flat["SCV"] = scv
    if isinstance(rcv, list):
        flat["RCV"] = "; ".join(rcv)
    else:
        flat["RCV"] = rcv

    # Germline Classification
    germline = record.get("Germline Classification", {})
    flat["Germline Description"] = germline.get("Description", "N/A")
    flat["Last Evaluated"] = germline.get("Last Evaluated", "N/A")
    flat["Review Status"] = germline.get("Review Status", "N/A")
    traits = germline.get("Traits", "N/A")
    if isinstance(traits, list):
        flat["Traits"] = "; ".join(traits)
    else:
        flat["Traits"] = traits

    # Genes (flatten list into a string)
    genes_info = record.get("Genes", "N/A")
    if isinstance(genes_info, list):
        gene_list = []
        for g in genes_info:
            gene_list.append(f'{g.get("Symbol", "N/A")} ({g.get("GeneID", "N/A")})')
        flat["Genes"] = "; ".join(gene_list)
    else:
        flat["Genes"] = genes_info

    flat["DbBuild"] = record.get("DbBuild", "N/A")
    return flat

def query_clinvar_and_save_csv(gene_rsid_list, output_filename):
    """
    For each (gene, rsID) pair, queries ClinVar, extracts key info,
    and writes the results to a CSV file.

    Parameters:
        gene_rsid_list (list of tuples): Each tuple is (gene, rsid)
        output_filename (str): The output CSV file name.
    """
    fieldnames = [
        "Gene", "rsID", "Accession", "Accession Version", "Title", "Object Type",
        "Measure ID", "Variation Name", "cDNA Change", "Aliases",
        "SCV", "RCV", "Germline Description", "Last Evaluated", "Review Status", "Traits",
        "Genes", "DbBuild"
    ]

    records = []
    records.append(fieldnames)
    with open(output_filename, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for gene, rsid in gene_rsid_list:
            id_list = query_clinvar_for_pair(gene, rsid)

            # If no record is found, write a row with N/A values.
            if not id_list:
              pass
            else:
                for clinvar_id in id_list:
                    summary_data = get_clinvar_summary(clinvar_id)
                    extracted_info_list = extract_clinvar_info(summary_data)
                    for extracted_record in extracted_info_list:
                        flat_record = flatten_extracted_record(extracted_record, gene, rsid)
                        writer.writerow(flat_record)
                        records.append(flat_record)
                    time.sleep(0.4)  # Pause between requests
    return records

# Example usage:

if __name__ == "__main__":
    from utilis import read_csv_as_string_table
    import google.generativeai as genai
    gene_rsid_pairs_nn = [
    ("APOE", "rs429358"),
    ("TOMM40", "rs2075650"),
    ("APOC1", "rs4420638"),
    ("APOE", "rs769449"),
    ("APOE", "rs7412"),
    ("TOMM40", "rs157582"),
    ("TREM2", "rs75932628"),
    ("TOMM40", "rs157580"),
    ("NECTIN2", "rs6859"),
    ("BCAM", "rs28399637"),
    ("LOC105373605", "rs6733839"),
    ("APOE", "rs405509"),
    ("SORL1", "rs11218343"),
    ("TOMM40", "rs1160985"),
    ("NECTIN2", "rs6857"),
    ("TOMM40", "rs405697"),
    ("NCK2", "rs143080277"),
    ("TBC1D9", "rs10004266"),
    ("LOC105377557", "rs10005919"),
    ("FAM193B-DT", "rs1001530"),
    ("SPOCK3", "rs10034594"),
    ("CERT1", "rs10038689"),
    ("CD2AP", "rs1004173"),
    ("MCTP1", "rs10051466"),
    ("CNTN2", "rs1006054"),
    ("LINC02161", "rs10060770"),
    ("OSMR-DT", "rs10061009"),
    ("LINC01484", "rs10064888"),
    ("COL23A1", "rs10065484"),
    ("CERT1", "rs10068901"),
    ("SLC9A3", "rs10076162"),
    ("KLF7", "rs10084301"),
    ("FAM135B", "rs10087496"),
    ("TUSC3", "rs10093681"),
    ("LOC105375721", "rs10093803"),
    ("PTK2B", "rs10097651"),
    ("CSMD1", "rs10100516"),
    ("MAL2", "rs10101864"),
    ("MAL2-AS1", "rs10101864"),
    ("LOC102724977", "rs10131280"),
    ("SLC24A4", "rs10135174"),
    ("GPR68", "rs10135388"),
    ("LINC00519", "rs10148013"),
    ("PLA2R1", "rs10168568"),
    ("MARCHF7", "rs10171236"),
    ("INPP5D", "rs10171658"),
    ("ADAM17", "rs10179642"),
    ("PLA2R1", "rs10188272"),
    ("BIN1", "rs10194375"),
    ("BIN1", "rs10200967"),
    ("INPP5D", "rs10202748"),
    ("EPHA1-AS1", "rs10224310"),
    ("CARD11", "rs10225022"),
    ("CARD11-AS1", "rs10225022"),
    ("EPHA1-AS1", "rs10226151"),
    ("EPHA1-AS1", "rs10228407"),
    ("CADPS2", "rs10241928"),
    ("CARD11", "rs10242979"),
    ("CARD11-AS1", "rs10242979"),
    ("JAZF1-AS1", "rs10264306"),
    ("CADPS2", "rs10268511"),
    ("TNS3", "rs10279935"),
    ("TOMM40", "rs1038026"),
    ("BCL3", "rs10401176"),
    ("RTN2", "rs10401270"),
    ("ADAMTS10", "rs10401300"),
    ("MYPOP", "rs10403030"),
    ("CEACAM16", "rs10403682"),
    ("CEACAM16-AS1", "rs10403682"),
    ("PPP1R37", "rs10405086"),
    ("SIGLEC11", "rs10405621"),
    ("PPP1R37", "rs10405859")
]


    gene_rsid_pairs_mc = [
    ("LINC00943", "rs117394726"),
    ("D2HGDH", "rs1106639"),
    ("TOMM40", "rs157583"),
    ("MAG", "rs62109573"),
    ("TOMM40", "rs34095326"),
    ("TOMM40", "rs157581"),
    ("PLCG2", "rs72824905"),
    ("TRAPPC6A", "rs150685845"),
    ("AAK1", "rs79644719"),
    ("TOMM40", "rs112019714"),
    ("LOC107986595", "rs60755019"),
    ("LOC105376942", "rs79816395"),
    ("ANKRD20A8P", "rs141342242"),
    ("APOC1", "rs4420638"),
    ("GNB1", "rs11260622"),
    ("EPHA5", "rs28660482"),
    ("ZNF618", "rs76920613"),
    ("LOC105377949", "rs72942081"),
    ("UNC13C", "rs75607373"),
    ("GPC6", "rs9516245"),
    ("NECTIN2", "rs6857"),
    ("OR52N4", "rs4910844"),
    ("EMP3", "rs4893"),
    ("TPM4", "rs60340772"),
    ("GPC6", "rs60257121"),
    ("GPC6-AS2", "rs60257121"),
    ("RAB3C", "rs13159296"),
    ("LOC105373605", "rs34779859"),
    ("SORL1", "rs12280714"),
    ("CR1", "rs10779336"),
    ("TRIM59-IFT80", "rs77393335"),
    ("CSMD1", "rs1714682"),
    ("ABI3", "rs616338"),
    ("APOE", "rs1081105"),
    ("TREM2", "rs75932628"),
    ("TOMM40", "rs77301115"),
    ("LINC01824", "rs73143087"),
    ("VWF", "rs61750595"),
    ("UBE2L3", "rs138727474"),
    ("SORL1", "rs117807585"),
    ("CCDC83", "rs116136578"),
    ("APOC1", "rs12721046"),
    ("TREM2", "rs143332484"),
    ("TRIP4", "rs74615166"),
    ("IMMP2L", "rs79936819"),
    ("AFTPH", "rs145469924"),
    ("TEX36", "rs77920798"),
    ("EPHX2", "rs142885341"),
    ("NECTIN2", "rs41289512"),
    ("ARID1B", "rs144363587"),
    ("TOMM40", "rs116881820"),
    ("CLPTM1", "rs116949436"),
    ("SORT1", "rs141749679"),
    ("MRPL23", "rs147239461"),
    ("NECTIN2", "rs12972156"),
    ("FERMT2", "rs113575650"),
    ("APOE", "rs769449"),
    ("UST", "rs114072388"),
    ("ZYX", "rs182263486"),
    ("NCK2", "rs143080277"),
    ("LOC101928219", "rs114673780"),
    ("SND1", "rs117240937"),
    ("LINC02147", "rs72778387"),
    ("NAALADL2", "rs186123476"),
    ("CD2AP", "rs34746412"),
    ("INPP5D", "rs141658619"),
    ("NME8", "rs59976014"),
    ("DSG1-AS1", "rs7241860"),
    ("DSG1-DSG4", "rs7241860"),
    ("DSG1-DSG4", "rs7241860"),
    ("MEF2C-AS1", "rs117958293"),
    ("TOMM40", "rs11556505"),
    ("ARHGAP45", "rs61242726"),
    ("FAM209B", "rs6069777"),
    ("SCAPER", "rs3812908"),
    ("ZAN", "rs72364644"),
    ("FPGS", "rs10760502"),
    ("SLC6A18", "rs7447815"),
    ("ZNF419", "rs2074071")
]


    # output_csv = "saved_file/clinvar_results_nn.csv"
    # records_nn = query_clinvar_and_save_csv(gene_rsid_pairs_nn, output_csv)
    #
    # print(f"Results saved to {output_csv}")
    #
    # output_csv = "saved_file/clinvar_results_mc.csv"
    # records_mc = query_clinvar_and_save_csv(gene_rsid_pairs_mc, output_csv)
    #
    # print(f"Results saved to {output_csv}")

    records_nn_string = read_csv_as_string_table(csv_filename="saved_file/clinvar_results_nn.csv")
    records_mc_string = read_csv_as_string_table(csv_filename="saved_file/clinvar_results_mc.csv")



    # prompt = generate_prompt_analy(Disease="Alzheimer's Disease", Gene=records_mc_string)
    # print(prompt)

    records_mc_geneenrich = read_csv_as_string_table(csv_filename="saved_file/unsup_pip/pathway_enrichment_mc.csv")
    records_nn_geneenrich = read_csv_as_string_table(csv_filename="saved_file/sup_pip/pathway_enrichment_nn.csv")

    prompt = generate_prompt_geneenrich( data_entries=records_mc_geneenrich, Disease="Alzheimer's Disease")
    print(prompt)
    # #
    genai.configure(api_key="AIzaSyA9zIygboCR5nrj3DJWFjg1fUBQbluv0Yo")
    model = genai.GenerativeModel("gemini-1.5-flash")
    response = model.generate_content(prompt)
    print(response.text)




