import csv
import json


def generate_prompt_PGS(disease_category, top_n_variants, rank_result=None):
    # Define the prompt components
    content = (f"You are an expert genetic data analyst specializing for {disease_category}.")

    chain_of_thought = (
        "Please perform a detailed, s step-by-step analysis for the following top-ranked variants and their associated genes. "
        "Ensure that each step of your reasoning process is thoroughly explained.")

    task_instructions = (f"**Task:**\n\n"
                         f"- Integrate your knowledge with the ranks of each genetic variant obtained from given data.\n"
                         f"- Identify the top {top_n_variants} unique genes from the ranking results, ensuring no repetitions in identified genes, in accordance with your professional expertise.\n"
                         f"- Provide a detailed summary explaining how the final ranks were determined."
                         f"- Analysis your results and provide the common knowledge to verify your results, listing at least 2 evidences for each identified gene.")
# list common knowledge AD drive genes
    output_format = ("**Output Format:**\n\n"
                     "Represent your answer in JSON format, adhering to the following structure:\n\n"
                     "{\n"
                     '  "inference": [\n'
                     '    {\n'
                     '      "Variant_rsID": "<variant_id>",\n'
                     '      "gene": "<gene_name>",\n'
                     '    },\n'
                     '    ...\n'
                     '  ],\n'
                     '  "summary": "<Detailed explanation of the function of the top ranked gene related to the disease.>"\n'
                     '  "Verification": "<Detailed explanation of the verification and support evidence of the results.>"\n'
                     "}\n\n"
                     "Ensure that the JSON is properly formatted and valid.")

    # Read data from CSV file
    data_section = "**Data:**\n\nVariant and associated gene rankings from the aggregation methods:\n\n[\n"
    variants = rank_result

    print(variants)
    # Convert variants list to JSON string with proper indentation
    data_section += variants
    data_section += "\n]\n\n*Note:* Lower numerical ranks indicate higher importance (Rank 1 is the highest)."

    # Combine all parts to form the final prompt
    prompt = f"{content}\n\n{chain_of_thought}\n\n{task_instructions}\n\n{output_format}\n\n{data_section}"

    return prompt


# Example usage:
if __name__ == "__main__":
    import google.generativeai as genai
    disease_category = "Alzheimer's disease"
    # csv_filename = "variant_ranks.csv"  # Replace with your CSV file name
    top_n_variants = 10  # Specify the number of top variants you want
    rank_result = """"
    rs117394726  & LINC00943       & rs117394726  & LINC00944        \\
rs1106639    & D2HGDH          & rs157583     & TOMM40           \\
rs62109573   & MAG             & rs34095326   & TOMM40           \\
rs157581     & TOMM40         & rs72824905   & PLCG2            \\
rs150685845  & TRAPPC6A        & rs79644719   & AAK1             \\
rs112019714  & TOMM40          & rs60755019   & LOC107986595     \\
rs79816395   & LOC105376942    & rs141342242  & ANKRD20A8P       \\
rs4420638    & APOC1           & rs11260622   & GNB1             \\
rs28660482   & EPHA5           & rs76920613   & ZNF618           \\
rs72942081   & LOC105377949    & rs75607373   & UNC13C           \\
rs9516245    & GPC6            & rs6857       & NECTIN2          \\
rs4910844    & OR52N4          & rs4893       & EMP3             \\
rs60340772   & TPM4            & rs60257121   & GPC6-AS2         \\
rs60257121   & GPC6            & rs13159296   & RAB3C            \\
rs34779859   & LOC105373605    & rs12280714   & SORL1            \\
rs10779336   & CR1             & rs77393335   & TRIM59-IFT80     \\
rs1714682    & CSMD1           & rs616338     & ABI3             \\
rs1081105    & APOE            & rs75932628   & LOC105375056     \\
rs75932628   & TREM2           & rs77301115   & TOMM40           \\
rs73143087   & LINC01824       & rs61750595   & VWF              \\
rs138727474  & UBE2L3          & rs117807585  & SORL1            \\
rs116136578  & CCDC83          & rs12721046   & APOC1            \\
rs143332484  & LOC105375056    & rs143332484  & TREM2            \\
rs74615166   & TRIP4           & rs79936819   & IMMP2L           \\
rs145469924  & LOC105374773    & rs145469924  & AFTPH            \\
rs77920798   & LOC102724793    & rs77920798   & TEX36            \\
rs142885341  & EPHX2           & rs41289512   & NECTIN2          \\
rs144363587  & ARID1B          & rs116881820  & TOMM40           \\
rs116949436  & CLPTM1          & rs141749679  & SORT1            \\
rs147239461  & MRPL23          & rs12972156   & NECTIN2          \\
rs113575650  & LOC105370500    & rs113575650  & FERMT2           \\
rs769449     & APOE            & rs114072388  & UST              \\
rs182263486  & ZYX             & rs143080277  & NCK2             \\
rs114673780  & LOC101928219    & rs117240937  & SND1             \\
rs72778387   & LINC02147       & rs186123476  & NAALADL2         \\
rs34746412   & CD2AP           & rs141658619  & INPP5D           \\
rs59976014   & NME8            & rs7241860    & DSG1-AS1         \\
rs7241860    & DSG4            & rs117958293  & MEF2C-AS1        \\
rs11556505   & TOMM40          & rs61242726   & ARHGAP45         \\
rs6069777    & FAM209B         & rs3812908    & SCAPER           \\
rs72364644   & ZAN             & rs10760502   & FPGS             \\
rs7447815    & SLC6A18         & rs2074071    & ZNF419           \\"""

    rank_result_1 = """
            rs429358 & APOE & rs2075650 & TOMM40 \\
            rs4420638 & APOC1 & rs769449 & APOE \\
            rs7412 & APOE & rs157582 & TOMM40 \\
            rs75932628 & TREM2 & rs157580 & TOMM40 \\
            rs6859 & NECTIN2 & rs28399637 & BCAM \\
            rs6733839 & LOC105373605 & rs405509 & APOE \\
            rs11218343 & SORL1 & rs1160985 & TOMM40 \\
            rs6857 & NECTIN2 & rs405697 & TOMM40 \\
            rs143080277 & NCK2 & rs10004266 & TBC1D9 \\
            rs10005919 & LOC105377557 & rs1001530 & FAM193B-DT \\
            rs10034594 & SPOCK3 & rs10038689 & CERT1 \\
            rs1004173 & CD2AP & rs10051466 & MCTP1 \\
            rs1006054 & CNTN2 & rs10060770 & LINC02161 \\
            rs10061009 & OSMR-DT & rs10064888 & LINC01484 \\
            rs10065484 & COL23A1 & rs10068901 & CERT1 \\
            rs10076162 & SLC9A3 & rs10084301 & KLF7 \\
            rs10087496 & FAM135B & rs10093681 & TUSC3 \\
            rs10093803 & LOC105375721 & rs10097651 & PTK2B \\
            rs10100516 & CSMD1 & rs10101864 & MAL2(-AS1) \\
            rs10131280 & LOC102724977 & rs10135174 & SLC24A4 \\
            rs10135388 & GPR68 & rs10148013 & LINC00519 \\
            rs10168568 & PLA2R1 & rs10171236 & MARCHF7 \\
            rs10171658 & INPP5D & rs10179642 & ADAM17 \\
            rs10188272 & PLA2R1 & rs10194375 & BIN1 \\
            rs10200967 & BIN1 & rs10202748 & INPP5D \\
            rs10224310 & EPHA1-AS1 & rs10225022 & CARD11(-AS1) \\
            rs10226151 & EPHA1-AS1 & rs10228407 & EPHA1-AS1 \\
            rs10241928 & CADPS2 & rs10242979 & CARD11(-AS1) \\
            rs10264306 & JAZF1-AS1 & rs10268511 & CADPS2 \\
            rs10279935 & TNS3 & rs1038026 & TOMM40 \\
            rs10401176 & BCL3 & rs10401270 & RTN2 \\
            rs10401300 & ADAMTS10 & rs10403030 & MYPOP \\
            rs10403682 & CEACAM16(-AS1) & rs10405086 & PPP1R37 \\
            rs10405621 & SIGLEC11 & rs10405859 & PPP1R37 \\
    """
    prompt = generate_prompt_PGS(disease_category, top_n_variants, rank_result_1)
    print(prompt)

    # #
    genai.configure(api_key="AIzaSyA9zIygboCR5nrj3DJWFjg1fUBQbluv0Yo")
    model = genai.GenerativeModel("gemini-1.5-flash")
    response = model.generate_content(prompt)
    print(response.text)