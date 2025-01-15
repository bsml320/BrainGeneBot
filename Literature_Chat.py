import csv
import json
import re
from Bio import Entrez
from tabulate import tabulate

def generate_prompt_literature(**keywords):
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
    You are a skilled genetic data analyst with expertise in researching relevant literature.

    Design a detailed search strategy using artificial intelligence to refine and optimize genetic data retrieval from public databases. 
    Outline the approach to identify key search terms, the selection of databases, and the formulation of query parameters tailored for the E-utilities API.
    

    Data:
    {{
        {data_entries}
    }}

    Output Format:

    Represent your answer in JSON format, adhering to the following structure:

    {{
        "inference": [
            {{
                "db": "<db used for eutils>",
                "terms": "<terms generated for search>"
            }}
        ]
    }}

    Ensure that the JSON is properly formatted and valid.
    {notes_section}
    """
    return prompt_template

def parse_json_like_string(json_string):
    """
    Parses a JSON-like string to extract database names and their corresponding search terms using regex,
    handling more complex structures, cleaning terms to include only textual content, and removing duplicates.

    Args:
    json_string (str): A string that mimics JSON formatting.

    Returns:
    dict: A dictionary with database names as keys and lists of their respective search terms as values,
          cleaned and deduplicated.
    """
    data = {}
    # Regex to extract db and terms, paying attention to nested structures
    pattern = r'"db":\s*"([^"]+)"[^[]+\[([\s\S]+?)\](?=\s*[,}])'
    matches = re.finditer(pattern, json_string)

    # Process matches
    for match in matches:
        db = match.group(1)
        terms_block = match.group(2)

        # Extract terms and handle nested structures
        terms = set()  # Using a set to avoid duplicates
        term_pattern = r'"((?:[^"\\]|\\.)*)"'
        term_matches = re.finditer(term_pattern, terms_block)
        for term_match in term_matches:
            # Clean the term
            term = clean_term(term_match.group(1))
            terms.add(term)

        # Convert set back to list to preserve insertion order and store in dictionary
        data[db] = list(terms)

    return data

def clean_term(term):
    """
    Cleans a term to remove any JSON-specific encoding, extra quotes, or other non-textual elements.

    Args:
    term (str): A term possibly containing JSON-specific characters or formatting.

    Returns:
    str: A cleaned text-only term.
    """
    # Decode JSON escapes and remove any surrounding quotes
    term = term.replace('\\\\"', '"')  # Handle escaped double quotes
    term = term.replace("\\\\", "\\")  # Handle escaped backslashes
    term = re.sub(r'\\[tnr]', '', term)  # Remove escaped tabs, newlines, and returns
    return term


def extract_article_details(article_record):
    """
    Extracts details from the article record returned by PubMed.
    Handles cases where author information might be incomplete.

    Args:
    article_record (Entrez Parser Object): Parsed XML from NCBI E-utilities.

    Returns:
    dict: A dictionary containing 'PMID', 'Title', 'Authors', 'Abstract'.
    """
    try:
        article = article_record['PubmedArticle'][0]
        pmid = str(article['MedlineCitation']['PMID'])
        title = article['MedlineCitation']['Article']['ArticleTitle']

        # Handling potentially missing author details
        author_list = article['MedlineCitation']['Article'].get('AuthorList', [])
        authors = []
        for author in author_list:
            last_name = author.get('LastName', 'No LastName')
            fore_name = author.get('ForeName', '')
            authors.append(f"{last_name} {fore_name}".strip())

        # Concatenate author names into a single string
        authors_str = ', '.join(authors) if authors else "No authors listed"

        abstract = \
        article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', ['No abstract available'])[0]

        return {
            'PMID': pmid,
            'Title': title,
            'Authors': authors_str,
            'Abstract': abstract
        }
    except IndexError as e:
        print(f"Error processing article details: {e}")
        return {}


def print_article_details(articles):
    """
    Prints detailed information about each article in a readable format.

    Args:
    articles (list of dicts): List containing article details.
    """
    headers = ["PMID", "Title", "Authors", "Abstract"]
    table = []

    for article in articles:
        # Extracting each field from the article dictionary
        pmid = article['PMID']
        title = article['Title']
        authors = article['Authors']
        abstract = article['Abstract'][:500] + "..."  # Display only the first 500 characters of the abstract

        # Append to the table list
        table.append([pmid, title, authors, abstract])

    # Print the table using tabulate for better readability
    print(tabulate(table, headers=headers, tablefmt="grid"))


def search_articles(db, search_terms, api_key=None):
    """
    Search articles using NCBI E-utilities for the given search terms, optionally using an API key.

    Args:
    db (str): Database to search, e.g., 'pubmed'.
    search_terms (list): A list of search terms.
    api_key (str, optional): API key for using NCBI E-utilities.

    Returns:
    list: PMIDs of articles found.
    """
    Entrez.email = 'gang.qu@uth.tmc.edu'  # Replace with your actual email
    if api_key:
        Entrez.api_key = api_key  # Set the API key if provided

    if api_key:
        Entrez.api_key = api_key

    articles = []
    for term in search_terms:
        handle = Entrez.esearch(db=db, term=term, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        pmid_list = record['IdList']
        print(pmid_list)
        # Fetch details for each article
        for pmid in pmid_list:
            try:
                handle = Entrez.efetch(db=db, id=pmid, rettype="abstract", retmode="xml")
                article_record = Entrez.read(handle)
                handle.close()

                if not article_record or 'PubmedArticle' not in article_record:
                    print(f"Failed to retrieve data for PMID {pmid}")
                    continue

                article_details = extract_article_details(article_record)
                articles.append(article_details)
            except Exception as e:
                print(f"Error fetching details for PMID {pmid}: {e}")

    return articles

if __name__ == "__main__":
    import google.generativeai as genai

    Gene_list_1 = ['APOC1',
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
    api_key = 'f8e150b485cca1f9c1e3ac9f8149518d9c08'
    prompt = generate_prompt_literature(Disease="Alzheimer's Disease", Gene=Gene_list_2)
    print(prompt)

    # #
    genai.configure(api_key="AIzaSyA9zIygboCR5nrj3DJWFjg1fUBQbluv0Yo")
    model = genai.GenerativeModel("gemini-1.5-flash")
    response = model.generate_content(prompt)
    print(response.text)
    print(parse_json_like_string(response.text))
    data = parse_json_like_string(response.text)

     # Optional; replace with your actual NCBI API key or None

    # Parse the JSON-like string and search for articles
    for db, terms in data.items():
        if db == "PubMed":
            articles = search_articles("pubmed", terms, api_key=api_key)
            print_article_details(articles)