from typing import Dict

import requests
import json
import os


from constants import excluded_structures, base_url


def deslugify_allele(allele_slug:str) -> str:
    """
    This function turns an allele slug (for HLA) back into an IMGT allele number

    Args:
        allele_slug (str): the slugified version of the allele numer e.g. hla_a_02_01

    Returns:
        str: the allele number corresponding to the allele_slug e.g. HLA-A*02:01
    """
    allele_components = allele_slug.split('_')
    allele_number = f"{allele_components[0]}-{allele_components[1]}*{allele_components[2]}:{allele_components[3]}"
    return allele_number.upper()


def fetch_structure_info(sql_query:str) -> Dict:
    """
    This function fetches a set of rows matching a SQL query from a datasette instance (datasette.histo.fyi)

    Args:
        sql_query (str): the sql_query to fetch from datasette (should return data in a JSON format)

    Returns:
        Dict: the dictionary of structures which contains the following fields in a 'structures' dict (pdb_code, locus, allele, peptide, resolution and pmhc_key - a compound key of allele_slug and lowercased peptide sequence)
    """
    query_url = 'https://datasette.histo.fyi/core.json?sql=' + sql_query
    r = requests.get(query_url)

    rows = [row for row in r.json()['rows'] if row[0] not in excluded_structures]

    structure_info = {'metadata':{}, 'structures':{}}

    for row in rows:
        structure_info['structures'][row[0]] = {
            'pdb_code': row[0],
            'locus': row[1],
            'allele': deslugify_allele(row[2]),
            'allele_slug': row[2],
            'peptide': row[3],
            'pmhc_key': f"{row[2]}_{row[3].lower()}",
            'resolution': row[4]
        }

    return structure_info


# Define the query variables, conditions and order_by clause for the SQL query
# TODO this may be overkill for a simple script and may make it harder to maintain/reuse for other queries
query_variables = ['pdb_code', 'locus', 'allele_slug', 'peptide_sequence', 'resolution']
conditions = {
    'peptide_length':'9', 
    'peptide_features':'correct_sequence_and_length', 
    'complex_type':'class_i_with_peptide'
}
order_by = 'resolution asc'
loci = ['hla-a', 'hla-b', 'hla-c']

sql_query = "select " + ", ".join(query_variables) + " from core where "

for condition in conditions:
    sql_query += f"{condition}='{conditions[condition]}' and "

or_query = " or ".join([f"locus='{locus}'" for locus in loci])
sql_query += '(' + or_query + ')'

sql_query += f" order by {order_by}"

# run the query against datasette and fetch the structure information
structure_info = fetch_structure_info(sql_query)

# add the metadata to the structure_info dictionary
structure_info['metadata']['sql_query'] = sql_query
structure_info['metadata']['query_variables'] = query_variables
structure_info['metadata']['conditions'] = conditions
structure_info['metadata']['order_by'] = order_by
structure_info['metadata']['loci'] = loci


with open('output/structure_information/all.json', 'w') as filehandle:
    json.dump(structure_info, filehandle, indent=4)

if not os.path.exists('structures'):
    os.mkdir('structures')

i = 0

downloaded = []
cached = []
for pdb_code in structure_info['structures']:
    structure_filepath = f"structures/{pdb_code}_peptide.pdb"

    if not os.path.exists(structure_filepath):

        url = f"{base_url}/{pdb_code}_1_peptide.pdb"

        r = requests.get(url)
        if r.status_code == 200:
            structure_data = r.text

            with open(structure_filepath, 'w') as filehandle:
                filehandle.write(structure_data)

            allele = structure_info[pdb_code]['allele']
            peptide = structure_info[pdb_code]['peptide']
            pmhc_key = structure_info[pdb_code]['pmhc_key']
            resolution = structure_info[pdb_code]['resolution']
            print(f"Downloaded {pdb_code.upper()} - {allele} binding {peptide} at {resolution}Å resolution")

            downloaded.append(pdb_code)
    else:
        cached.append(pdb_code)
    i += 1

print(f"{i} structures ready for analysis. {len(downloaded)} downloaded, {len(cached)} previously downloaded")