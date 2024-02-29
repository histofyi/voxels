from typing import Dict, List, Union

import requests
import os

from biopandas.pdb import PandasPdb

from functions.helpers import deslugify_allele


def fetch_structure_info(sql_query:str, excluded_structures:List, only_highest_resolution:bool=False) -> Dict:
    """
    This function fetches a set of rows matching a SQL query from a datasette instance (datasette.histo.fyi)

    Args:
        sql_query (str): the sql_query to fetch from datasette (should return data in a JSON format)
        excluded_structures (List): a list of PDB codes to exclude from the results

    Returns:
        Dict: the dictionary of structures which contains the following fields in a 'structures' dict (pdb_code, locus, allele, peptide, resolution and pmhc_key - a compound key of allele_slug and lowercased peptide sequence)
    """

    # build the full query URL for the datasette query
    # TODO should the datasette url be a config variable?
    query_url = 'https://datasette.histo.fyi/core.json?sql=' + sql_query

    # fetch the data from the datasette instance
    r = requests.get(query_url)

    # filter out any excluded structures from the results
    rows = [row for row in r.json()['rows'] if row[0] not in excluded_structures]

    # create the structure_info dictionary and populate the metadata and structures fields
    structure_info = {'metadata':{'excluded_structures':excluded_structures, 'sql_query':sql_query}, 'structures':{}}

    # populate the structures field with the data from the datasette query
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

    # and add the total structure count
    structure_info['metadata']['original_structure_count'] = len(structure_info['structures'])

    # if only_highest_resolution is True, filter out any structures that are not the highest resolution for their pmhc_key
    if only_highest_resolution:

        # create a new dictionary to store the truncated structure_info
        truncated_structure_info = {'metadata':{}, 'structures':{}}
        
        # create a list of used pmhc_keys
        used_pmhc_keys = []

        # brind over the metadata from the original structure_info
        truncated_structure_info['metadata'] = structure_info['metadata']
        truncated_structure_info['metadata']['only_highest_resolution'] = only_highest_resolution

        # iterate through the structures and add them to the truncated_structure_info if they are the highest resolution for their pmhc_key
        for structure in structure_info['structures']:
            if structure_info['structures'][structure]['pmhc_key'] not in used_pmhc_keys:
                truncated_structure_info['structures'][structure] = structure_info['structures'][structure]
                used_pmhc_keys.append(structure_info['structures'][structure]['pmhc_key'])
        
        # add the count of structures in the truncated_structure_info
        truncated_structure_info['metadata']['only_highest_resolution_count'] = len(truncated_structure_info['structures'])

    else:

        structure_info['metadata']['only_highest_resolution_count'] = None
        truncated_structure_info = structure_info

    # add in the value of the only_highest_resolution boolean to the metadata
    truncated_structure_info['metadata']['only_highest_resolution'] = only_highest_resolution

    # return the appropriate structure_info dictionary
    return truncated_structure_info


def download_structure(pdb_code:str, config:Dict, verbose=False) -> Union[str, bool, bool, str]:
    """
    This function downloads a PDB file from the histo.fyi website and saves it to the structures directory

    Args:
        pdb_code (str): The PDB code of the structure to be downloaded, e.g. '1hhk'.
        config (Dict): The configuration dictionary which contains the base_url for the Histo website
        verbose (bool): Whether to print verbose output
    
    Returns:
        str: The structure data as a string
        bool: Whether the structure was downloaded
        bool: Whether the structure was already cached
        str: An error message, or None if no error occurred
        
    """

    # lowercase the pdb code for both the filename and the URL
    pdb_code = pdb_code.lower()
    downloaded = False
    cached = False
    error = None

    # create the filename and the URL
    filename = f"structures/{pdb_code}_peptide.pdb"
    url = f"{config['base_url']}/{pdb_code}_1_peptide.pdb"

    if not os.path.exists('structures'):
        os.makedirs('structures')
    
    if os.path.exists(filename):
        if verbose:
            print (f"{pdb_code} already exists")
        structure_data = open(filename, 'r').readlines()
        cached = True

    else:
        if verbose:
            print (f"Downloading {pdb_code} from {url}")
        # download the file and save it to the structures directory
        r = requests.get(url)
        if r.status_code == 200:
            structure_data = r.text

            with open(filename, 'w') as filehandle:
                filehandle.write(structure_data)
            downloaded = True
        else:
            error = f"PDB code {pdb_code} has failed with status {r.status_code}"

    return structure_data, downloaded, cached, error


def load_pdb_file_to_dataframe(pdb_code:str, domain:str, config:Dict) -> Dict:
    """
    Loads a PDB file from the Histo website and returns a dataframe containing the data.

    Args:
        pdb_code (str): The PDB code of the structure to be loaded, e.g. '1hhk'.
        domain (str): The domain of the structure to be loaded, e.g. 'peptide'.

    """
    # lowercase the pdb code for both the filename and the URL
    pdb_code = pdb_code.lower()

    # create the filename and the URL
    filename = f"structures/{pdb_code}_{domain}.pdb"
    url = f"{config['base_url']}/{pdb_code}_1_{domain}.pdb"

    # if the file exists, read it into a dataframe
    if os.path.exists(filename):
        structure_df = PandasPdb().read_pdb(filename)

    # otherwise, download the file and then read it into a dataframe
    else:
    
        r = requests.get(url)
        if r.status_code == 200:
            structure_data = r.text

            with open(filename, 'w') as filehandle:
                filehandle.write(structure_data)

            structure_df = PandasPdb().read_pdb(filename)
        # if the request fails, set structure_df None and print an error message
        else:
            print (f"PDB code {pdb_code} has failed with status {r.status_code}")
            structure_df = None

    # if the structure_df is not None, return the dataframe
    if structure_df:        
        return structure_df.df
    else:
        return None