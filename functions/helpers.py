from typing import List, Dict

import json
import os

def load_config() -> Dict:
    """
    Loads a JSON configuration file and returns the data as a dictionary.


    Returns:
        dict: The data from the JSON configuration file.
    """
    config_filename = 'config.json'
    # open the config json file and load the data into a dictionary
    if os.path.exists(config_filename):
        try:
            return json.load(open(config_filename, 'r'))
        except json.JSONDecodeError:
            print (f"Config file {config_filename} is not a valid JSON file")
            return None
    else:
        print (f"Config file {config_filename} not found")
        return None



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