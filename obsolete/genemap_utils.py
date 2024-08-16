import os

# Define the global variable for the data path
#DATA_PATH = '/home/nbatada/data/geneinfo/genemaps/'
DATA_PATH = os.path.join( os.path.expanduser("~"), 'data', 'geneinfo','maps') 

# Example data handle information
data_handles_info = {
    'ensgeneid_to_genesymbol': {
        'filename': 'gene_aliases/genemap_geneid_to_genesymbol.tsv',
        'info': 'parsed from gencode (v46)',
        'species': 'human',
        'source': 'https://www.gencodegenes.org/human/release_46.html'
    },
    'genesymbol_to_geneid': {
        'filename': os.path(DATA_PATH, 'hgnc','genesymbols2geneid_hgnc.tsv')
        'info': 'parsed from hgnc'
        'species': 'human',
        'source': 'https://www.genenames.org/download/archive/ go down to Current tab separated hgnc_complete_set file'
    },

}

def getAvailableData():
    """
    Returns a dictionary of data_handles and their info from the data_handles_info dictionary.
    """
    return data_handles_info.keys()

def getFile(data_handle):
    """
    Returns the full file path for the given data_handle.
    
    Args:
    data_handle (str): The key of the data handle to look up.
    
    Returns:
    str: The full file path for the file associated with the data_handle.
    """
    if data_handle in data_handles_info:
        filename = data_handles_info[data_handle]['filename']
        return os.path.join(DATA_PATH, filename)
    else:
        raise ValueError(f"Data handle '{data_handle}' not found.")

# Example usage
if __name__ == "__main__":
    print("Available gene maps:")
    available_data = getAvailableData()
    for key, value in available_data.items():
        print(f"{key}: {value['info']}")

    try:
        path = getFilePath('ensgeneid_to_genesymbol')
        print(f"File path for 'ensgeneid_to_genesymbol': {path}")
    except ValueError as e:
        print(e)
        
