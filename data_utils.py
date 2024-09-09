import re
import pandas as pd
import importlib.util
from pathlib import Path
import shutil
from datetime import datetime
import os
HOME=os.path.expanduser("~")
import sys
sys.path.insert(0,f'{HOME}/lib')
import utils
import logging
database = {}
DATA_DIR = os.path.join(HOME, 'data', 'rnaseq')
DATABASE_FILE_PATH = os.path.join(os.path.dirname(__file__), 'data_utils_database.txt')

def __load_database(file_path=DATABASE_FILE_PATH):
    global database
    database = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            key, value = line.split('=', 1)
            try:
                data_handle, sub_key = key.split('.', 1)
            except ValueError:
                print(f"Skipping invalid line: {line}")
                continue

            data_handle = data_handle.upper()

            if data_handle not in database:
                database[data_handle] = {}
            database[data_handle][sub_key] = value

__load_database(DATABASE_FILE_PATH)

DATA_ORDER = {
    'experiment.species': 1,
    'experiment.tissue': 2,
    'experiment.disease': 3,
    'experiment.disease_list': 4,
    'experiment.celltype': 5,
    'experiment.celltype_definition': 6,
    'experiment.design': 7,
    'experiment.nsamples': 8,
    'assay.type': 9,
    'assay.platform': 10,
    'assay.ncells': 11,
    'files.raw.counts': 12,
    'files.raw.meta': 13,
    'files.raw.source': 14,
    'files.processed.counts': 15,
    'files.processed.normalized': 16,
    'files.processed.meta': 17,
    'files.processed.rds': 18,
    'files.processed.source': 19,
    'files.processed.qc': 20,
    'files.processed.deg': 21,
    'files.processed.notes':22,
    'files.reprocessed.counts': 23,
    'files.reprocessed.normalized': 24,
    'files.reprocessed.meta': 25,
    'files.reprocessed.scrnaseq_obj': 26,
    'files.reprocessed.source': 27,
    'files.reprocessed.qc': 28,
    'files.reprocessed.deg': 29,
    'files.reprocessed.notes':30,
    'manuscript.info': 31,
    'manuscript.doi': 32,
    'manuscript.url_supplementary_data': 33,
    'manuscript.notes': 34,
    'manuscript.pdf':35,
    'files.raw_dir': 36,
    'files.processed_dir': 37,
    'files.reprocessed_dir': 38,
    'colname_celltype':39,
    'colname_celltype_value_treg':40, # name of the column value in the cell type
    'colname_tissue':41, # name of the column where tissue information is given
    'colname_treg_subtype':42,
    'colname_group':43,
    'files.analysis_ready':44,
}

def __update_database(file_path=DATABASE_FILE_PATH):
    if os.path.isfile(file_path):
        try:
            with open(file_path, 'r+'):
                pass
        except IOError:
            input(f"The file {file_path} is currently open. Please close it and press Enter to continue.")

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file_path = f"{file_path}.{timestamp}.bak"
    shutil.copy(file_path, backup_file_path)

    with open(file_path, 'w') as file:
        for data_handle, sub_dict in database.items():
            data_handle = data_handle.upper()
            for sub_key, value in sub_dict.items():
                line = f"{data_handle}.{sub_key}={value}\n"
                file.write(line)

def getInfo(data_handle):
    data_handle = data_handle.upper()
    info = []

    data_full = database.get(data_handle, {})
    data  ={}
    for key, value in data_full.items():
        data[key] = value
    return(data)


def deleteData(data_handle):
    data_handle = data_handle.upper()
    if data_handle not in database:
        print(f"Data handle '{data_handle}' not found.")
        return

    non_empty_count = sum(1 for v in database[data_handle].values() if v.strip())

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file_path = f"{DATABASE_FILE_PATH}.{timestamp}.bak"
    shutil.copy(DATABASE_FILE_PATH, backup_file_path)

    del database[data_handle]
    __update_database()
    
    # print(f"Data handle '{data_handle}' has been deleted and the database has been updated.")

def listAvailableData():
    """Prints out the keys in the database dictionary."""
    return(list(database.keys()))


def newDataTemplate(data_handle):
    """Print a template for the specified data handle using DATA_ORDER."""
    data_handle = data_handle.upper()
    def print_nested_dict(prefix, d):
        """Recursively print nested dictionary keys with prefix."""
        for k, v in d.items():
            if isinstance(v, dict):
                print_nested_dict(f"{prefix}.{k}", v)
            else:
                print(f"{prefix}.{k}=")

    for key in sorted(DATA_ORDER, key=DATA_ORDER.get):
        parts = key.split('.')
        if len(parts) > 1:
            print_nested_dict(data_handle, {parts[0]: {k: '' for k in parts[1:]}})
        else:
            print(f"{data_handle}.{key}=")

def help():
    import inspect

    module = inspect.getmodule(inspect.currentframe())

    functions = [name for name, obj in inspect.getmembers(module, inspect.isfunction)]

    print("Exported Functions:")
    i = 1
    for func_name in functions:
        if func_name.startswith('help') or func_name.startswith('__'): continue
        func = getattr(module, func_name)
        signature = inspect.signature(func)
        docstring = inspect.getdoc(func)

        args = [f"{param.name}={param.default}" if param.default != inspect.Parameter.empty else param.name for param in signature.parameters.values()]
        args_str = ", ".join(args)
        print(f"{i})  - {func_name}({args_str})")
        i += 1
        if docstring:
            print(f"          {docstring}")
            print("\n")

            
def updateInfo(data_handle=None, help=False, **kwargs):
    """Update the information for a given data handle with key=value pairs. provide it as **{key:value}. Keys are defined in DATA_ORDER and should be full (i.e. files.processed_dir not processed_dir)"""
    
    if help:
        print("Usage: updateInfo(data_handle, key1=value1, key2=value2, ...)")
        print("Keys available for update:", '\n'.join(DATA_ORDER.keys()))
        return

    if not data_handle:
        print("Error: 'data_handle' must be provided unless 'help=True' is specified.")
        return

    data_handle = data_handle.upper()
    
    # Initialize if data_handle doesn't exist in the database
    if data_handle not in database:
        database[data_handle] = {key: '' for key in DATA_ORDER}
    
    # Check each provided key and warn if it's not valid
    valid_keys = DATA_ORDER.keys()
    for key, value in kwargs.items():
        if key not in valid_keys:
            # Suggest similar valid keys
            closest_match = [k for k in valid_keys if key in k]
            suggestion = f" Instead, specify '{closest_match[0]}'" if closest_match else ""
            print(f"Warning: Key '{key}' not found in DATA_ORDER.{suggestion}")
        else:
            # Update the database if the key is valid
            database[data_handle][key] = value

    # Persist the database updates
    __update_database()


def build_key_mapping(data_order):
    """
    Build a dictionary to map partial keynames to full keynames from DATA_ORDER.
    """
    key_mapping = {}
    for key in data_order:
        if key.startswith('files.'):
            continue
        
        parts = key.split('.')
        for i in range(len(parts)):
            partial_key = '.'.join(parts[i:])
            if partial_key not in key_mapping:
                key_mapping[partial_key] = []
            key_mapping[partial_key].append(key)
    
    # Remove ambiguous mappings by checking if a partial key maps to multiple full keys
    for partial_key in key_mapping:
        if len(set(key_mapping[partial_key])) > 1:
            print(f"Warning: Partial key '{partial_key}' maps to multiple full keys: {set(key_mapping[partial_key])}")
    
    return key_mapping

def searchDatabase(query, restrict_to_key=None, restrict_to_keyvalue=None):
    """
    Search for data_handles in the database based on the query and optional restrictions.
    """
    if query.lower() == 'help':
        help_searchDatabase()
        return {}

    result = {}
    query = query.lower()
    query_upper = query.upper()  # Convert query to uppercase for data_handle comparison

    # Build key mapping dictionary
    key_mapping = build_key_mapping(DATA_ORDER)

    # Parse restrict_to_keyvalue
    keyname = None
    value = None
    if restrict_to_keyvalue:
        keyname, value = (restrict_to_keyvalue.split('=', 1) if '=' in restrict_to_keyvalue else (None, None))
        keyname = keyname.strip().lower() if keyname else None
        value = value.strip().lower() if value else None

    # Define keys to ignore
    ignore_keys = [k for k in DATA_ORDER if k.startswith('files.')]

    for data_handle, data in database.items():
        data_handle = data_handle.upper()  # Ensure data_handle is uppercase
        key_matches = []

        # Check if query matches the data_handle
        if query_upper in data_handle:
            result.setdefault(data_handle, [])

        if restrict_to_key:
            restrict_to_key = restrict_to_key.strip().lower()
            restrict_to_key_parts = restrict_to_key.split('.')
            full_keyname = '.'.join(restrict_to_key_parts)
            partial_keyname = restrict_to_key_parts[-1]

            if full_keyname in ignore_keys:
                continue

            # Resolve partial keynames
            resolved_keys = key_mapping.get(partial_keyname, [])
            if full_keyname in data:
                if query in str(data[full_keyname]).lower():
                    key_matches.append(full_keyname)
            elif resolved_keys:
                for resolved_key in resolved_keys:
                    if query in str(data.get(resolved_key, '')).lower():
                        key_matches.append(resolved_key)

            # Handle restrict_to_keyvalue
            if key_matches:
                if keyname:
                    keyname_matches = any(k.endswith(keyname) for k in data)
                    if keyname_matches:
                        if value and (any(str(data.get(k)).lower() == value for k in data if keyname in k.lower())):
                            result.setdefault(data_handle, []).extend(key_matches)
                else:
                    result.setdefault(data_handle, []).extend(key_matches)

        else:
            matched_keys = [k for k in data if k not in ignore_keys and query in str(data[k]).lower()]
            if matched_keys:
                if keyname:
                    if keyname in data and (str(data.get(keyname)).lower() == value):
                        result[data_handle] = matched_keys
                else:
                    result[data_handle] = matched_keys

    if not result:
        print(f"No data record found for the query '{query}' and/or restrict_to_key '{restrict_to_key}' and/or restrict_to_keyvalue '{restrict_to_keyvalue}'")

    return result

def help_searchDatabase():
    """Print usage information for the searchDatabase function."""
    print("Usage:")
    print("    searchDatabase(query, restrict_to_key=None, restrict_to_keyvalue=None)")
    print()
    print("Parameters:")
    print("    query (str): Search term for keys and values in the database.")
    print("    restrict_to_key (str, optional): Limit search to this specific keyname.")
    print("    restrict_to_keyvalue (str, optional): Restrict search to data_handles with this key=value.")
    print()
    print("Examples:")
    print("    searchDatabase('cancer')")
    print("    searchDatabase(query='cancer')")
    print("    searchDatabase(query='cancer', restrict_to_key='experiment.disease')")
    print("    searchDatabase(query='cancer', restrict_to_key='celltype', restrict_to_keyvalue='experiment.species=human')")
    print("    searchDatabase(query='cancer', restrict_to_key='celltype', restrict_to_keyvalue='species=human')")





def prepareData(data_handle, GENES):
    """
    Prepares RNA-seq data from a Seurat object and processes gene expression data.

    Parameters:
    - data_handle: The key to retrieve files from data_utils.
    - GENES: List of genes to process.

    Expected fields in getInfo object:
    - colname_celltype: Column name for cell type classification.
    - colname_celltype_value_treg: Specific cell type value to identify Treg cells.
    - colname_tissue: Column name for tissue information.
    - colname_treg_subtype: Column name for Treg subclustering.
    - colname_group: Column name for group information.
    - script_path: Path to the script for running Seurat object gene extraction.
    - default_genes: List of default genes to be included if missing from GENES.

    Usage:
    df = prepareData(
      data_handle='CANCER_2020_32603654_LUOMA', 
      GENES=['CCR8', 'FOXP3', 'CD3E', 'CD4', 'CD8A', 'IFNG', 'IL17A', 'CTLA4', 'PDCD1', 'CD27', 'IKZF2', 'TNFRSF18', 'CD177', 'TNFRSF8', 'CD74', 'FUT7', 'IL1R2', 'IL12RB2', 'TNFRSF4', 'CXCR3', 'STAT1', 'IL10'], 
      colname_celltype='authors_cell_type_level_2__cd4_t_cell', 
      colname_celltype_value_treg='7__regulatory_t_cell_treg', 
      colname_tissue='tissue_column_name', 
      colname_treg_subtype='authors_cell_type_level_3__treg_subclustering', 
      colname_group='group_column_name'
    )

    """
    
    script_path = os.path.expanduser('~/bin')

    # Check if the data_handle exists in data_utils
    if data_handle not in database:
        raise ValueError(f"Data handle '{data_handle}' not found in data_utils.")

    d_info=getInfo(data_handle)
    
    # Ensure the required files exist in the data_handle
    required_files = ['files.processed.meta', 'files.processed_dir', 'files.processed.rds']
    missing_files = [f for f in required_files if f not in d_info]
    if missing_files:
        raise ValueError(f"Missing files for data_handle '{data_handle}': {', '.join(missing_files)}")
    
    # Ensure the gene list is valid
    if not GENES or not isinstance(GENES, list):
        raise ValueError("GENES must be a non-empty list of gene names.")
    
    colname_celltype = d_info.get('colname_celltype',None)
    if not colname_celltype:
        raise ValueError("getInfo doesn't define required parameter: colname_celltype.")
    
    colname_celltype_value_treg = d_info.get('colname_celltype_value_treg', None)
    if not colname_celltype_value_treg:
        raise ValueError("getInfo doesn't define required parameter: colname_celltype_value_treg.")
    
    colname_tissue = d_info.get('colname_tissue', None)
    if not colname_tissue:
        raise ValueError("getInfo doesn't define required parameter: colname_tissue.")
    
    colname_treg_subtype = d_info.get('colname_treg_subtype', None)
    if not colname_treg_subtype:
        raise ValueError("getInfo doesn't define required parameter: colname_treg_subtype.")
    
    colname_group = d_info.get('colname_group')
    if not colname_group:
        raise ValueError("getInfo doesn't define required parameter: colname_group.")
    
    # Run external Seurat processing script
    RDS = d_info['files.processed.rds']
    DIR = d_info['files.processed_dir']    
    os.chdir(DIR)
    os.system(f'{script_path}/seurat_object_get_data_one_gene.sh {RDS} True normalized {" ".join(GENES)}')
    
    # Load and clean data
    df = pd.read_csv(d_info['files.processed.meta'], sep='\t')
    df.columns = utils.clean_str_values(df.columns)
    df.rename({'unnamed__0': 'barcode'}, axis='columns', inplace=True)
    df.set_index('barcode', inplace=True)
    
    # Ensure required columns exist in the meta file
    required_columns = [colname_celltype, colname_tissue, colname_treg_subtype, colname_group]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in 'processed.meta': {', '.join(missing_columns)}.")
    
    # Add positivity columns dynamically for each gene
    for gene in GENES:
        gene_lower = gene.lower()  # Convert gene name to lowercase
        if gene_lower in df.columns:
            df[f'is_{gene_lower}_positive'] = (df[gene_lower] > 0).astype(int)
        else:
            logging.warning(f"Gene '{gene_lower}' not found in the data columns. Skipping.")
    
    # Clean string columns
    for c in df.columns:
        if df[c].dtype == 'object':
            df[c] = utils.clean_str_values(df[c])

    # Identify Treg cells
    df['nb_celltype']=df[colname_celltype]
    df['is_treg'] = 0
    IDX = df[colname_celltype] == colname_celltype_value_treg
    df.loc[IDX, 'is_treg'] = 1
    
    # Add additional columns
    df['nb_tissue'] = df[colname_tissue]
    df['nb_treg_subtype'] = df[colname_treg_subtype]
    df['nb_group'] = df[colname_group]
    
    # Save the final DataFrame (but what if I need to make further modifications?)

    
    #if d_info.get('files.analysis_ready', None):
    #    df.to_csv(d_info['files.analysis_ready'], index=True, header=True, sep='\t')
    #else:
    #    raise ValueError(f'getInfo for {data_handle}, needs to have "files.analysis_ready" parameter defined')    
    #print(f'Created analysis ready data: {d_info["files.analysis_ready"]}. Load it via data_utils.getInfo(data_handle)["files.analysis_ready"]')

    saveAnalysisReadyData(df, data_handle)

def getAnalysisReadyData(data_handle):
    d_info=getInfo(data_handle)
    if d_info.get('files.analysis_ready', None):
        df=pd.read_csv(d_info['files.analysis_ready'], sep='\t')
        return(df)
    else:
        raise ValueError(f'getInfo for {data_handle}, needs to have "files.analysis_ready" parameter defined')
    

def saveAnalysisReadyData(df, data_handle):
    d_info=getInfo(data_handle)
    if d_info.get('files.analysis_ready', None):
        df.to_csv(d_info['files.analysis_ready'], index=True, header=True, sep='\t')
        print(f'[Message] Analysis ready file saved {d_info["files.analysis_ready"]}')
    else:
        raise ValueError(f'getInfo for {data_handle}, needs to have "files.analysis_ready" parameter defined')
    
