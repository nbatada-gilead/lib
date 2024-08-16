# Nizar Batada
# 8 August 2024
# Gilead

# NOTE: This script relies on the associated data_utils_database.txt file where all the data resides

import re
import pandas as pd
import importlib.util
import os

DATA_DIR = os.path.join( os.path.expanduser("~"), 'data', 'rnaseq') 

#------------ This will be executed when module is loaded 

def __load_database(file_path):
    """Load the database dictionary from a two-column file and initialize the global dictionary."""
    global database
    database = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            if '=' not in line:
                continue
            
            # Split the line into key and value
            key, value = line.strip().split('=', 1)
            
            # Extract data_handle and nested keys
            parts = key.split('.')
            data_handle = parts[0]
            nested_key = '.'.join(parts[1:])
            
            # Initialize data_handle entry if not exists
            if data_handle not in database:
                database[data_handle] = {}
            
            # Initialize nested dictionary structure
            current_level = database[data_handle]
            for part in parts[1:-1]:
                if part not in current_level:
                    current_level[part] = {}
                current_level = current_level[part]
            
            # Set the final value
            current_level[parts[-1]] = value

# Initialize database with the file
DATABASE_FILE_PATH = os.path.join(os.path.dirname(__file__), 'data_utils_database.txt')
__load_database(DATABASE_FILE_PATH)

#------------------------------------------------------------
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
    'files.processed.scrnaseq_obj': 18,
    'files.processed.source': 19,
    'files.processed.qc': 20,
    'files.reprocessed.counts': 21,
    'files.reprocessed.normalized': 22,
    'files.reprocessed.meta': 23,
    'files.reprocessed.scrnaseq_obj': 24,
    'files.reprocessed.source': 25,
    'files.reprocessed.qc': 26,
    'manuscript.info': 27,
    'manuscript.doi': 28,
    'manuscript.url_supplementary_data': 29,
}



def listAvailableData():
    """Prints out the keys in the database dictionary."""
    print(list(database.keys()))

def queryData(query, keyword=None):
    """Search for data_handles in the database that meet certain criteria."""
    query = query.lower()
    result = []

    # Check if the query is a prefix for data_handle
    for data_handle in database:
        if query in data_handle.lower():
            result.append(data_handle)
        else:
            data = database[data_handle]
            if keyword and keyword in data:
                for v in data[keyword].values():
                    if query in str(v).lower():
                        result.append(data_handle)
                        break
            else:
                for subdict in data.values():
                    for v in subdict.values():
                        if query in str(v).lower():
                            result.append(data_handle)
                            break
    if len(result)==1:
        result=result[0] # just the string if it's a unique data_handle
    return result

def getFiles(data_handle, tsv=False):
    """Retrieve file paths for a specific data handle."""

    files_info = {}
    data = database.get(data_handle, {}).get('files', {})

    for category in ['raw', 'processed', 'reprocessed']:
        category_info = data.get(category, {})
        if category_info:
            files_info[category] = {
                'counts': category_info.get('counts', ''),
                'meta': category_info.get('meta', ''),
                'source': category_info.get('source', ''),
            }

    if tsv:
        # Convert to pandas DataFrame
        rows = []
        for category, content in files_info.items():
            for key, path in content.items():
                if path:  # Include only non-empty paths
                    rows.append([category, key, path])
        df = pd.DataFrame(rows, columns=['Category', 'Key', 'File Path'])
        return df
    else:
        # Create dictionary with file paths
        flattened_files = {}
        for category, content in files_info.items():
            for key, filename in content.items():
                if filename:  # Include only non-empty paths
                    #flattened_files[f"{category}.{key}"] = path
                    # if it's a note or a link, it will put in DATA_DIR to those as well
                    # flattened_files[f"{category}.{key}"] = os.path.join(DATA_DIR, filename)
                    flattened_files[f"{category}.{key}"] = filename
        flattened_files['DATA_DIR']=os.path.join(DATA_DIR, data_handle)
        return flattened_files
    
    
def getInfo(data_handle, tsv_format=True):
    """Retrieve and format data from the database."""
    info = []

    def flatten_dict(d, parent_key=''):
        for k, v in d.items():
            new_key = f"{parent_key}.{k}" if parent_key else k
            if isinstance(v, dict):
                flatten_dict(v, new_key)
            else:
                info.append((new_key, v))

    # Assuming 'database' is the dictionary holding the data
    data = database.get(data_handle, {})

    flatten_dict(data)

    # Sort info according to DATA_ORDER
    sorted_info = sorted(info, key=lambda x: DATA_ORDER.get(x[0], float('inf')))

    if tsv_format:
        # Convert to pandas DataFrame and transpose
        df = pd.DataFrame(sorted_info, columns=['Key', 'Value'])
        df_transposed = df.set_index('Key').T
        df_transposed.index = [data_handle]
        df_transposed.columns.name = 'data_handle'
        return df_transposed
    else:
        # Print in two-column format
        for key, value in sorted_info:
            print(f"{key}\t{value}")


def newDataTemplate(data_handle):
    """Print a template for the specified data handle using DATA_ORDER."""
    
    def print_nested_dict(prefix, d):
        """Recursively print nested dictionary keys with prefix."""
        for k, v in d.items():
            if isinstance(v, dict):
                # Recursive call for nested dictionaries
                print_nested_dict(f"{prefix}.{k}", v)
            else:
                # Print key-value pairs for non-dictionary values
                print(f"{prefix}.{k}=")
    
    for key in sorted(DATA_ORDER, key=DATA_ORDER.get):
        parts = key.split('.')
        if len(parts) > 1:
            # Build the key path from parts
            print_nested_dict(data_handle, {parts[0]: {k: '' for k in parts[1:]}})
        else:
            # Handle top-level keys
            print(f"{data_handle}.{key}=")

def help():
    import inspect
    
    # Get the current module
    module = inspect.getmodule(inspect.currentframe())
    
    # Get the names of all functions
    functions = [name for name, obj in inspect.getmembers(module, inspect.isfunction)]
    
    # Print exported functions with arguments and comments
    print("Exported Functions:")
    i=1
    for func_name in functions:
        if func_name.startswith('help') or func_name.startswith('__'): continue
        func = getattr(module, func_name)
        # Get function signature
        signature = inspect.signature(func)
        # Get the function's docstring (comment right after definition)
        docstring = inspect.getdoc(func)
        
        # Print function name and arguments
        args = [f"{param.name}={param.default}" if param.default != inspect.Parameter.empty else param.name for param in signature.parameters.values()]
        args_str = ", ".join(args)
        print(f"{i})  - {func_name}({args_str})")
        i+=1
        # Print the docstring if available
        if docstring:
            print(f"          {docstring}")
            print("\n")

# Add this to your module to use the help() function
if __name__ == "__main__":
    help()
    
