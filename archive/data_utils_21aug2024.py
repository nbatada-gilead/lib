# Nizar Batada
# 8 August 2024
# Gilead

# NOTE: This script relies on the associated data_utils_database.txt file where all the data resides

import re
import pandas as pd
import importlib.util
import os
from pathlib import Path
import shutil
from datetime import datetime

database={}
DATA_DIR = os.path.join( os.path.expanduser("~"), 'data', 'rnaseq') 
DATABASE_FILE_PATH = os.path.join(os.path.dirname(__file__), 'data_utils_database.txt')


#------------ This will be executed when module is loaded 

def __load_database(file_path=DATABASE_FILE_PATH):
    """Parse the database file and create a nested dictionary."""
    
    global database
    database = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Strip whitespace and ignore empty lines
            line = line.strip()
            if not line:
                continue

            # Split the line into key and value
            key, value = line.split('=', 1)
            
            # Extract the data_handle and sub_key
            try:
                data_handle, sub_key = key.split('.', 1)
            except ValueError:
                print(f"Skipping invalid line: {line}")
                continue
            
            # Ensure the database structure is correct
            if data_handle not in database:
                database[data_handle] = {}
            database[data_handle][sub_key] = value

# Initialize database with the file
__load_database(DATABASE_FILE_PATH)

## BE AWARE, IF THE DATABASE file is opened by someone else (say I'm adding additional infomration) and if I do
## __update_database, then the currently loaded informatio will be written to the file and all modificaitons will be lost.
## if I check that the file is not currently open and, if yes, ask the user to close it first before rerunning __update_database (which is done automatically by createFiles
##
import os
import shutil
from datetime import datetime

def __update_database(file_path=DATABASE_FILE_PATH):
    """Write the database dictionary back to the file in the original format."""

    # Check if the file is open
    if os.path.isfile(file_path):
        try:
            with open(file_path, 'r+'):
                pass
        except IOError:
            input(f"The file {file_path} is currently open. Please close it and press Enter to continue.")
    
    # Create a backup with a timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file_path = f"{file_path}.{timestamp}.bak"
    shutil.copy(file_path, backup_file_path)

    # Write the updated database to the original file
    with open(file_path, 'w') as file:
        for data_handle, sub_dict in database.items():
            if not isinstance(sub_dict, dict):
                raise ValueError(f"Expected dictionary for {data_handle}, got {type(sub_dict).__name__}.")
            
            # Ensure all keys in DATA_ORDER are present in the database
            for key in DATA_ORDER:
                if key not in sub_dict:
                    sub_dict[key] = ""

            for sub_key, value in sub_dict.items():
                if isinstance(value, str) and sub_key.startswith('files.'):
                    # Remove DATA_DIR and data_handle prefixes if necessary
                    value = os.path.basename(value)
                
                line = f"{data_handle}.{sub_key}={value}\n"
                file.write(line)


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
    'files.processed.deg':21,
    'files.reprocessed.counts': 22,
    'files.reprocessed.normalized': 23,
    'files.reprocessed.meta': 24,
    'files.reprocessed.scrnaseq_obj': 25,
    'files.reprocessed.source': 26,
    'files.reprocessed.qc': 27,
    'files.reprocessed.deg':28,
    'manuscript.info': 29,
    'manuscript.doi': 30,
    'manuscript.url_supplementary_data': 31,
    'manuscript.notes':32,
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
    if len(result)==1:
        result=result[0] # just the string if it's a unique data_handle
    return result


def getFiles(data_handle):
    """Retrieve file paths for a specific data handle."""

    out={}
    from pathlib import Path
    for key, value in database[data_handle].items():
        if key.startswith('files.'):
            if Path(value).suffix in {'.csv','.tsv','.txt'} and not value.startswith('http'):
                value = os.path.join(DATA_DIR, data_handle, value)

            key=key.removeprefix('files.')
            out[key]=value
        
    return(out)

    
def createFileName(data_handle): 
    ''' Will update database, will write to the database.txt'''

    UPDATE_DATABASE=False
    files=getFiles(data_handle)
    if files['reprocessed.counts']=='':
        database[data_handle]['files.reprocessed.counts']= os.path.join(DATA_DIR, data_handle, data_handle + '_reprocessed_expr_counts.csv')
        UPDATE_DATABASE=True
    else:
        print('reprocessed.counts file already exists: ', files['reprocessed.counts'])
    if files['reprocessed.meta']=='':
        database[data_handle]['files.reprocessed.meta']= os.path.join(DATA_DIR, data_handle, data_handle + '_reprocessed_meta.csv')
        UPDATE_DATABASE=True
        
    else:
        print('reprocessed.meta file already exists: ', files['reprocessed.meta'])
    if UPDATE_DATABASE:
        # BUT I NEED TO NOT ADD THE FULL PATH
        __update_database()
    
    # write the database to the file
    
    
    
    
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

def addNewData(data_handle, **kwargs):
    """
    Add new information to the database for a given data handle.

    Args:
        data_handle (str): The handle to add to the database.
        **kwargs: Key-value pairs to initialize (must match keys in DATA_ORDER).
    
    Returns:
        None: Updates are saved directly to the database file.
    """
    global database

    # Check if data_handle exists in the database
    if data_handle in database:
        print(f"Error: Data handle '{data_handle}' already exists in the database. \nUse updateInfo() instead. ")
        return

    # Create an entry in the database
    database[data_handle] = {}

    for key, value in kwargs.items():
        if key not in DATA_ORDER:
            print(f"Error: Key '{key}' is not valid. Valid keys are: {', '.join(DATA_ORDER.keys())}")
            continue
        
        database[data_handle][key] = value

    # Update the database
    __update_database()
    print(f"[Message] New data added for data handle: '{data_handle}'")


def updateInfo(data_handle=None, help=False, **kwargs):
    """
    Update the information for a given data handle with key=value pairs.

    Args:
        data_handle (str): The handle to update in the database.
        help (bool): If True, provide usage information.
        **kwargs: Key-value pairs to update (must match keys in DATA_ORDER).
    
    Returns:
        None: Updates are saved directly to the database file.
    """
    if help:
        print("Usage: updateInfo(data_handle, key1=value1, key2=value2, ...)")
        print("Keys available for update:", '\n'.join(DATA_ORDER.keys()))
        return

    if not data_handle:
        print("Error: 'data_handle' must be provided unless 'help=True' is specified.")
        return

    # Check if data_handle exists in the database
    if data_handle not in database:
        print(f"Error: Data handle '{data_handle}' not found in the database.")
        return

    # Iterate over the provided key-value pairs
    for key, value in kwargs.items():
        if key not in DATA_ORDER:
            print(f"Error: Key '{key}' is not valid. Valid keys are: {', '.join(DATA_ORDER.keys())}")
            continue

        # Check if the value is actually being updated
        current_value = database[data_handle].get(key)
        if current_value == value:
            print(f"Info: No change needed for key '{key}' (value is already '{current_value}').")
            continue

        # Update the database
        old_value = database[data_handle].get(key, 'Not set')
        database[data_handle][key] = value

        # Print the update message
        print(f"\nKey: {key}")
        print(f"Old value: {old_value}")
        print(f"New value: {value}\n")
    
    # Write changes to the database file
    __update_database()
    print(f"[Message] Database updated for data handle: '{data_handle}'")
    


    
