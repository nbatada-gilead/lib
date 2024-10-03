# genemap_utils.py
# Nizar Batada
# 16 Aug 2024
# Module to provide annotations to genesymbols


import re
import pandas as pd
import os
from pathlib import Path
import shutil
from datetime import datetime

database = {}
DATA_DIR = os.path.join( os.path.expanduser("~"), 'data', 'geneinfo','genemaps') 
DATABASE_FILE_PATH = os.path.join(os.path.dirname(__file__), 'genemap_utils_database.txt')

#------------ This will be executed when module is loaded 

def __load_database(file_path=DATABASE_FILE_PATH):
    """Parse the database file and create a nested dictionary."""
    
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
            
            if data_handle not in database:
                database[data_handle] = {}
            database[data_handle][sub_key] = value

__load_database(DATABASE_FILE_PATH)

def __update_database(file_path=DATABASE_FILE_PATH):
    """Write the database dictionary back to the file in the original format."""
    
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
            if not isinstance(sub_dict, dict):
                raise ValueError(f"Expected dictionary for {data_handle}, got {type(sub_dict).__name__}.")
            
            for sub_key, value in sub_dict.items():
                if isinstance(value, str) and sub_key.startswith('files.'):
                    value = os.path.basename(value)
                
                line = f"{data_handle}.{sub_key}={value}\n"
                file.write(line)

# Data Order as in data_utils.py (you can adjust if needed)
DATA_ORDER = {
    'file.raw': 1,
    'file.map': 2,
    'source': 3,
    'species': 4,
    'script': 5,
    'info': 6
}

def listAvailableData():
    """Prints out the keys in the database dictionary."""
    return(list(database.keys()))

def _queryData(query, keyword=None):
    """Search for data_handles in the database that meet certain criteria."""
    query = query.lower()
    result = []

    for data_handle in database:
        if query in data_handle.lower():
            result.append(data_handle)
        else:
            data = database[data_handle]
            if keyword and keyword in data:
                if query in str(data[keyword]).lower():
                    result.append(data_handle)
                    break
    if len(result) == 1:
        result = result[0]
    return result
def queryData(query, keyword=None):
    """Search for data_handles in the database that meet certain criteria."""
    query = query.lower()
    result = []

    for data_handle, details in database.items():
        if query in data_handle.lower():
            result.append(data_handle)
        elif keyword:
            # Check if the keyword exists in the data
            value = details.get(keyword, '').lower()
            if query in value:
                result.append(data_handle)
        else:
            # Check all values if no specific keyword is given
            if any(query in str(value).lower() for value in details.values()):
                result.append(data_handle)

    if len(result) == 1:
        result = result[0]
    return result

def getFiles(data_handle):
    """Retrieve file paths for a specific data handle."""
    out = {}
    for key, value in database[data_handle].items():
        if key.startswith('file.'):
            if Path(value).suffix in {'.csv', '.tsv', '.txt','.map'} and not value.startswith('http'):
                value = os.path.join(DATA_DIR, data_handle, value)

            key = key.removeprefix('file.')
            out[key] = value
        
    return out



def addNewData(data_handle, file_raw, file_map, source='human protein atlas', species='human', script='create_map.py', create_dir=True):
    """Add new information to the database and save it back to file."""
    global database
    
    # Construct the folder path
    folder_path = os.path.join(DATA_DIR, data_handle)
    
    # Check if the folder exists
    if not os.path.exists(folder_path):
        if create_dir:
            # Create the folder if it doesn't exist and create_dir is True
            os.makedirs(folder_path)
            print(f'[Message] Created directory: {folder_path}')
        else:
            # Print a message and exit the function if create_dir is False
            print(f'[Warning] Did not update the database because the directory does not exist. Use create_dir=True to create the directory.')
            return
    
    # Add new data to the database
    if data_handle not in database:
        database[data_handle] = {}
    
    database[data_handle]['file.raw'] = file_raw
    database[data_handle]['file.map'] = file_map
    database[data_handle]['source'] = source
    database[data_handle]['species'] = species
    database[data_handle]['script'] = script

    # Update the database
    __update_database()
    print('[Message] Added new data: ', getInfo(data_handle))


def getInfo(data_handle, tsv=False):
    """Retrieve and format data from the database.
    tsv=True if you want to copy and paste into excel"""
    info = []

    def flatten_dict(d, parent_key=''):
        for k, v in d.items():
            new_key = f"{parent_key}.{k}" if parent_key else k
            if isinstance(v, dict):
                flatten_dict(v, new_key)
            else:
                info.append((new_key, v))

    data = database.get(data_handle, {})
    flatten_dict(data)
    sorted_info = sorted(info, key=lambda x: DATA_ORDER.get(x[0], float('inf')))

    if tsv:
        df = pd.DataFrame(sorted_info, columns=['Key', 'Value'])
        df_transposed = df.set_index('Key').T
        df_transposed.index = [data_handle]
        df_transposed.columns.name = 'data_handle'
        return df_transposed
    else:
        for key, value in sorted_info:
            print(f"{key}\t{value}")

def newDataTemplate(data_handle):
    """Print a template for the specified data handle using DATA_ORDER."""
    def print_nested_dict(prefix, d):
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
        if func_name.startswith('help') or func_name.startswith('__'):
            continue
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


def map(data_handle, input_data):
    """Map gene symbols in the input_data to the values in data_handle.
    
    Args:
        data_handle (str): The data handle to find the file.map path.
        input_data (list or pd.Series): List of gene symbols or a column from a DataFrame.

    Returns:
        list: A list of gene IDs or empty strings, matching the length of input_data.
    """
    # Ensure input_data is a pandas Series if it's a column from a DataFrame
    if isinstance(input_data, pd.Series):
        input_data = input_data.tolist()
    
    # Check if the data_handle exists and retrieve the map file path
    if data_handle not in database:
        raise ValueError(f"Data handle '{data_handle}' not found in the database.")
    
    d_files = getFiles(data_handle)
    
    if 'map' not in d_files:
        raise ValueError(f"'file.map' not found for data handle '{data_handle}' in the database.")
    
    map_file_path = d_files['map']
    
    # Load the mapping file
    try:
        map_df = pd.read_csv(map_file_path, sep='\t', header=None, names=['genesymbol', 'geneid'])
    except FileNotFoundError:
        raise FileNotFoundError(f"Mapping file '{map_file_path}' not found.")
    
    # Create a dictionary for quick lookup
    map_dict = pd.Series(map_df.geneid.values, index=map_df.genesymbol).to_dict()
    
    # Map input_data using the dictionary
    output = [map_dict.get(gene, '0') for gene in input_data]
    # always string because will be working with other categories of data
    return output

def updateInfo(data_handle=None, help=False, **kwargs):
    """
    Update the information for a given data_handle with key=value pairs.
    
    Args:
        data_handle (str): The handle to update in the database.
        help (bool): If True, provide usage information.
        **kwargs: Key-value pairs to update (must match keys in DATA_ORDER).
    
    Returns:
        None: Updates are saved directly to the database file.
    """
    if help:
        print("Usage: updateInfo(data_handle, key1=value1, key2=value2, ...)")
        print("Keys available for update:", ', '.join(DATA_ORDER.keys()))
        print("If no data_handle is provided, the function will print this help message.")
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
