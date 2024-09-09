import re
import pandas as pd
import importlib.util
import os
from pathlib import Path
import shutil
from datetime import datetime

database = {}
DATA_DIR = os.path.join(os.path.expanduser("~"), 'data', 'rnaseq')
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
    'files.processed.scrnaseq_obj': 18,
    'files.processed.source': 19,
    'files.processed.qc': 20,
    'files.processed.deg': 21,
    'files.processed.notes':22.
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
            for sub_key, value in sub_dict.items():
                if isinstance(value, str) and sub_key.startswith('files.'):
                    value = os.path.basename(value)
                line = f"{data_handle}.{sub_key}={value}\n"
                file.write(line)


def getInfo(data_handle, tsv_format=True):
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

    if tsv_format:
        df = pd.DataFrame(sorted_info, columns=['Key', 'Value'])
        df_transposed = df.set_index('Key').T
        df_transposed.index = [data_handle]
        df_transposed.columns.name = 'data_handle'
        return df_transposed
    else:
        for key, value in sorted_info:
            print(f"{key}\t{value}")


def deleteData(data_handle):
    if data_handle not in database:
        print(f"Data handle '{data_handle}' not found.")
        return

    non_empty_count = sum(1 for v in database[data_handle].values() if v.strip())

    if non_empty_count > 0:
        confirmation = input(f"The data handle '{data_handle}' contains {non_empty_count} non-empty values. Are you sure you want to delete it? (yes/no): ")
        if confirmation.lower() != 'yes':
            print("Deletion cancelled.")
            return

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file_path = f"{DATABASE_FILE_PATH}.{timestamp}.bak"
    shutil.copy(DATABASE_FILE_PATH, backup_file_path)

    del database[data_handle]
    __update_database()

    print(f"Data handle '{data_handle}' has been deleted and the database has been updated.")




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
    if len(result) == 1:
        result = result[0]  # just the string if it's a unique data_handle
    return result

def getFiles(data_handle):
    """Retrieve file paths for a specific data handle."""
    out = {}
    from pathlib import Path
    for key, value in database[data_handle].items():
        if key.startswith('files.'):
            if Path(value).suffix in {'.csv', '.tsv', '.txt'} and not value.startswith('http'):
                value = os.path.join(DATA_DIR, data_handle, value)

            key = key.removeprefix('files.')
            out[key] = value

    return out


def getInfo(data_handle):
    """Retrieve and format data from the database."""
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

    # Sort info according to DATA_ORDER
    sorted_info = sorted(info, key=lambda x: DATA_ORDER.get(x[0], float('inf')))

    df = pd.DataFrame(sorted_info, columns=['Key', 'Value'])
    df_transposed = df.set_index('Key').T
    df_transposed.index = [data_handle]
    df_transposed.columns.name = 'data_handle'
    return df_transposed


def newDataTemplate(data_handle):
    """Print a template for the specified data handle using DATA_ORDER."""
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

def addNewData(data_handle, **kwargs):
    """Add new information to the database for a given data handle."""
    global database

    if data_handle in database:
        print(f"Error: Data handle '{data_handle}' already exists in the database. \nUse updateInfo() instead.")
        return

    database[data_handle] = {}

    for key, value in kwargs.items():
        if key not in DATA_ORDER:
            print(f"Error: Key '{key}' is not valid. Valid keys are: {', '.join(DATA_ORDER.keys())}")
            continue

        database[data_handle][key] = value

    __update_database()
    print(f"[Message] New data added for data handle: '{data_handle}'")

def updateInfo(data_handle=None, help=False, **kwargs):
    """Update the information for a given data handle with key=value pairs."""
    if help:
        print("Usage: updateInfo(data_handle, key1=value1, key2=value2, ...)")
        print("Keys available for update:", '\n'.join(DATA_ORDER.keys()))
        return

    if not data_handle:
        print("Error: 'data_handle' must be provided unless 'help=True' is specified.")
        return

    if data_handle not in database:
        print(f"Error: Data handle '{data_handle}' not found in the database.")
        return

    for key, value in kwargs.items():
        if key not in DATA_ORDER:
            print(f"Error: Key '{key}' is not valid. Valid keys are: {', '.join(DATA_ORDER.keys())}")
            continue

        current_value = database[data_handle].get(key)
        if current_value == value:
            print(f"Info: No change needed for key '{key}' (value is already '{current_value}').")
            continue

        old_value = database[data_handle].get(key, 'Not set')
        database[data_handle][key] = value

        print(f"\nKey: {key}")
        print(f"Old value: {old_value}")
        print(f"New value: {value}\n")

    __update_database()
    print(f"[Message] Database updated for data handle: '{data_handle}'")



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

def search_database(query, restrict_to_key=None, restrict_to_keyvalue=None):
    """
    Search for data_handles in the database based on the query and optional restrictions.
    """
    if query.lower() == 'help':
        help_search_database()
        return {}

    result = {}
    query = query.lower()

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
        key_matches = []

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
        print(f"No data record found with for the query '{query}' and/or restrict_to_key '{restrict_to_key}' and/or restrict_to_keyvalue '{restrict_to_keyvalue}'")

    return result

def help_search_database():
    """Print usage information for the search_database function."""
    print("Usage:")
    print("    search_database(query, restrict_to_key=None, restrict_to_keyvalue=None)")
    print()
    print("Parameters:")
    print("    query (str): Search term for keys and values in the database.")
    print("    restrict_to_key (str, optional): Limit search to this specific keyname.")
    print("    restrict_to_keyvalue (str, optional): Restrict search to data_handles with this key=value.")
    print()
    print("Examples:")
    print("    search_database('cancer')")
    print("    search_database(query='cancer')")
    print("    search_database(query='cancer', restrict_to_key='experiment.disease')")
    print("    search_database(query='cancer', restrict_to_key='celltype', restrict_to_keyvalue='experiment.species=human')")
    print("    search_database(query='cancer', restrict_to_key='celltype', restrict_to_keyvalue='species=human')")

