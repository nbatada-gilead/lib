import re
import pandas as pd
from collections.abc import Iterable
import plotnine as p9
import numpy as np

from datetime import datetime

import subprocess


def field_filter(columns, field_names, filters):
    """
    Filters a list of column names based on fields in their structure,
    with validations for field count and order consistency.

    Parameters:
    columns (list): List of column names to filter.
    field_names (pandas.Index or list): Field names corresponding to the fields in the column names.
    filters (dict): A dictionary where keys are field names and values are filtering criteria.

    Returns:
    list: Filtered list of column names that match the criteria.

    Usage:

    # Corrected function call
    filtered_columns = field_filter(columns=df_counts.columns, field_names=df_meta.columns, filters={'ccr8': 'pos'}    )

    where
    
    df_counts.columns
        Index(['1_neg_10000_d1_ss', '2_neg_10000_d2_ss', '3_neg_1000_d1_ds',
       '4_neg_1000_d2_ds', '5_neg_5000_d1_ds', '6_neg_5000_d1_ss',
       '7_neg_5000_d2_ds', '8_neg_5000_d2_ss', '9_pos_10000_d1_ss',
       '10_pos_10000_d2_ss', '11_pos_1000_d1_ds', '12_pos_1000_d2_ds',
       '13_pos_5000_d1_ds', '14_pos_5000_d1_ss', '15_pos_5000_d2_ds',
       '16_pos_5000_d2_ss'],
      dtype='object')

    df_meta.columns
        Index(['ccr8', 'n_cells', 'donor', 'protocol'], dtype='object')
    
    """
    # Convert field_names to a list if it's a pandas.Index
    if isinstance(field_names, pd.Index):
        field_names = field_names.tolist()

    # Split a column name to determine the number of fields
    sample_fields = columns[0].split("_")

    # Validation 1: Check if the number of fields matches
    if len(field_names) != len(sample_fields):
        raise ValueError(
            f"Number of fields in field_names ({len(field_names)}) "
            f"does not match the number of fields in column names ({len(sample_fields)})."
        )

    # Validation 2: Check order consistency
    for field, filter_value in filters.items():
        if field not in field_names:
            raise ValueError(
                f"Field '{field}' in filters is not in field_names. Available fields are: {field_names}"
            )
        field_index = field_names.index(field)
        valid_values = set(column.split("_")[field_index] for column in columns)
        filter_values = set(filter_value.split("|"))

        if not filter_values.issubset(valid_values):
            raise ValueError(
                f"Filter values '{filter_value}' for field '{field}' are inconsistent with "
                f"the column structure. Valid values for this field are: {valid_values}"
            )

    # Generate regex pattern for filtering
    def generate_pattern():
        pattern_parts = []
        for idx, field in enumerate(field_names):
            if field in filters:
                # If a filter is provided, use it
                pattern_parts.append(f"({filters[field]})")
            else:
                # If no filter is provided, allow any value
                pattern_parts.append(r".*")
        return "_".join(pattern_parts)

    # Generate the pattern and filter columns
    pattern = generate_pattern()
    filtered_columns = [col for col in columns if pd.Series([col]).str.match(pattern).any()]
    
    return filtered_columns

def filter_rows(df, N=0, threshold_fraction=0.75):
    """
    Filter rows where values are greater than N in at least a specified fraction of columns.
    
    Parameters:
    - df (pd.DataFrame): DataFrame with numeric columns.
    - N (int, optional): Minimum value threshold. Default is 0.
    - threshold_fraction (float, optional): Fraction of columns that must exceed N. Default is 0.75.

    Returns:
    - pd.DataFrame: Filtered DataFrame.
    """
    min_nonzero_columns = int(threshold_fraction * df.shape[1])
    return df[df.gt(N).sum(axis=1) >= min_nonzero_columns]



def fields_reorder(column_or_string, sep='_', neworder=[]):
    neworder = [i - 1 for i in neworder]  # Convert 1-indexed to 0-indexed
    def reorder(string):
        parts = string.split(sep)
        if max(neworder) >= len(parts):
            raise IndexError(f"New order indices {neworder} exceed available fields in '{string}'")
        return sep.join([parts[i] for i in neworder])
    
    if isinstance(column_or_string, pd.Series):
        return column_or_string.apply(reorder)
    return reorder(column_or_string)

def field_insert(column_or_string, from_pattern, sep='-'):
    """
    Add a specified separator between groups matched by the from_pattern.

    Parameters:
    - column_or_string (str or pd.Series): Input string or pandas Series.
    - from_pattern (str): Regex pattern to match groups where the separator will be added.
    - sep (str): Separator to insert between the matched groups. Default is '-'.

    Returns:
    - str or pd.Series: Modified string or pandas Series with the separator added.

    # Single string input
    string = 'D1-Ds-1000Neg'
    result = add_separator(string, r'(\d)([A-Za-z])', sep='-')
    print(result)  # Output: 'D1-Ds-1000-Neg'
    
    # Pandas Series input
    data = {'colname': ['D1-Ds-1000Neg', 'D2-Df-5000Pos']}
    df = pd.DataFrame(data)
    df['modified_col'] = field_insert(df['colname'], r'(\d)([A-Za-z])', sep='-')
    print(df)
    
    """
    if isinstance(column_or_string, pd.Series):  # If input is a pandas column
        return column_or_string.str.replace(from_pattern, rf'\1{sep}\2', regex=True)
    elif isinstance(column_or_string, str):  # If input is a single string
        return re.sub(from_pattern, rf'\1{sep}\2', column_or_string)
    else:
        raise ValueError("Input must be a string or pandas Series.")

def suffix_remove(columns, suffix):
    return columns.str.rstrip(suffix) if isinstance(columns, pd.Index) else [col.rstrip(suffix) for col in columns]
    
def prefix_add(column, prefix='', mode='incremental'):
    """
    Add a prefix to a pandas DataFrame column.

    Parameters:
    - column (pd.Series): The pandas DataFrame column to modify.
    - prefix (str): The prefix string to add to each element.
    - mode (str, optional): Determines how the prefix is added.
        - 'fixed': Adds the same fixed prefix to all elements.
        - 'incremental': Adds an incremental index (e.g., 'prefix0_', 'prefix1_', ...) to each element.
        Defaults to 'fixed'.

    Returns:
    - pd.Series: A new Series with the prefix applied.

    Raises:
    - ValueError: If the mode is invalid.

    Examples:
    >>> import pandas as pd
    >>> df = pd.DataFrame({'colname': ['A', 'B', 'C']})
    >>> add_prefix_to_column(df['colname'], 'fixed_', 'fixed')
    0    fixed_A
    1    fixed_B
    2    fixed_C
    Name: colname, dtype: object

    >>> add_prefix_to_column(df['colname'], '*_', 'incremental')
    0    *_0_A
    1    *_1_B
    2    *_2_C
    dtype: object
    """
    if mode == 'fixed':
        return prefix + column.astype(str)  # Add fixed prefix
    elif mode == 'incremental':
        return [f"{i+1}_" + str(val) for i, val in enumerate(column)]
    else:
        raise ValueError("Invalid mode. Use 'fixed' or 'incremental'.")

def fields_rename(strings, position, from_value, to_value, sep='_'):
    if isinstance(strings, str):  # Single string input
        parts = strings.split(sep)
        if parts[position] == from_value:
            parts[position] = to_value
        return sep.join(parts)
    elif isinstance(strings, list):  # List of strings input
        return [rename_fields(s, position, from_value, to_value, sep) for s in strings]
    else:  # Pandas column input
        return strings.apply(lambda s: rename_fields(s, position, from_value, to_value, sep))

def fields_swap(strings, from_position, to_position, sep='_'):
    if isinstance(strings, str):  # Single string input
        parts = strings.split(sep)
        parts[from_position], parts[to_position] = parts[to_position], parts[from_position]
        return sep.join(parts)
    elif isinstance(strings, list):  # List of strings input
        return [swap_fields(s, from_position, to_position, sep) for s in strings]
    else:  # Pandas column input
        return strings.apply(lambda s: swap_fields(s, from_position, to_position, sep))

    
def run_cmd_capture_output(command):
    try:
        result = subprocess.run(command, shell=True, text=True, capture_output=True, check=True)
        return result.stdout  # Capture and return stdout as a string
    except subprocess.CalledProcessError as e:
        print("Error:", e.stderr)  # Print stderr if an error occurs
        return None


def current_time():
    now = datetime.now()
    #current_time_str = now.strftime("%H:%M:%S (%Y-%m-%d)")
    current_time_str = now.strftime("%H:%M (%m/%d)")
    return(current_time_str)

def tabulate(df, r, c, normalize_by=None):
    '''
    r and c are column names in df
    in the output: r-values will be in rownames/index and c-values will be colnames
    
    USAGE: tabulate(df_tregs, 'group', 'authors_cell_type_level_3__treg_subclustering')

    normalize_by=None, rowsum or colsum

    '''
    u = df.groupby([r, c]).size().reset_index(name='count') # long
    
    T = u.pivot(index=r, columns=c, values='count') # compact r x c 
    
    #T['rowsum'] = T.sum(axis=1).astype('int')
    rowsum = T.sum(axis=1) #.astype('int')
    
    #T.loc['colsum'] = T.sum(axis=0).astype('int')
    colsum = T.sum(axis=0) #.astype('int')


    T.index.name = None
    T.columns.name = None
    T=T.fillna(0)
    if (normalize_by=='rowsum') | (normalize_by=='row'):
        #T=x.div(T[normalize_by], axis=0)*100
        T=T.div(rowsum, axis=0)*100
    elif (normalize_by=='colsum')| (normalize_by=='column'):
        #T=x.div(T.loc[normalize_by], axis=1)*100
        T=T.div(colsum, axis=1)*100
    elif (normalize_by=='total'):
        T = T / T.to_numpy().sum()
    return T


def clean_str_values(x):
    '''
    Cleanup names from strings such as column names of a dataframe or values within a column of a dataframe.

    Usage:
          df.columns = utils.clean_str_values(df.columns)
          df[c] = utils.clean_str_values(df[c])

    '''
    if isinstance(x, pd.Series) or isinstance(x, pd.Index):
        return x.str.lower().str.replace(' ', '_').str.replace("'", '').str.replace('-', '').str.replace('(', '').str.replace(')', '').str.replace('.', '_').str.replace(':', '_')
    elif isinstance(x, list):
        return [str(item).lower().replace(' ', '_').replace("'", '').replace('-', '').replace('(', '').replace(')', '').replace('.', '_').replace(':', '_') for item in x]
    else:
        raise TypeError("Input must be a Pandas Series, Index, or a list.")



def cmp_natural(text):
    '''
    compare function for sorting sample or gene names with (string+number )

    USAGE: alist.sort(key=cmp_natural) sorts in human order
    '''
    def atoi(text): 
        return int(text) if text.isdigit() else text
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def jupyter_header():
    return '''
import pandas as pd
pd.set_option('display.max_columns', 120)
pd.set_option('max_colwidth', 400)
pd.set_option('display.max_rows', 100)

import warnings
warnings.filterwarnings('ignore')

#--
import numpy as np
import plotnine as p9
from collections import Counter
from importlib import reload
    
#---
import sys
import os    
HOME=f'{os.path.expanduser("~")}'
sys.path.insert(0,f'{HOME}/lib')
import data_utils
import utils
import ccr8_utils
import plotnine as p9    
    '''
