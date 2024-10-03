import re
import pandas as pd
from collections.abc import Iterable
import plotnine as p9
import numpy as np

from datetime import datetime

def current_time():
    now = datetime.now()
    #current_time_str = now.strftime("%H:%M:%S (%Y-%m-%d)")
    current_time_str = now.strftime("%H:%M (%m/%d)")
    return(current_time_str)

def dotplot(dfm, xcol='variable',ycol='value',sizecol='count'):
    # dfm = melted/long format
    g = (
        p9.ggplot(dfm, p9.aes(x=xcol, y=ycol, size=sizecol))
        + p9.geom_point(p9.aes(size=sizecol), color='blue', fill='lightblue', shape='o')
        + p9.scale_size(range = (1, 10))  # Adjust the range of point sizes
        + p9.theme(
            axis_text_x=p9.element_text(rotation=30, hjust=1),  
            figure_size=(3, 8) 
        )
    )
    return(g)

            
def barplot_rownames_as_xticks(T): 
    # Input: T is a tabulated data (i.e. utils.tabulate output)
    T_long=T.reset_index().melt(id_vars=['index'])

    T_long.columns=['colnames','rownames','value']

    FILL='rownames' # default incase it is not normalized

    # check both colsum and rowsum to determine what coloring to use in FILL
    if np.all(np.isin(T.sum(axis=0).astype('int'), [1, 100])): # could be 1 or 100
        FILL = 'colnames'
        #print('is normalized by column (down) sum')
    elif np.all(np.isin(T.sum(axis=1).astype('int'), [1, 100])):
        FILL = 'rownames'
        #print('is normalized by row (across) sum')
    g = (p9.ggplot(T_long, p9.aes(x='rownames', y='value', fill=FILL))
         + p9.geom_bar(stat='identity', position='dodge')
         + p9.facet_grid('colnames~') # alway colnames (stacked) because rownames is on xtick axis
         + p9.theme(axis_text_x=p9.element_text(rotation=25, hjust=1),
                    legend_position='none')  # Remove legend
         )
    
    return(g)


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
