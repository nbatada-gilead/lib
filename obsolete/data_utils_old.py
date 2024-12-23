import os
import re
from tabulate import tabulate

DATA_PATH = os.path.join( os.path.expanduser("~"), 'data', 'rnaseq') 
DATA_ORDER={
}

def summarizeData(data_handle):
    '''
    prints all the information in data_handles as tab separated row with fixed column order and information
    '''
    pass

database = {
    'd07_magnuson2018_30348759':{
        'experiment':{
            'species':'human',            
            'tissue':'colon',
            'disease_type':'cancer',
            'disease_list':'crc|normal'
            'celltype':'tregs',
            'design':'',
            'nsamples':'6 patients with healthy and tumor tissue + few more unmatched',
            
        },
        'assay':{            
            'type' :'bulk_sorted',
            'platform':'rnaseq',
            'ncells':'',
        },
        'files':{
            'raw':{
                'counts':'GSE116347_Gene_count_table_human.csv',
                'meta':'column header of files_raw expr_counts',
                'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116347',
            },
            'processed':{
                'counts':'',
                'normalized':'',
                'meta':'',
                'scrnaseq_obj':'',
                'source':'',
                'qc':'',
            },
            'reprocessed':{
                'counts': '',
                'normalized':'',
                'meta':'',
                'scrnaseq_obj':'',
                'source':'',
                'qc':'',
            },
        },
        'manuscript':{
            'info':'Magnuson et al. Identification and validation of a tumor-infiltrating Treg transcriptional signature conserved across species and tumor types. Proc Natl Acad Sci U S A. 2018.',
            'doi':'DOI: 10.1073/pnas.1810580115',
            'url_supplementary_data':'',
        },
    },
    
    'd06_elmentaite2021_34497389': {
        'keywords':['teichmann','unsorted','intestinal','healthy'], 
        'info':'https://www.nature.com/articles/s41586-021-03852-1 Full single cell RNA-seq dataset of 428K intestinal cells from fetal, pediatric, adult donors, and up to 11 intestinal regions.   https://www.gutcellatlas.org/',
        'data_type':'scrnaseq(10x)', 
        'tissue':'intestinal', 
        'experimental_design':'Human adult tissue was obtained by the Cambridge Biorepository of Translational Medicine from deceased transplant organ donors after ethical approval (reference 15/EE/0152, East of England - Cambridge South Research Ethics Committee) and informed consent from the donor families. Details of the ages and genders of donors are included as Supplementary Table 2. Samples were collected from 11 distinct locations, including the duodenum (two locations, DUO1 and DUO2, which were pooled in the analysis), jejunum (JEJ), ileum (two locations, ILE1 and ILE2, which were pooled in the analysis), appendix (APD), caecum (CAE), ascending colon (ACL), transverse colon (TCL), descending colon (DCL), sigmoid colon (SCL), rectum (REC) and mesenteric lymph nodes (mLN). Fresh mucosal intestinal tissue and lymph nodes from the intestinal mesentery were excised within 1 h of circulatory arrest; intestinal tissue was preserved in University of Wisconsin organ-preservation solution (Belzer UW Cold Storage Solution; Bridge to Life) and mLN were stored in saline at 4C until processing. Tissue dissociation was conducted within 2 h of tissue retrieval.',
        'disease_type':'healthy',
        'num_samples': 0,
        'celltype': 'unsorted',
        'data_url': {
            'url_expr': '',
            'url_meta': '', 
            'url_processed': 'https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Tcell_log_counts02_v2.h5ad'
            }
    },
    'd05_plitas2016_27851913': {
        'keywords':['bulk','tregs','brca','healthy','ccr8'], 
        'info':'Regulatory T Cells Exhibit Distinct Features in Human Breast Cancer Immunity 2016 https://www.ncbi.nlm.nih.gov/pubmed/27851913', 
        'data_type':'bulk_rnaseq', 
        'tissue':'primary_brca, healthy_breast_tissue', 
        'disease_type':'tumor_brca',
        'num_samples': '',
        'celltype': 'Tregs',
        'data_url': {
            'url_expr': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89225/suppl/GSE89225%5FIllumina%5Fcounts.csv.gz',
            'url_meta': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89225/matrix/GSE89225-GPL16791_series_matrix.txt.gz', 
            'url_misc': ''
            }
    },
    'd04_tabulasampies_35549404':{
        'keywords':['normal','tabulasapiens','human'],
        'info':'The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. Science 2022 500,000 cells from 24 different tissues and organs, many from the same donor.',
        'data_type': 'singlecell_rna',
        'tissue': 'multiple',
        'disease_type': 'normal',
        'num_samples': '',
        'celltype': 'multiple',
        'data_url': {'url_expr': '', 'url_meta': '', 'url_misc': ''}
    },
    'd03_yang2024_39047727': {
        'keywords':['pancancer','tumor','bcells','b cells','zemin','d03','yang2024','39047727'], 
        'info':'B cell atlas from tumor: Pan-cancer single-cell dissection reveals phenotypically distinct B cell subtypes. Cell 2024 DOI: 10.1016/j.cell.2024.06.038',
        'data_type': 'singlecell_rna',
        'tissue': 'multiple',
        'disease_type': 'cancer',
        'num_samples': '',
        'celltype': 'B cells',
        'data_url': {'url_expr': '', 'url_meta': '', 'url_misc': ''}
    },
    'd02_miragaia2019_30737144': {
        'keywords':['teichmann','colon','healthy','skin','pbmc','ss2'], 
        'info':'Single-Cell Transcriptomics of Regulatory T Cells Reveals Trajectories of Tissue Adaptation. Immunity 2019 DOI: 10.1016/j.immuni.2019.01.001',
        'data_type': 'singlecell_rna',
        'tissue': 'colon',
        'disease_type': 'normal',
        'num_samples': '',
        'celltype': 'Regulatory T cells',
        'data_url': {'url_expr': '', 'url_meta': '', 'url_misc': ''}
    },
    'd01_zheng2019_34914499': {
        'keywords': ['zheng2019', 'd01', '34914499', 'cd3', 'pancancer'],
        'info': 'Link to manuscript: https://example.com/zheng2019',
        'data_type': 'singlecell_rna',
        'tissue': 'multiple',
        'disease_type': 'cancer',
        'num_samples': '',
        'celltype': 'T cells',
        'data_url': {'url_expr': '', 'url_meta': '', 'url_misc': ''}
    }
}

def listAvailableData():
    '''
    List available data handles based on predefined information.

    Returns:
        list: A list of available data handles as strings.
    '''
    return list(data_handles_info.keys())


def queryData(keyword=None, data_type=None, tissue=None, disease_type=None):
    """
    Query data handles based on a keyword, data_type, tissue, and disease_type.

    Args:
        keyword (str): The keyword or pattern to search for.
        data_type (str): The type of data (e.g., 'singlecell_rna', 'bulk_rna').
        tissue (str): The tissue type (e.g., 'lung', 'colon').
        disease_type (str): The disease type (e.g., 'normal', 'cancer').

    Returns:
        list: A list of data handles that match the query parameters.
    """
    matches = []
    for handle, info in data_handles_info.items():
        if keyword:
            pattern = re.compile(re.escape(keyword), re.IGNORECASE)
            if not (pattern.search(handle) or any(pattern.search(kw) for kw in info["keywords"])):
                continue
        if data_type and info.get('data_type') != data_type:
            continue
        if tissue and info.get('tissue') != tissue:
            continue
        if disease_type and info.get('disease_type') != disease_type:
            continue
        matches.append(handle)
    return matches


def getFilePath(data_handle):
    '''
    Get the file paths for the given data handle.

    Args:
        data_handle (str): The data handle string.

    Returns:
        tuple: A tuple containing the paths to the expression and metadata files.
    '''
    if data_handle in data_handles_info:
        expr_file = os.path.join(DATA_PATH, data_handle, f'{data_handle}_expr.tsv')
        meta_file = os.path.join(DATA_PATH, data_handle, f'{data_handle}_meta.tsv')
        return {'expr_file':expr_file, 'meta_file':meta_file}
    else:
        raise ValueError(f"Data handle '{data_handle}' not found in predefined information.")


def getInfo(data_handle):
    '''
    Get detailed information about a data handle, excluding keywords, formatted as a table.

    Args:
        data_handle (str): The data handle string.

    Returns:
        str: Detailed information about the data handle formatted as a table.
    '''
    handle_info = data_handles_info.get(data_handle, {})
    if not handle_info:
        return 'No information available.'
    
    # Exclude the 'keywords' key from the information
    info_excluding_keywords = {key: value for key, value in handle_info.items() if key != 'keywords'}
    
    # Format the information as a table
    table = tabulate(info_excluding_keywords.items(), headers=['Field', 'Information'])
    
    return table
    

def getDataUrl(data_handle):
    '''
    Get the data URLs for the given data handle.

    Args:
        data_handle (str): The data handle string.

    Returns:
        dict: A dictionary containing the URLs for expression, metadata, and miscellaneous files.
    '''
    handle_info = data_handles_info.get(data_handle, {})
    if not handle_info:
        return 'No data URL available.'
    
    return handle_info.get('data_url', {})
