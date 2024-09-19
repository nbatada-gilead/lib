import os
HOME=f'{os.path.expanduser("~")}'
import sys
sys.path.insert(0,f'{HOME}/lib')
import utils
import pandas as pd
import plotnine as p9

def doublePostiveStats(df, gene_list, data_handle):
    gene_list = [g.lower() for g in gene_list]

    # Check if all genes are in the DataFrame
    missing_genes = [gene for gene in gene_list if f'is_{gene}_positive' not in df.columns]
    if missing_genes:
        raise ValueError(f"These genes are missing in the DataFrame: {', '.join(missing_genes)}")

    results = []
    results_ccr8 = []  # For CCR8+ alone
    results_single_positive = []  # For single positive stats

    # Loop over all unique tissues and groups in the DataFrame
    for tissue in df['nb_tissue'].unique():
        for group in df['nb_group'].unique():
            # CCR8+ alone counts
            IDX_ccr8 = (df['is_treg'] == 1) & (df['nb_group'] == group) & (df['nb_tissue'] == tissue)
            ccr8_counts = df[IDX_ccr8]['is_ccr8_positive'].value_counts(normalize=True) * 100
            if 1 in ccr8_counts.index:
                results_ccr8.append({
                    'tissue': tissue,
                    'group': group,
                    'ccr8_positive_percentage': ccr8_counts[1],
                    'reference': data_handle
                })

            # Calculate single positive percentages for each gene and CCR8
            for gene in gene_list + ['ccr8']:
                column_name = f'is_{gene}_positive' if gene != 'ccr8' else 'is_ccr8_positive'
                gene_counts = df[IDX_ccr8][column_name].value_counts(normalize=True) * 100
                if 1 in gene_counts.index:
                    results_single_positive.append({
                        'tissue': tissue,
                        'group': group,
                        'gene': gene.upper(),
                        'positive_percentage': gene_counts[1],
                        'reference': data_handle
                    })

            # Double-positive calculation for each gene
            for gene in gene_list:
                IDX = (df['is_treg'] == 1) & (df['nb_group'] == group) & (df['nb_tissue'] == tissue)
                T = utils.tabulate(df[IDX], 'is_ccr8_positive', f'is_{gene}_positive', normalize_by='total') * 100

                try:
                    if T.shape == (2, 2):
                        T.index = ['ccr8_neg', 'ccr8_pos']
                        T.columns = [f'{gene}_neg', f'{gene}_pos']
                        T.index.name = 'ccr8_status'

                        Tm = T.reset_index().melt(id_vars=['ccr8_status'])
                        Tm['tissue'] = tissue
                        Tm['group'] = group
                        Tm['reference'] = data_handle
                        results.append(Tm)
                    else:
                        print(f"Skipping gene {gene} for group {group} and tissue {tissue} due to insufficient data.")
                except ValueError as e:
                    print(f"Error setting index/columns for gene {gene}, group {group}, tissue {tissue}")
                    print(f"T matrix: \n{T}")
                    raise e

    # Concatenate results for all genes, groups, and tissues
    df_stats = pd.concat(results, ignore_index=True)
    df_stats_dp = df_stats[(df_stats.variable.str.endswith('_pos')) & (df_stats.ccr8_status.str.endswith('_pos'))]
    df_stats_dp['gene'] = df_stats_dp.variable.str.removesuffix('_pos').str.upper()

    # CCR8+ alone statistics
    df_stats_ccr8_alone = pd.DataFrame(results_ccr8)

    # Single positive statistics
    df_stats_single_positive = pd.DataFrame(results_single_positive)

    # Normalize double positive stats by CCR8+ percentage
    df_stats_dp = df_stats_dp.merge(df_stats_ccr8_alone, on=['tissue', 'group', 'reference'], how='left')
    df_stats_dp['normalized_value'] = df_stats_dp['value'] / df_stats_dp['ccr8_positive_percentage']

    return df_stats_dp, df_stats_ccr8_alone, df_stats_single_positive, df_stats


def permutation_test_for_coexpression(group, is_gene1_present, is_gene2_present, n_permutations=1000):
    import numpy as np

    # observed_diff, p_value = permutation_test(group, is_gene1_present, is_gene2_present)
    # Calculate co-expression
    def calc_coexpression(g1, g2):
        return np.mean(np.array(g1) & np.array(g2))
    
    # Split the data into control and tumor groups
    control_idx = np.where(np.array(group) == 0)[0]
    tumor_idx = np.where(np.array(group) == 1)[0]
    
    control_coexp = calc_coexpression(np.array(is_gene1_present)[control_idx], 
                                      np.array(is_gene2_present)[control_idx])
    tumor_coexp = calc_coexpression(np.array(is_gene1_present)[tumor_idx], 
                                    np.array(is_gene2_present)[tumor_idx])
    
    # Observed difference in co-expression
    observed_diff = tumor_coexp - control_coexp
    
    # Permutation test
    combined = np.array(is_gene1_present) & np.array(is_gene2_present)
    perm_diffs = []
    
    for _ in range(n_permutations):
        perm_group = np.random.permutation(group)
        control_idx_perm = np.where(np.array(perm_group) == 0)[0]
        tumor_idx_perm = np.where(np.array(perm_group) == 1)[0]
        
        control_perm = calc_coexpression(combined[control_idx_perm], combined[control_idx_perm])
        tumor_perm = calc_coexpression(combined[tumor_idx_perm], combined[tumor_idx_perm])
        
        perm_diffs.append(tumor_perm - control_perm)
    
    perm_diffs = np.array(perm_diffs)
    p_value = np.mean(perm_diffs >= observed_diff)
    
    return observed_diff, p_value




def plotDoublePositiveStats(df_stats_dp, version=None):
    #input is a output of doublePostiveStats()
    #g=p9.ggplot(df_stats_dp, p9.aes(x='gene', y='value', fill='group')) + p9.geom_bar(stat='identity', position='position_dodge') + p9.facet_grid('tissue~reference') + p9.ylab('Percent double positive with CCR8+')
    
    if version:
        g=p9.ggplot(df_stats_dp, p9.aes(x='tissue', y='normalized_value', fill='tissue')) + p9.geom_bar(stat='identity', position='position_dodge') + p9.facet_grid('gene~group' ) + p9.ylab('Percent double positive with CCR8+') + p9.theme_classic() + p9.theme(axis_text_x=p9.element_text(rotation=30, hjust=1))
    else:
        g=p9.ggplot(df_stats_dp, p9.aes(x='reference', y='normalized_value', fill='group')) + p9.geom_bar(stat='identity', position='position_dodge') + p9.facet_grid('gene~tissue') + p9.ylab('Percent double positive with CCR8+') + p9.theme_classic() + p9.theme(axis_text_x=p9.element_text(rotation=30, hjust=1))
        
    return(g)

def ccr8BispecificGenes():
    # TNFRSF8 (aka CD30)
    # FUT7 (aka CD15s)
    return 'CD177 TNFRSF8 CD74 FUT7 IL1R2'.split()
