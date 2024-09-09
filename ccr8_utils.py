import os
HOME=f'{os.path.expanduser("~")}'
import sys
sys.path.insert(0,f'{HOME}/lib')
import utils

def doublePostiveStats(df, gene_list, data_handle):
    gene_list=[g.lower() for g in gene_list]
    # Check if all genes are in the DataFrame
    missing_genes = [gene for gene in gene_list if f'is_{gene}_positive' not in df.columns]
    if missing_genes:
        raise ValueError(f"These genes are missing in the DataFrame: {', '.join(missing_genes)}")

    results = []

    # Loop over all unique tissues and groups in the DataFrame
    for tissue in df['nb_tissue'].unique():
        for group in df['nb_group'].unique():
            for gene in gene_list:
                IDX = (df['is_treg'] == 1) & (df['nb_group'] == group) & (df['nb_tissue'] == tissue)
                T = utils.tabulate(df[IDX], 'is_ccr8_positive', f'is_{gene}_positive', normalize_by='total')*100
                print(tissue, group, gene, T)
                T.index = ['ccr8_neg', 'ccr8_pos']
                T.columns = [f'{gene}_neg', f'{gene}_pos']
                T.index.name = 'ccr8_status'
                
                Tm = T.reset_index().melt(id_vars=['ccr8_status'])
                Tm['tissue'] = tissue
                Tm['group'] = group
                Tm['reference'] = data_handle
                results.append(Tm)

    # Concatenate results for all genes, groups, and tissues
    df_stats= pd.concat(results, ignore_index=True)
    df_stats_dp=df_stats[(df_stats.variable.str.endswith('_pos')) & (df_stats.ccr8_status.str.endswith('_pos'))]
    df_stats_dp['gene']=df_stats_dp.variable.str.removesuffix('_pos').str.upper()

    return(df_stats_dp)

def plotDoublePositiveStats(df_stats_dp):
    #input is a output of doublePostiveStats()
    g=p9.ggplot(df_stats_dp, p9.aes(x='gene', y='value', fill='group')) + p9.geom_bar(stat='identity', position='position_dodge') + p9.facet_grid('tissue~reference') + p9.ylab('Percent double positive with CCR8+')
    return(g)

def ccr8BispecificGenes():
    # TNFRSF8 (aka CD30)
    # FUT7 (aka CD15s)
    return 'CD177 TNFRSF8 CD74 FUT7 IL1R2'.split()
