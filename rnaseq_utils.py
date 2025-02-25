import re
import pandas as pd
from collections.abc import Iterable
import plotnine as p9
import numpy as np
    

def do_diffexp(df, groupA='control', groupB='colitis', subtypeA=['3__treg3','4__treg4'], subtypeB=['1__treg1'], col_subtype = 'authors_cell_type_level_3__treg_subclustering', RDS='luoma2020_tregs_cluster_bioturing.rds',THRESHOLD_PVALADJ=1e-2):
    '''
    Does DE via Seurat FindMarkers
    '''
    import os
    import plotnine as p9
    CWD='/Users/nbatada/data/rnaseq/CANCER_2020_32603654_LUOMA/bioturing/treg_cluster_over_allcells'
    os.chdir(CWD)

    # check that subtypeA and subtypeB are values in the col_subtype

    col_group='group'
    dfA=df[df[col_group]==groupA] 
    dfB=df[df[col_group]==groupB] 

    barcodesA= dfA[ dfA[col_subtype].isin(subtypeA)].reset_index()['barcode']
    barcodesB= dfB[ dfB[col_subtype].isin(subtypeB)].reset_index()['barcode']

    barcodesA.to_csv('barcodesA.txt', index=False, header=False)
    barcodesB.to_csv('barcodesB.txt', index=False, header=False)

    # run the diff exp
    results_filename='results_diffexpr_A_vs_B.tsv'
    os.system(f'/Users/nbatada/bin/find_diff_expr_genes.sh {RDS} barcodesA.txt barcodesB.txt {results_filename}')

    # read results
    df_deg=pd.read_csv(results_filename,sep='\t', header=0)
    
    # only keep significant ones
    THRESHOLD_PVALADJ=1e-2

    df_deg = df_deg[ (df_deg['p_val_adj'] < THRESHOLD_PVALADJ)]
    #g=p9.ggplot(df_deg, p9.aes(x='avg_log2FC', y='p_val_adj')) + p9.geom_point(alpha=0.1) + p9.xlim(-5,5) + p9.scale_y_log10()
    #print(g)
    os.system(f'/Users/nbatada/bin/plot_heatmap_of_degs.sh {RDS} barcodesA.txt barcodesB.txt {results_filename}')
    print('open heatmap_output.pdf')
    return(df_deg)



def volcano_plot(df, title=''):
    # Reset index to access gene names
    df = df.reset_index()
    
    # Calculate -log10(p_val_adj) for the y-axis
    df['neg_log10_pval_adj'] = -np.log10(df['p_val_adj'])
    
    # Get top 10 significant genes for avg_log2FC > 0
    top10_up = df[df['avg_log2FC'] > 0].nsmallest(10, 'p_val_adj')

    # Get top 10 significant genes for avg_log2FC < 0
    top10_down = df[df['avg_log2FC'] < 0].nsmallest(10, 'p_val_adj')

    # Combine both lists
    top_genes = pd.concat([top10_up, top10_down])

    # Determine y-axis limits
    ymin = 0
    ymax = np.ceil(df['neg_log10_pval_adj'].max())
    # Create the volcano plot
    g = (
        p9.ggplot(df, p9.aes(x='avg_log2FC', y='neg_log10_pval_adj')) +
        p9.geom_point(alpha=0.2) +
        p9.scale_x_continuous(limits=(-5, 5)) +
        p9.scale_y_continuous(limits=(ymin, ymax)) +
        p9.geom_vline(xintercept=1, linetype='dashed', color='red') +
        p9.geom_vline(xintercept=-1, linetype='dashed', color='red') +
        p9.geom_hline(yintercept=-np.log10(1e-2), linetype='dashed', color='blue') +
        p9.geom_text(p9.aes(label='index'), 
                     data=top_genes, 
                     nudge_y=0.2, 
                     size=8, 
                     ha='center')  +
        p9.ggtitle(f'{title}')
    )

    return g


