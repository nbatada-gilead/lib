import os
HOME=f'{os.path.expanduser("~")}'
import sys
sys.path.insert(0,f'{HOME}/lib')
import utils
import pandas as pd
import plotnine as p9

def singlePositiveStats(df, gene_list, data_handle):
    gene_list = [g.lower() for g in gene_list]

    # Check if all genes are in the DataFrame
    missing_genes = [gene for gene in gene_list if f'is_{gene}_positive' not in df.columns]
    if missing_genes:
        raise ValueError(f"These genes are missing in the DataFrame: {', '.join(missing_genes)}")

    results_single_positive = []

    # Loop over all unique tissues and groups in the DataFrame
    for tissue in df['nb_tissue'].unique():
        for group in df['nb_group'].unique():
            # Calculate single positive percentages for each gene
            for gene in gene_list:
                column_name = f'is_{gene}_positive'
                IDX = (df['is_treg'] == 1) & (df['nb_group'] == group) & (df['nb_tissue'] == tissue)
                gene_counts = df[IDX][column_name].value_counts(normalize=True) * 100
                if 1 in gene_counts.index:
                    results_single_positive.append({
                        'tissue': tissue,
                        'group': group,
                        'gene': gene.upper(),
                        'positive_percentage': gene_counts[1],
                        'reference': data_handle
                    })

    # Single positive statistics
    df_stats_single_positive = pd.DataFrame(results_single_positive)

    return df_stats_single_positive


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

def plotSinglePositiveStats(df_stats):
    import plotnine as p9
    g=p9.ggplot(df_stats, p9.aes(x='gene', y='positive_percentage', fill='group')) + p9.geom_bar(stat='identity', position='position_dodge')  + p9.facet_grid('.~tissue') + p9.theme(figure_size=(10,5), axis_text_x=p9.element_text(rotation=30, hjust=1))
    return(g)


compute_contingency_table <- function(sobj, gene) {
  # Ensure the gene exists in the Seurat object
  if (!gene %in% rownames(sobj@assays$RNA@data)) {
    stop(paste("Gene", gene, "not found in the Seurat object."))
  }
  
  # Map both "normal" and "control" to "control"
  sobj@meta.data <- sobj@meta.data %>%
    mutate(group_rds = ifelse(group_rds %in% c("normal", "control"), "control", group_rds))

  # Create the group_combined column if it doesn't exist
  if (!"group_combined" %in% colnames(sobj@meta.data)) {
    sobj@meta.data <- sobj@meta.data %>%
      mutate(CCR8_expr_group = ifelse(sobj@assays$RNA@data["CCR8", ] > 0, "CCR8_positive", "CCR8_zero"),
             group_combined = paste(group_rds, CCR8_expr_group, sep = "_"))
  }

  # Calculate gene positive based on expression > 0
  sobj@meta.data <- sobj@meta.data %>%
    mutate(gene_positive = ifelse(sobj@assays$RNA@data[gene, ] > 0, 1, 0))

  # Subset for CCR8 positive/zero and tumor/control cells
  tumor_ccr8_positive <- sobj@meta.data %>%
    filter(group_combined == "tumor_CCR8_positive")
  control_ccr8_positive <- sobj@meta.data %>%
    filter(group_combined == "control_CCR8_positive")
  
  tumor_ccr8_zero <- sobj@meta.data %>%
    filter(group_combined == "tumor_CCR8_zero")
  control_ccr8_zero <- sobj@meta.data %>%
    filter(group_combined == "control_CCR8_zero")

  # Create contingency table for the gene
  gene_tumor_positive <- sum(tumor_ccr8_positive$gene_positive == 1)
  gene_tumor_negative <- sum(tumor_ccr8_positive$gene_positive == 0)
  gene_control_positive <- sum(control_ccr8_positive$gene_positive == 1)
  gene_control_negative <- sum(control_ccr8_positive$gene_positive == 0)

  # Calculate percentages
  tumor_pct_positive <- gene_tumor_positive / (gene_tumor_positive + gene_tumor_negative) * 100
  control_pct_positive <- gene_control_positive / (gene_control_positive + gene_control_negative) * 100
  fold_change <- tumor_pct_positive / control_pct_positive
  
  # Combine into a contingency table
  gene_table <- matrix(c(gene_tumor_positive, gene_tumor_negative, gene_control_positive, gene_control_negative),
                       nrow = 2,
                       byrow = TRUE,
                       dimnames = list(c("Tumor", "Control"), c(paste0(gene, "_Positive"), paste0(gene, "_Negative"))))

  # Return results
  return(list(table = gene_table, tumor_pct_positive = tumor_pct_positive, 
              control_pct_positive = control_pct_positive, fold_change = fold_change))
}

# Given a gene, a stat test for increase will be done for CCR8+ tumor versus CCR8+ normal tregs (this will be from compute_contingency_table)
perform_stat_test <- function(sobj, gene) {
  # Compute the contingency table and percentages
  results <- compute_contingency_table(sobj, gene)
  gene_table <- results$table
  tumor_pct_positive <- results$tumor_pct_positive
  control_pct_positive <- results$control_pct_positive
  fold_change <- results$fold_change

  # Perform statistical test
  if (any(gene_table == 0)) {
    gene_test <- fisher.test(gene_table)
    test_type <- "Fisher's Exact Test"
  } else {
    gene_test <- chisq.test(gene_table)
    test_type <- "Chi-Square Test"
  }

  # Print results
  print(paste(gene, test_type))
  print(gene_test)
  print(paste("Tumor % positive:", tumor_pct_positive))
  print(paste("Control % positive:", control_pct_positive))
  print(paste("Fold change (tumor/control):", fold_change))

  # Return test results
  return(list(test = gene_test, tumor_pct_positive = tumor_pct_positive, 
              control_pct_positive = control_pct_positive, fold_change = fold_change))
}

# Given a gene,data from compute_contingency_table will be used to create a plot

create_pie_chart <- function(sobj, gene) {
    # Compute the contingency table and percentages
    results <- compute_contingency_table(sobj, gene)
    tumor_pct_positive <- results$tumor_pct_positive
    control_pct_positive <- results$control_pct_positive
    tumor_zero_pct <- 100 - tumor_pct_positive
    control_zero_pct <- 100 - control_pct_positive
    
    # Prepare data for pie chart
    pie_data <- data.frame(
        Group = factor(c("Tumor", "Tumor", "Control", "Control"),
                       levels = c("Control", "Tumor")),
        Status = c("Positive", "Negative", "Positive", "Negative"),
        Percent = c(tumor_pct_positive, tumor_zero_pct, control_pct_positive, control_zero_pct)
    )
    
    # Add labels for percentages, only for Positive slices
    pie_data$label <- ifelse(pie_data$Status == "Positive", paste0(round(pie_data$Percent, 1), "%"), "")
    
    # Create pie chart
    pie_chart <- ggplot(pie_data, aes(x = "", y = Percent, fill = Status)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar("y") +
        facet_wrap(~ Group, ncol = 2) +
        labs(title = paste(gene, "Co-expression Pie Chart (CCR8 Positive Only)"),
             fill = "Gene Expression") +
        geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
        theme_void() +
        scale_fill_manual(values = c("Positive" = "dodgerblue", "Negative" = "gray80"))
    
    # Print pie chart
    print(pie_chart)
    
    # Return plot object
    return(pie_chart)
}


