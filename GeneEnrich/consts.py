"""Constants used in GeneEnrich."""

# GeneEnrich supports Differential Expression Genes file from Scanpy or Seurat(v3 & v4)
SCANPY_TITLE = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'pct_nz_group', 'pct_nz_reference']
SEURATV3_TITLE = ['gene_id', 'p_val', 'avg_logFC' ,'pct.1', 'pct.2', 'p_val_adj']

# Default P-value cutoff for genes in Differential Expression Genes file
DEG_P_VALUE_CUTOFF = 0.05

# Default MIN PCT cutoff for genes in Differential Expression Genes file
DEG_MIN_PCT_CUTOFF = 0.1

# Default LOG FOLDCHANGES cutoff for genes in Differential Expression Genes file
DEG_LOG_FC_CUTOFF = 0.25
