from GeneEnrich.enricher import enricher_internal
from GeneEnrich.preprocess import get_gene_info
from GeneEnrich.run_kegg.kegg_utilities import get_pathway_mapping, get_gene_pathway_mapping
import pandas as pd
import sys


def enrichkegg(
    gene,
    organism="hsa",
    pvalueCutoff=0.05,
    pAdjustMethod="fdr_bh",
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff=0.2,
):

    pathway_name_df, gene_pathway_mapping_df = get_data_from_KEGG_db(organism)

    res = enricher_internal(gene,
                            pvalueCutoff=pvalueCutoff,
                            pAdjustMethod=pAdjustMethod,
                            minGSSize=minGSSize,
                            maxGSSize=maxGSSize,
                            qvalueCutoff=qvalueCutoff,
                            pathway_name_df=pathway_name_df,
                            gene_pathway_mapping_df=gene_pathway_mapping_df)

    if res is not None and not res.empty:
        gene2convert = {geneid: genesymbol for genesymbol, geneid in get_gene_info().items()}
        res['geneID'] = res['geneID'].apply(lambda x: '/'.join([gene2convert.get(i, i) for i in x.split('/')]))
    else:
        sys.exit('No enrichment analysis result can be found, check input file')

    return res


def get_data_from_KEGG_db(species: str):

    pathway_name = get_pathway_mapping(
        organism=species
    )
    pathway_name_df = pd.DataFrame({
        'path': pathway_name.keys(),
        'name': pathway_name.values()
    })

    gene_pathway_mapping = get_gene_pathway_mapping(
        organism=species
    )
    gene_pathway_mapping_df = pd.DataFrame({
        'ind': [path0[1] for path0 in gene_pathway_mapping],
        'values': [path0[0] for path0 in gene_pathway_mapping],
    })
    gene_pathway_mapping_df['values'] = gene_pathway_mapping_df['values'].str.replace(f'{species}:', '')
    gene_pathway_mapping_df['ind'] = gene_pathway_mapping_df['ind'].str.replace('path:', '')

    return pathway_name_df.dropna().drop_duplicates(), gene_pathway_mapping_df.dropna().drop_duplicates()

