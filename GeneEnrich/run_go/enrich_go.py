from GeneEnrich.run_go.go_utilities import get_go_pathway_mapping, get_go_gene_pathway_mapping
from GeneEnrich.enricher import enricher_internal
from GeneEnrich.preprocess import get_gene_info
import pandas as pd
import sys


def enrichgo(
    gene,
    organism="hsa",
    pvalueCutoff=0.05,
    pAdjustMethod="fdr_bh",
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff=0.2,
):

    res= pd.DataFrame()
    for onto in ['BP', 'MF', 'CC']:
        pathway_name_df, gene_pathway_mapping_df = get_data_from_GO_db(organism, onto)
        res_tmp = enricher_internal(gene,
                                pvalueCutoff=pvalueCutoff,
                                pAdjustMethod=pAdjustMethod,
                                minGSSize=minGSSize,
                                maxGSSize=maxGSSize,
                                qvalueCutoff=qvalueCutoff,
                                pathway_name_df=pathway_name_df,
                                gene_pathway_mapping_df=gene_pathway_mapping_df)
        if res_tmp is not None and not res_tmp.empty:
            res_tmp['ONTOLOGY'] = onto
            res = pd.concat([res, res_tmp])
        else:
            pass

    if len(res) != 0:
        gene2convert = {geneid: genesymbol for genesymbol, geneid in get_gene_info().items()}
        res['geneID'] = res['geneID'].apply(lambda x: '/'.join([gene2convert.get(i, i) for i in x.split('/')]))
        res = res.sort_values(by='pvalue')
    else:
        sys.exit('No enrichment analysis result can be found, check input file')

    return res


def get_data_from_GO_db(species: str, ont: str):

    pathway_name = get_go_pathway_mapping(
        organism=species
    )
    pathway_name_df = pd.DataFrame({
        'path': pathway_name.keys(),
        'name': pathway_name.values()
    })

    gene_pathway_mapping = get_go_gene_pathway_mapping(
        organism=species,
        ont=ont
    )
    gene_pathway_mapping_df = pd.DataFrame({
        'ind': [path0[1] for path0 in gene_pathway_mapping],
        'values': [path0[0] for path0 in gene_pathway_mapping],
    })
    gene_pathway_mapping_df['values'] = gene_pathway_mapping_df['values'].str.replace(f'{species}:', '')
    gene_pathway_mapping_df['ind'] = gene_pathway_mapping_df['ind'].str.replace('path:', '')

    return pathway_name_df, gene_pathway_mapping_df
