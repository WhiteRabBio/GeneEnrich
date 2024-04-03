import pandas as pd
import os
import pickle
import numpy as np
from GeneEnrich.run_kegg.kegg_utilities import save, load

dir_, _ = os.path.split(__file__)
database_dir = os.path.abspath(os.path.join(dir_, '..', 'database'))
download_dir = os.path.abspath(os.path.join(dir_, '..', 'extract_db_from_R'))


def update_go_db(organism='hsa', force=False, version='R'):
    if version == 'python':
        pass
    else:
        for ont in ['BP', 'MF', 'CC']:
            cmd = f'Rscript {download_dir}/extract_godb_from_R.r ' \
                  f'--species {organism} ' \
                  f'--output {download_dir}/GO_db  ' \
                  f'--ont {ont}'
            os.system(cmd)

            annotation_df = pd.read_table(f'{download_dir}/GO_db/TermID_GO_{ont}_{organism}_df.xls')
            annotation_tuples = annotation_df.apply(
                lambda annotation_df: (annotation_df['geneid'], annotation_df['pathid']), axis=1
            ).tolist()
            save(tuples=annotation_tuples, file=f'{database_dir}/go/{organism}_{ont}_annotations.txt')

        pathway_df = pd.read_table(f'{download_dir}/GO_db/PathwayName_GO_{organism}_df.xls')
        pathway_df = pathway_df[np.isin(pathway_df['pathid'], 'all', invert=True)]
        pathway_tuples = pathway_df.apply(
            lambda pathway_df: (pathway_df['pathid'], pathway_df['pathname']), axis=1
        ).tolist()
        save(tuples=pathway_tuples, file=f'{database_dir}/go/{organism}_pathway.txt')

def get_go_pathway_mapping(organism='hsa'):
    pathway_tuples = load(f'{database_dir}/go/{organism}_pathway.txt')
    return {pathway_id: name for pathway_id, name in pathway_tuples}


def get_go_gene_pathway_mapping(organism='hsa', ont='BP'):
    annotation_tuples = load(f'{database_dir}/go/{organism}_{ont}_annotations.txt')
    return [(g, p) for g, p in annotation_tuples]
