from GeneEnrich.parse import Request
import pandas as pd
import pickle
import os

dir_, _ = os.path.split(__file__)
database_dir = os.path.abspath(os.path.join(dir_, '..', 'database'))


def get_gene_id_mapping(organism='hsa'):
    keggd = load(f'{database_dir}/kegg/{organism}_geneid.txt')
    Kgid = {}
    for line in keggd:
        Kgid[line[0].split(':')[1]] = line[1]
    return Kgid


def get_pathway_mapping(organism='hsa'):
    pathway_tuples = load(f'{database_dir}/kegg/{organism}_pathway.txt')
    return {pathway_id: name for pathway_id, name in pathway_tuples}


def get_gene_pathway_mapping(organism='hsa'):
    annotation_tuples = load(f'{database_dir}/kegg/{organism}_annotations.txt')
    return [(g, p) for g, p in annotation_tuples]


def update_kegg_db(organism='hsa', force=False, version='python'):
    keggd_tuples = Request('conv', organism, 'ncbi-geneid', force=force)
    save(tuples=keggd_tuples, file=f'{database_dir}/kegg/{organism}_geneid.txt')
    if version == 'python':
        pathway_tuples = Request('list', 'pathway', organism, force=force)
        annotation_tuples = Request('link', 'path', organism, force=force)
    else:
        download_dir = os.path.abspath(os.path.join(dir_, '..', 'extract_db_from_R'))
        cmd = f'Rscript {download_dir}/extract_keggdb_from_R.r --species {organism} --output {download_dir}/KEGG_db'
        os.system(cmd)
        pathway_df = pd.read_table(f'{download_dir}/KEGG_db/PathwayName_KEGG_{organism}_df.xls')
        pathway_tuples = pathway_df.apply(
            lambda pathway_df: (pathway_df['pathid'], pathway_df['pathname']), axis=1
        ).tolist()
        annotation_df = pd.read_table(f'{download_dir}/KEGG_db/TermID_KEGG_{organism}_df.xls')
        annotation_tuples = annotation_df.apply(
            lambda annotation_df: (annotation_df['geneid'], annotation_df['pathid']), axis=1
        ).tolist()
    save(tuples=pathway_tuples, file=f'{database_dir}/kegg/{organism}_pathway.txt')
    save(tuples=annotation_tuples, file=f'{database_dir}/kegg/{organism}_annotations.txt')


def save(tuples=None, file=None):
    with open(file, 'wb') as f:
        pickle.dump(tuples, f)


def load(file=None):
    with open(file, 'rb') as f:
        tuples = pickle.load(f)
    return tuples