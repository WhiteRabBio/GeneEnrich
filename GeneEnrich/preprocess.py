import pandas as pd
import numpy as np
import os
from GeneEnrich.run_kegg.kegg_utilities import get_gene_id_mapping, load
from GeneEnrich.consts import SCANPY_TITLE, SEURATV3_TITLE
from GeneEnrich.consts import DEG_P_VALUE_CUTOFF, DEG_MIN_PCT_CUTOFF, DEG_LOG_FC_CUTOFF


def get_gene_info():
	dir_, _ = os.path.split(__file__)
	database_dir = os.path.abspath(os.path.join(dir_, 'database'))
	symbol_tuple = load(f'{database_dir}/geneinfo.txt')
	return {genesymbol: str(geneid) for geneid, genesymbol in symbol_tuple}


def preprocess_deg(deg_path, type='deg', updown='up', organism='hsa', database='KEGG'):

	if type == 'scanpy_deg':
		df = pd.read_table(deg_path, sep="\t")
		assert set(SCANPY_TITLE).issubset(set(df.columns)) , \
            "Not a valid Scanpy DEG File."
		if updown == 'up':
			df = df[df.logfoldchanges >= DEG_LOG_FC_CUTOFF]
			df = df[df.pct_nz_group >= df.pct_nz_reference]
		else:
			df = df[df.logfoldchanges <= -DEG_LOG_FC_CUTOFF]
			df = df[df.pct_nz_group <= df.pct_nz_reference]
		df = df.query(
			f' ( pct_nz_group>={DEG_MIN_PCT_CUTOFF} | pct_nz_reference>={DEG_MIN_PCT_CUTOFF} ) & '
			f' ( pvals<={DEG_P_VALUE_CUTOFF} )'
		).sort_values(
			'logfoldchanges', ascending=False
		).reset_index(
			drop=True
		)
		df['old_names'] = df['names']
	elif type == 'seuratv3_deg':
		df = pd.read_table(deg_path, sep="\t")
		assert set(SEURATV3_TITLE).issubset(set(df.columns)), \
			"Not a valid Seurat v3 version DEG File."
		df = df.rename(columns=lambda x: x.replace('.', '_'))
		if updown == 'up':
			df = df[df.avg_logFC >= DEG_LOG_FC_CUTOFF]
		else:
			df = df[df.avg_logFC <= -DEG_LOG_FC_CUTOFF]
		df = df.query(
			f' ( pct_1>={DEG_MIN_PCT_CUTOFF} | pct_2>={DEG_MIN_PCT_CUTOFF} ) & '
			f' ( p_val<={DEG_P_VALUE_CUTOFF} )'
		).sort_values(
			'avg_logFC', ascending=False
		).reset_index(
			drop=True
		)
		df['old_names'] = df['gene_id']
	else:
		print('Not a valid-format input file. ')
		quit()

	gene2convert = get_gene_info()
	gene_id_mapping = get_gene_id_mapping(organism)
	id_lst = [gene2convert.get(g, np.nan) for g in df['old_names'].tolist()]
	if database == 'KEGG':
		id_lst = [gene_id_mapping[g].replace(f'{organism}:', '')
				  if g in gene_id_mapping.keys() else np.nan for g in id_lst]
	df['new_names'] = id_lst
	print(f"{'%.2f%%' % ((1 -len(df.dropna()) / len(df)) * 100)} skipped, not found in database")
	df = df.dropna()

	return list(df['new_names'])
