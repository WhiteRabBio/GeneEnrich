import numpy as np
import pandas as pd
import scipy.stats


def enricher_internal(
    gene,
    pathway_name_df,
    gene_pathway_mapping_df,
    pvalueCutoff=0.05,
    pAdjustMethod="fdr_bh",
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff=0.2,
):

    gene = [str(g)
            for g in np.unique(gene)]
    genePathwayMappingDF = gene_pathway_mapping_df[
        gene_pathway_mapping_df['values'].isin(gene)
    ]
    idx = genePathwayMappingDF.groupby(
        by='values'
    ).count()
    not_zero_idx = idx[idx > 0].index.tolist()
    qExtID2Path = genePathwayMappingDF[
        genePathwayMappingDF['values'].isin(not_zero_idx)
    ].drop_duplicates()
    if len(qExtID2Path) == 0:
        print("--> No gene can be mapped....")
        return None
    geneSets = qExtID2Path.copy()

    extID = gene_pathway_mapping_df['values'].unique().tolist()
    qTermID = qExtID2Path['ind'].unique().tolist()
    termID2ExtID = gene_pathway_mapping_df[
        gene_pathway_mapping_df['ind'].isin(qTermID)
    ]
    idx = get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
    termID2ExtID = termID2ExtID[
        termID2ExtID['ind'].isin(idx)
    ]
    qExtID2Path = qExtID2Path[
        qExtID2Path['ind'].isin(idx)
    ]
    qTermID = qExtID2Path['ind'].unique()

    k = qExtID2Path.groupby("ind").count().loc[qTermID]
    M = termID2ExtID.groupby("ind").count().loc[qTermID]
    N = [len(extID)] * len(M)
    n = [len(geneSets['values'].unique())] * len(M)
    df = pd.DataFrame({'numWdrawn': (k-1)['values'].values.tolist(),
                       'numW': M['values'].values.tolist(),
                       'numB': N,
                       'numDrawn': n}, index=qTermID)
    pvalues = scipy.stats.hypergeom.sf(
        df['numWdrawn'].tolist(),
        df['numB'].tolist(),
        df['numDrawn'].tolist(),
        df['numW'].tolist())
    pvalues = pd.DataFrame(pvalues, index=qTermID)

    GeneRatio = pd.DataFrame(['/'.join([str(j) for j in list(i)])
                              for i in zip(k['values'].values.tolist(), n)],
                             index=qTermID)
    BgRatio = pd.DataFrame(['/'.join([str(j) for j in list(i)])
                            for i in zip(M['values'].values.tolist(), N)],
                           index=qTermID)
    Over = pd.concat([pvalues, GeneRatio, BgRatio], axis=1).reset_index(drop=False)
    Over.columns = ['ID', 'pvalue', 'GeneRatio', 'BgRatio']

    import statsmodels.stats.multitest as sm
    adjusted_pvals = sm.multipletests(Over['pvalue'].tolist(),
                                      alpha=0.05,
                                      method=pAdjustMethod)[1]
    adjusted_pvals = pd.DataFrame(
        {'p.adjust': list(adjusted_pvals)},
        index=qTermID)

    geneID = pd.DataFrame(
        {'geneID': qExtID2Path.groupby('ind')['values'].apply(
            lambda x: '/'.join(x)).loc[qTermID].values.tolist()
         },
        index=qTermID)

    pathwayNameDf = pathway_name_df[pathway_name_df['path'].isin(qTermID)]
    pathwayNameDf = pd.DataFrame(
        {'Description': pathwayNameDf['name'].tolist()},
        index=pathwayNameDf['path']
    )

    count_df = pd.DataFrame(
        {'Count': [len(g.split("/")) for g in geneID['geneID']]},
        index=qTermID)

    qvalues = qvalue(Over['pvalue'].tolist(), lambda_=[0.05], pi0_method='bootstrap', pi0=None)
    qvalue_df = pd.DataFrame(
        {'qvalue': list(qvalues)},
        index=qTermID)

    Over.index = qTermID
    Over = pd.concat([Over, geneID, adjusted_pvals, pathwayNameDf, count_df, qvalue_df], axis=1)

    Over = Over[
        (Over['pvalue'] < pvalueCutoff)
        &
        (Over['p.adjust'] < pvalueCutoff)
        &
        (Over['qvalue'] < qvalueCutoff)
        ]
    Over = Over[['ID', 'Description', 'GeneRatio', 'BgRatio', 'pvalue', 'p.adjust', 'qvalue', 'geneID', 'Count']]
    Over = Over.sort_values(by='pvalue')

    return Over


def pi0est(
    p,
    pi0_method='bootstrap',
    lambda_=np.arange(0.05, 1, 0.05),
):

    p = [p_tmp for p_tmp in p if not np.isnan(p_tmp)]
    m = len(p)
    if len(lambda_) > 1:
        lambda_ = np.sort(lambda_)
    else:
        lambda_ = np.array(lambda_)

    ll = len(lambda_)
    if (np.min(p) < 0) | (np.max(p) > 1):
        print('Error: p-values not in valid range `[0, 1]`.')
        quit()
    if (ll > 1) and (ll < 4):
        print('Error: If length of lambda greater than 1, you need at least 4 values.')
        quit()
    if (np.min(lambda_) < 0) | (np.max(lambda_) >= 1):
        print('Error: Lambda must be within `[0, 1)`.')
        quit()
    if np.max(p) < np.max(lambda_):
        print('Error: maximum p-value is smaller than lambda range. Change the range of lambda.')
        quit()

    if ll == 1:
        pi0 = np.mean(
            np.array(p) > lambda_[0]
        ) / (
                1 - lambda_[0]
        )
        pi0_lambda_ = pi0
        pi0 = np.min([pi0, 1])
        pi0smooth = None
    else:
        ind = np.arange(len(lambda_)-1, -1, -1)
        search_sort = np.searchsorted(
            lambda_, np.array(p), side='right'
        ) - 1
        pi0 = np.cumsum(
            np.bincount(
                [s for s in search_sort if s > 0]
            )[ind]
        ) / (
                len(p) * (1 - lambda_)[ind]
        )
        pi0 = pi0[ind]
        pi0_lambda_ = pi0
        if pi0_method == 'bootstrap':
            minpi0 = np.quantile(pi0_lambda_, .10)
            w = np.array([np.sum(p >= l) for l in lambda_])
            mse = (w / (m ** 2 * (1 - lambda_) ** 2)) * (1 - w / m) + (pi0 - minpi0) ** 2
            pi0 = np.min([pi0[mse == np.min(mse)], 1])[0]
            pi0smooth = None
        else:
            print('Error: pi0_method only support `bootstrap`.')

    if pi0 <= 0:
        print('Warning: The estimated pi0 <= 0. '
              'Setting the pi0 estimate to be 1. '
              'Check that you have valid p-values or use a different range of lambda.')
        pi0 = 1
        pi0_lambda_ = 1
        pi0smooth = 0
        lambda_ = 0

    return pi0, pi0_lambda_, lambda_, pi0smooth


def qvalue(p, pi0, **kwargs):

    p = [p_tmp for p_tmp in p if not np.isnan(p_tmp)]
    if (np.min(p) < 0) | (np.max(p) > 1):
        print('Error: p-values not in valid range `[0, 1]`.')
        quit()

    pi0s = {}
    if not pi0:
        pi0s['pi0'], pi0s['pi0_lambda_'], pi0s['lambda_'], pi0s['pi0smooth'] = pi0est(p, **kwargs)
    else:
        if (pi0 > 0) and (pi0 <= 1):
            pi0s['pi0'] = pi0
        else:
            print('Error, pi0 is not `(0, 1]`')
            quit()

    m = len(p)
    i = np.arange(m, 0, -1)
    o = np.argsort(p)[::-1]
    ro = np.argsort(o)
    qvals = pi0s['pi0'] * np.minimum(1, np.minimum.accumulate(np.array(p)[o] * m / i))[ro]

    return qvals


def get_geneSet_index(
    geneSets,
    minGSSize=10,
    maxGSSize=500,
):
    geneSet_size = geneSets.groupby(
        "ind"
    ).count()
    geneSets = geneSet_size[
        (geneSet_size['values'] >= minGSSize)
        &
        (geneSet_size['values'] <= maxGSSize)
        ]
    return geneSets.index


def organism_mapper(species):

    if species == "human":
        organism = "hsa"
    elif species == "mouse":
        organism = "mmu"
    else:
        raise ValueError(f'Only human (`organism`="hsa") or mouse (`organism` = "mmu") pathways are implemented')

    return organism
