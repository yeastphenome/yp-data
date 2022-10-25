import sys
import re
import numpy as np
import pandas as pd

from scipy import stats

from yp_translate import *

path_to_genes = '../../Utils/yp_sgd_features.txt'
path_to_consensus_tested = '../../Utils/yp_consensus_tested_20200901.txt'
path_to_sgd_phenotypes = '../../Utils/yp_sgd_phenotypes.txt'


def normalize_phenotypic_scores(df, has_tested=False):
    
    if not has_tested:
        
        genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
        genes = genes.reset_index().set_index('systematic_name', drop=False)
        
        yp_orfs = pd.read_csv(path_to_consensus_tested, header=None)
        yp_orfs = yp_orfs[1].values
        consensus_tested = [tuple(genes.loc[orf,['id','systematic_name']]) for orf in yp_orfs]
        
        df_index = [tuple(x) for x in df.index]
        
        consensus_tested = list(set(consensus_tested + df_index))
        
        df = df.reindex(index=consensus_tested, fill_value=0)
        
    df_norm = z_transform_mode(df)
    
    return df_norm


def z_transform_mode(data):

    """
    Transform data into Z-like scores relative to the mode
    :param data: a
    :return:
    """

    data_z = data.copy()

    for i in data.columns:
        md = mode_kde(data[i])
        sigma = np.sqrt(np.nanmean((data[i] - md) ** 2))
        data_z[i] = (data[i] - md) / sigma

    return data_z


def mode_kde(data):
    bw = 0.25

    data = data[pd.notnull(data)]
    data = data.astype(float)

    if np.any(~np.isnan(data)):

        kde = stats.gaussian_kde(data)
        kde.set_bandwidth(bw)
        x = np.linspace(min(data), max(data), 100)

        md = x[np.argmax(kde(x))]

    else:

        md = np.nan

    return md


def clean_orf(lst):

    LST = [l.upper() for l in lst]
    LST = [re.sub("[^a-zA-Z0-9\-]", "", L) for L in LST]

    return LST


def clean_genename(lst):

    LST = [l.upper() for l in lst]
    LST = [re.sub("[^a-zA-Z0-9\',\-]", "", L) for L in LST]

    # Fix some common genename issues in yeast
    for idx, L in enumerate(LST):
        L = re.sub('ADE5-*7', 'ADE5,7', L)
        L = re.sub('ARG5-*6', 'ARG5,6', L)
        L = re.sub('DUR1-*2', 'DUR1,2', L)
        L = re.sub('MF-*ALPHA-*1', 'MF(ALPHA)1', L)
        L = re.sub('MF-*ALPHA-*2', 'MF(ALPHA)2', L)
        L = re.sub('AI5-*ALPHA', 'AI5_ALPHA', L)
        L = re.sub('AI5-*BETA', 'AI5_BETA', L)

        # See if correcting for a missing hyphens makes it look like an ORF (e.g., YAL064CA)
        if not looks_like_orf(L) and len(L) == 8:
            L0 = L[:7] + '-' + L[7]
            if looks_like_orf(L0):
                L = L0

        LST[idx] = L

    return LST


def looks_like(qq, patterns, case_sensitive=True):

    # Transform the input into a list to facilitate processing
    if isinstance(qq, pd.Series):
        qq_list = qq.tolist()
    elif isinstance(qq, str):
        qq_list = [qq]
    elif isinstance(qq, float) or isinstance(qq, int) or isinstance(qq, np.int64):
        return False
    else:
        qq_list = qq

    if not case_sensitive:
        qq_list = [q.upper() for q in qq_list]

    ww = [False] * len(qq_list)
    for pattern in patterns:
        regex = re.compile(pattern)
        pattern_match = [bool(regex.match(q)) if pd.notnull(q) else False for q in qq_list]
        ww = [a or b for a, b in zip(ww, pattern_match)]

    # Transform the translation into the type of the input
    if isinstance(qq, pd.Series):
        ww = pd.Series(ww, index=qq.index)
    elif isinstance(qq, list):
        ww = np.array(ww)
    elif isinstance(qq, str) or isinstance(qq, float) or isinstance(qq, int):
        ww = ww[0]

    return ww


def looks_like_orf(qq, case_sensitive=True):
    patterns = ['^Y[A-P][RL][0-9]{3}[CW](-[A-H])*$',
                'Q[0-9]{4}$']
    return looks_like(qq, patterns, case_sensitive)


def is_essential(lst):
   
    genes = pd.read_csv(path_to_genes, sep='\t', index_col='primary_sgdid')

    # Note: na_filter necessary to prevent pandas recognizing the "null" mutation as a missing value
    phenotypes = pd.read_csv(path_to_sgd_phenotypes, delimiter='\t', header=None, na_filter=False)
    phenotypes.index = phenotypes[3]

    phenotypes = phenotypes.loc[
        (phenotypes[6] == 'null') &
        (phenotypes[8] == 'S288C') &
        (phenotypes[5] == 'systematic mutation set') &
        (phenotypes[9] == 'inviable')
    ]

    phenotypes['orf'] = genes['systematic_name']
    lst = pd.Series(lst)

    return lst.isin(phenotypes['orf'])
