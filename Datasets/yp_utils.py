import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *
from Math.normalize import z_transform_mode

path_to_genes = '../../Private-Utils/datasets_gene2.txt'
path_to_consensus_tested = '../../Private-Utils/yp_2020-09-01_orfs.txt'


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