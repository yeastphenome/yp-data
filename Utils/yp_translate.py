import sys
import string
import random
import pickle
import re
import numpy as np
import pandas as pd

from itertools import chain

from os.path import join
from pandas.errors import MergeError


def translate_sc(qq, frm='', to='orfs', untranslated='keep'):

    data_folder = '../../Utils/'

    try:
        with open(join(data_folder, 'entrez3.p'), 'rb') as handle:
            entrez = pickle.load(handle)
    except AttributeError:
        # Works with pandas 1.0.3
        entrez = pd.read_pickle(join(data_folder, 'entrez3.p'))

    return translate(qq, translation_map=entrez, frm=frm, to=to, untranslated=untranslated)


def translate(qq, translation_map, frm='', to='ensembl', untranslated='keep', untranslated_return=False):

    # --- PREPARATION ---

    if isinstance(qq, str):
        qq_original = qq
        qq = pd.Series(qq)
    elif ~isinstance(qq, pd.Series):
        qq_original = qq.copy()
        qq = pd.Series(qq)
    else:
        qq_original = qq.copy()

    to = to.lower()
    ref_id = 'entrezid'

    entrez_all2id = translation_map['entrez_all2id']
    entrez_id2all = translation_map['entrez_id2all']

    switcher = {
        'ensembl': 'EnsemblID',
        'symbol': 'Symbol',
        'synonym': 'Synonyms',
        'uniprot': 'UniProt',
        'entrez': 'EntrezID',
        'orf': 'LocusTag',
        'db_object_id': 'DB_Object_ID',
        'omim': 'OMIM',
    }

    try:
        to = switcher[to]
    except KeyError:
        print('Unknown destination format. The available options are: %s.' % ', '.join(list(switcher.keys())))
        return

    if frm:
        try:
            frm = switcher[frm]
        except KeyError:
            print('Unknown source format. The available options are: %s.' % ', '.join(list(switcher.keys())))
            return

    # --- SETUP ---

    qq = pd.DataFrame(data={'label': qq, 'label_type': frm}, index=qq.index)
    qq.index.name = 'index_input'
    
    # If EnsemblID are involved, drop the version number (.XX)
    qq['label'] = qq['label'].apply(lambda x: x.split('.')[0] if looks_like_ensemblid(x) else x)

    qq = qq.reset_index()
    
    if frm:

        # Make sure that the frm label_type exists in the translation map (unlikely but possible scenario)
        all_label_types = entrez_all2id.index.get_level_values(level='label_type').unique().values
        if frm not in all_label_types:
            print('%s is not a valid source label.' % frm)

        # Switch EntrezIDs to strings, if the user hasn't done so already
        if frm == 'EntrezID':
            qq['label'] = qq['label'].astype(str)

        qq.set_index(['label', 'label_type'], drop=True, inplace=True)

        # Map everything to the reference ID
        translation = qq.reset_index().merge(entrez_all2id.reset_index(), on=['label', 'label_type'], how='left')

    else:

        qq.set_index('label', drop=True, inplace=True)
        entrez_all2id = entrez_all2id.reset_index()
        grp = entrez_all2id.groupby('label')
        t1 = grp['label_type'].apply(list)
        t2 = grp['entrezid'].apply(list).apply(lambda x: list(chain(*x)))
        entrez_all2id = pd.concat([t1, t2], axis=1)

        # Map everything to the reference ID
        try:
            translation = qq.reset_index().merge(entrez_all2id.reset_index(), on=['label'], how='left',
                                                 validate='many_to_one')
            translation.rename(columns={'label_type_y': 'label_type'}, inplace=True)
        except MergeError:
            print('The translation map contains more than one entry for items on your list. Specify frm.')
            return
        
    # --- TRANSLATION ---
    
    # Map everything to the reference ID

    # Choose one entrezid if multiple are available
    # Logic:
    # - if the multiple entrez_ids are due to an ambiguous source, use a priority list;
    # - if not, just pick the first one
    def source_priority(v):
        sort_dict = {'Symbol': 0, 'Synonyms': 1}
        if v in sort_dict.keys():
            p = sort_dict[v]
        else:
            p = 2
        return p

    def choose_entrezid(x):
        if not isinstance(x['entrezid'], list):    # single entrezid
            y = x['entrezid']
        else:
            if not isinstance(x['label_type'], list):   # multiple entrezids, but single source label_type
                y = x['entrezid'][0]
            else:    # multiple entrezids, multiple source label_types
                # Sort according to label_type, so that Symbol, Synonyms come first and everything else second
                sorted_entrezids = [x for _, x in sorted(zip(x['label_type'], x['entrezid']),
                                                         key=lambda pair: source_priority(pair[0]))]
                y = sorted_entrezids[0]
        return y

    translation[ref_id] = translation.apply(lambda x: choose_entrezid(x), axis=1)
    
    # Can't have NaNs in the ref_id column (needed for join) --> replace NaNs with random strings
    is_untranslated = pd.isnull(translation[ref_id])    # first, keep a note about which entries were null

    def mask_nan_w_random_string(x): return ''.join(np.random.choice(list(string.ascii_uppercase), size=10)) if pd.isnull(x) else x
    translation[ref_id] = translation[ref_id].apply(lambda x: mask_nan_w_random_string(x))    # second, replace the null entries with random strings

    if to == 'EntrezID':
        translation['label_to'] = translation[ref_id]
        translation['label_frm'] = translation['label']
    else:
        # Now, map reference IDs to the right "to" format
        translation['label_type'] = to
        translation = translation.merge(entrez_id2all.reset_index(), on=[ref_id, 'label_type'], how='left', suffixes=('_frm', '_to'))
        translation['label_to'] = translation['label_to'].apply(lambda x: x[0] if isinstance(x, list) else x)

    translation = translation.reset_index().set_index('index_input', drop=True)
    translation.index.name = 'index'
    
    # Deal with the untranslated items
    is_untranslated = is_untranslated | pd.isnull(translation['label_to'])
    if untranslated == 'keep':
        translation.loc[is_untranslated, 'label_to'] = translation.loc[is_untranslated, 'label_frm']

    if isinstance(qq_original, pd.Series):
        ww = translation['label_to']
    elif isinstance(qq_original, str):
        ww = translation['label_to'].values[0]
    else:
        ww = translation['label_to'].values

    if untranslated_return:
        return ww, is_untranslated
    else:
        return ww


def looks_like_ensemblid(qq):
    patterns = ['^ENSG[0-9]{11}']
    return looks_like(qq, patterns)


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


# translate(str(sys.argv[1]))

