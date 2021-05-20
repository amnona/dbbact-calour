from collections import defaultdict
from logging import getLogger

import numpy as np
import pandas as pd
from .db_access import _get_common_seqs_for_terms

import calour as ca

logger = getLogger(__name__)


def _get_sequences(db, annotation_ids, min_appearance=2):
    '''Get all sequences appearing in at least min_appearance of the annotation ids

    Parameters
    ----------
    db: DBAcess
        the dbaccess class to interface with dbBact
    anootation_ids: list of int
        the annotation ids to query
    min_appearance: int, optional
        the minimal number of annotations the sequence appears in, in order to include the sequence in the output

    Returns
    -------
    seqs: list of str
        the sequences appearing in the annotations
    seq_list: list of str
        same as seqs, but each sequence appears X times, where X is the number of annotations it appears in
    '''
    seqs_ids = defaultdict(int)
    seqs_list = []
    for cid in annotation_ids:
        res = db.get_annotation_sequences(cid)
        for cseq in res.json()['seqids']:
            seqs_ids[cseq] += 1
    # keep only sequences appearing in enough annotations
    ok_seqs = [ck for ck, cv in seqs_ids.items() if cv >= min_appearance]
    res = db.get_seq_id_info(ok_seqs)
    seqs = [x['seq'].upper() for x in res]
    # create the inflated list of sequences (each sequence appears according to the number of annotations it is found in)
    for idx, x in enumerate(res):
        seqs_list.extend([x['seq'].upper()] * seqs_ids[ok_seqs[idx]])
    logger.debug('found %d sequences associated with annotations' % len(seqs))
    return seqs, seqs_list


def _background_enrich(db, terms, exp, ignore_exp=True, min_appearance=2, include_shared=True, alpha=0.1):
    '''Find dbbact term enrichment comparing experiment sequences (COMMON bacteria) to all sequences in dbbact associated with the terms

    Parameters
    ----------
    db: DBAcess
        the database interface (for accessing dbBact)
    terms: list of str
        the terms to get the dbbact sequences for (AND). Retrieves only COMMON annotations
    exp: calour.AmpliconExperiment
        COMMON sequences (prev>0.5) are compared to the dbbact background sequences
    ignore_exp: None or bool or list of int, optional
        int: dbBact Experiment IDs to ignore when fetching sequences associated with the terms
        True: ignore annotations from the current experiment
        None: do not ignore any annotations
    min_appearance: int, optional
        keep only sequences associated with the terms in at least min_appearance annotations
    include_shared: bool, optional
        True: keep sequences present in both background and experiment (add them to both groups)
        False: ignore sequences present in both background and experiment (remove from both groups)

    Returns
    -------
    Pandas.DataFrame of enriched terms for common terms in current experiment compared to background dbBact experiments.
    Positive is higher in dbbact background seqs, negative is higher in the experiment seqs
    '''
    logger.info('getting term sequences')
    anno = _get_common_seqs_for_terms(db, terms)
    seqs, seqs_all = _get_sequences(anno, min_appearance=min_appearance)
    logger.info('preparing data')

    # remove shared bacteria between query and background if include_shared is False
    # NOT IMPLEMENTED
    if not include_shared:
        raise ValueError('not implemented yet')
        shared_seqs = set(seqs).intersection(set(exp.feature_metadata.index.values))
        print('removing %d shared sequences' % len(shared_seqs))
        seqs = list(set(seqs).difference(shared_seqs))
        print('%d bg seqs remaining' % len(seqs))
        exp = exp.filter_ids(list(shared_seqs), negate=True)

    exp = exp.filter_prevalence(0.5)
    logger.info('%d COMMON experiment sequences remaining' % len(exp.feature_metadata))
    features = pd.DataFrame(seqs, columns=['_feature_id'])
    features = features.set_index('_feature_id', drop=False)
    samples = exp.sample_metadata.copy()
    aa = ca.AmpliconExperiment(data=np.ones([len(samples), len(seqs)]), sample_metadata=samples, feature_metadata=features).normalize()
    bb = exp.join_experiments(aa, field='pita', prefixes=['e1', 'e2'])
    logger.info('getting annotations for %d sequences' % len(bb.feature_metadata))
    bb = bb.add_terms_to_features('dbbact')
    bb.info = exp.info.copy()
    logger.debug('Calculating diff. abundance')
    cc = db.enrichment(bb, seqs_all, bg_features=exp.feature_metadata.index.values, ignore_exp=ignore_exp, alpha=alpha)
    return cc
