from collections import defaultdict
from logging import getLogger, basicConfig

import numpy as np

from calour.dsfdr import dsfdr

logger = getLogger(__name__)
basicConfig(format='%(levelname)s:%(message)s')


def sample_enrichment(exp_features, group1_data, group2_data, db_annotations, db_sequence_annotations,
                      term_type='term', ignore_exp=None, min_appearances=3, fdr_method='dsfdr', score_method='all_mean', freq_weight='log', alpha=0.1):
    '''Get the list of enriched terms for all bacteria between two groups.

    It is equivalent to multiplying the (freq_weight transformed) feature X sample matrix by the database derived term X feature matrix
    (where the numbers are how strong is the term associated with the feature based on database annotations using score_method).
    A differntial abundance test (using dsFDR multiple hypothesis correction) is then performed on the resulting sample X term matrix.

    Parameters
    ----------
    exp_features: list of str
        list of sequences (same order as columns in the two data arrays)
    group1_data: np.array
        number of reads of group1 samples (rows) X features (columns)
    group2_data: np.array
        number of reads of group2 samples (rows) X features (columns)
    term_type : str or None (optional)
        The type of annotation data to test for enrichment
        options are:
        'term' - ontology terms associated with each feature.
        'parentterm' - ontology terms including parent terms associated with each feature.
        'annotation' - the full annotation strings associated with each feature
        'combined' - combine 'term' and 'annotation'
    ignore_exp: list of int or None
        List of experiments to ignore in the analysis
        None (default) to use annotations from all experiments including the current one
    min_appearances : int (optional)
        The minimal number of times a term appears in order to include in output list.
    fdr_method : str (optional)
        The FDR method used for detecting enriched terms (permutation test). options are:
        'dsfdr' (default): use the discrete FDR correction
        'bhfdr': use Benjamini Hochbert FDR correction
    score_method : str (optional)
        The score method used for calculating the term score. Options are:
        'all_mean' (default): mean over each experiment of all annotations containing the term
        'sum' : sum of all annotations (experiment not taken into account)
        'card_mean': use a null model keeping the number of annotations per each bacteria
    freq_weight : str (optional)
        How to incorporate the frequency of each feature (in a sample) into the term count for the sample. options:
            'linear' : use the frequency
            'log' : use the log2 of the frequency
            'binary' : use the presence/absence of the feature
    alpha : float (optional)
        the FDR level desired (0.1 means up to 10% of results can be due to null hypothesis)

    Returns
    -------
    pandas.DataFrame with  info about significantly enriched terms. columns:
        feature : str the feature
        pval : the p-value for the enrichment (float)
        odif : the effect size (float)
        observed : the number of observations of this term in group1 (int)
        expected : the expected (based on all features) number of observations of this term in group1 (float)
        frac_group1 : fraction of total terms in group 1 which are the specific term (float)
        frac_group2 : fraction of total terms in group 2 which are the specific term (float)
        num_group1 : number of total terms in group 1 which are the specific term (float)
        num_group2 : number of total terms in group 2 which are the specific term (float)
        description : the term (str)
numpy.Array where rows are features (ordered like the dataframe), columns are features and value is score
        for term in feature
    pandas.DataFrame with info about the features used. columns:
        group: int the group (1/2) to which the feature belongs
        sequence: str
    '''

    if term_type == 'term':
        feature_terms = _get_all_term_counts(exp_features, db_sequence_annotations, db_annotations, ignore_exp=ignore_exp, score_method=score_method)
    elif term_type == 'parentterm':
        pass
    elif term_type == 'annotation':
        feature_terms = _get_all_annotation_string_counts(exp_features, db_annotations, db_sequence_annotations, ignore_exp=ignore_exp)
    elif term_type == 'combined':
        feature_terms = _get_all_term_counts(exp_features, db_sequence_annotations, db_annotations, ignore_exp=ignore_exp, score_method=score_method)
        feature_annotations = _get_all_annotation_string_counts(exp_features, db_annotations, db_sequence_annotations, ignore_exp=ignore_exp)
        for cfeature,cvals in feature_annotations.items():
            if cfeature not in feature_terms:
                feature_terms[cfeature]=[]
            feature_terms[cfeature].extend(cvals)
    else:
        raise ValueError('term_type %s not supported for dbbact. possible values are: "term", "parentterm", "annotation", "combined"')

    # get all terms. store the index position for each term
    terms = {}
    term_features = defaultdict(int)
    cpos = 0
    for cfeature, ctermlist in feature_terms.items():
        for cterm, ccount in ctermlist:
            if cterm not in terms:
                terms[cterm] = cpos
                cpos += 1
            term_features[cterm] += 1

    logger.debug('found %d terms for all features' % len(terms))

    # create the term x sample array
    num_samples = group1_data.shape[0] + group2_data.shape[0]
    fs_array = np.zeros([num_samples, len(terms)])

    data = np.vcat([group1_data, group2_data])
    labels = np.hcat([np.zeros(group1_data.shape[0]), np.ones(group1_data.shape[0])])

    # how to weight the frequency of each bacteria
    if freq_weight == 'log':
        data[data < 1] = 1
        data = np.log2(data)
    elif freq_weight == 'binary':
        data = data>0
    elif freq_weight == 'linear':
        pass
    else:
        raise ValueError('unknown freq_weight option %s. Can use ["log", "binary","linear"].')

    # iterate over all features and add to all terms associated with the feature
    for idx, cfeature in enumerate(exp_features):
        fterms = feature_terms[cfeature]
        for cterm,cval in fterms:
            fs_array[:, terms[cterm]] += cval * data[:, idx]

    # get the differentially abundant terms between the two sample groups
    keep, odif, pvals, qvals = dsfdr(fs_array.transpose(), labels, transform_type=None, alpha=alpha, fdr_method=fdr_method)


def _get_all_term_counts(features, feature_annotations, annotations, ignore_exp=None, score_method='all_mean'):
    '''
    Get the term scores for each feature in features

    Parameters
    ----------
    features: list of str
        list of feature sequences to use for the score calculation
        feature_annotations: dict {str: list of int}
            per feature (key) list of annotation ids
        annotations : dict of {int: dict}
            key is annotationID, dict is the annotation details (see XXX)
        ignore_exp : list of int or None (optional)
            None (default) to use all experiments
            list of experimentIDs to not include annotations from these experiments in the score
        score_method: str (optional)
            The method to use for score calculation:
            'all_mean' : score is the mean (per experiment) of the scoring for the term out of all annotations that have the term
            'sum' : score is sum of all annotations where the term appears (experiment is ignored)

    Returns
    -------
    dict of {feature:str, list of tuples of (term:str, score:float)}
        key is feature, value is list of (term, score)
    '''
    # set up the list of annotations per experiment if needed
    if score_method == 'all_mean':
        exp_annotations =_get_exp_annotations(annotations)
    elif score_method == 'sum':
        exp_annotations = None
    else:
        logger.warn('score_method %s not supported' % score_method)
        return None

    # print('annotations:')
    # print(feature_annotations)

    # calculate the per-feature score for each term
    feature_terms = {}
    for cfeature in features:
        annotation_list = []
        for cannotation in feature_annotations[cfeature]:
            # test if we should ignore the annotation (if experiment is in ignore_exp)
            if ignore_exp is not None:
                if annotations[cannotation]['expid'] in ignore_exp:
                    continue
            annotation_list.append(annotations[cannotation])
        # print('---------------------')
        # print(cfeature)
        feature_terms[cfeature] = get_annotation_term_counts(annotation_list, exp_annotations=exp_annotations, score_method=score_method)
        # print(feature_terms)
    return feature_terms


def _get_all_annotation_string_counts(features, db_annotations, db_sequence_annotations, ignore_exp=None):
    feature_annotations = {}
    for cseq, annotations_list in db_sequence_annotations.items():
        if cseq not in features:
            continue
        newdesc = []
        for cannotation in annotations_list:
            if ignore_exp is not None:
                if db_annotations[cannotation]['expid'] in ignore_exp:
                    continue
            cdesc = get_annotation_string(db_annotations[cannotation])
            newdesc.append((cdesc, 1))
        feature_annotations[cseq] = newdesc
    return feature_annotations


def _get_exp_annotations(annotations):
    '''
    Get all annotations for each experiment in annotations

    Parameters
    ----------
    annotations : dict of {annotationid:int : annotation:dict}


    Returns
    -------
    dict of {expid:int : dict of {term:str : total:int}}
    '''
    # exp_annotations = {}
    exp_annotations = defaultdict(lambda: defaultdict(int))
    for cannotation in annotations.values():
        cexpid = cannotation['expid']
        if cannotation['annotationtype'] == 'contamination':
            exp_annotations[cexpid]['contamination'] += 1
            continue
        for cdetail in cannotation['details']:
            ctype = cdetail[0]
            cterm = cdetail[1]
            if ctype == 'low':
                cterm = '-' + cterm
            exp_annotations[cexpid][cterm]+=1
    return exp_annotations


def get_annotation_term_counts(annotations, exp_annotations=None, score_method='all_mean'):
    '''Get the annotation type corrected count for all terms in annotations

    Parameters
    ----------
    annotations : list of dict
        list of annotations where the feature is present
    dict of {expid:int : dict of {term:str : total:int}}
        from self._get_exp_annotations()
    score_method: str (optional)
            The method to use for score calculation:
            'all_mean' : score is the mean (per experiment) of the scoring for the term out of all annotations that have the term
            'sum' : score is sum of all annotations where the term appears (experiment is ignored)

    Returns
    -------
        list of tuples (term, count)
    '''
    term_count = defaultdict(int)

    for cannotation in annotations:
        # print('****%d' % cannotation['annotationid'])
        details = cannotation['details']
        annotation_type = cannotation['annotationtype']
        if annotation_type == 'common':
            cscore = 1
        elif annotation_type == 'dominant':
            cscore = 1
        elif annotation_type == 'other':
            cscore = 0.5
        elif annotation_type == 'contamination':
            cscore = 1
            details = [('all', 'contamination')]
        elif annotation_type == 'other':
            cscore = 0
        elif annotation_type == 'diffexp':
            cscore = None
        else:
            logger.warn('unknown annotation type %s encountered (annotationid %d). skipped' % (annotation_type, cannotation['annotationid']))
            continue
        for cdetail in details:
            ctype = cdetail[0]
            cterm = cdetail[1]
            if ctype == 'all':
                cscore = 1
            elif ctype == 'high':
                cscore = 2
            elif ctype == 'low':
                # if low - change the term to "-term" - i.e. lower in term
                cterm = '-' + cterm
                cscore = 1
            else:
                logger.warn('unknown detail type %s encountered for annotation %d' % (ctype, cannotation['annotationid']))

            if score_method == 'all_mean':
                scale_factor = exp_annotations[cannotation['expid']][cterm]
                if scale_factor == 0:
                    print('term %s, exp %d, annotationid %d, score %d, scale %d' % (cterm, cannotation['expid'], cannotation['annotationid'], cscore, scale_factor))
                    scale_factor = 1
            else:
                scale_factor = 1

            # fix for differential abundance
            term_count[cterm] += cscore / scale_factor

    res = []
    for k, v in term_count.items():
        # flip and add '-' to term if negative
        if v < 0:
            k = '-' + k
            v = -v
        res.append((k, v))
    return res


def get_annotation_string(cann):
    '''Get nice string summaries of annotation

    Parameters
    ----------
    cann : dict
        items of the output of get_seq_annotations()

    Returns
    -------
    desc : str
        a short summary of the annotation
    '''
    cdesc = ''
    if cann['description']:
        cdesc += cann['description'] + ' ('
    if cann['annotationtype'] == 'diffexp':
        chigh = []
        clow = []
        call = []
        for cdet in cann['details']:
            if cdet[0] == 'all':
                call.append(cdet[1])
                continue
            if cdet[0] == 'low':
                clow.append(cdet[1])
                continue
            if cdet[0] == 'high':
                chigh.append(cdet[1])
                continue
        cdesc += ' high in '
        for cval in chigh:
            cdesc += cval + ' '
        cdesc += ' compared to '
        for cval in clow:
            cdesc += cval + ' '
        cdesc += ' in '
        for cval in call:
            cdesc += cval + ' '
    elif cann['annotationtype'] == 'isa':
        cdesc += ' is a '
        for cdet in cann['details']:
            cdesc += 'cdet,'
    elif cann['annotationtype'] == 'contamination':
        cdesc += 'contamination'
    else:
        cdesc += cann['annotationtype'] + ' '
        for cdet in cann['details']:
            cdesc = cdesc + ' ' + cdet[1] + ','
    return cdesc

