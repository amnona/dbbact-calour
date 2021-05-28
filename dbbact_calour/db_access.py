'''
db_access (:mod:`dbbact_calour.db_access`)
====================================

.. currentmodule:: dbbact_calour.dbbact

Functions
^^^^^^^^^
.. autosummary::
   :toctree: generated

   DBAccess
   DBAccess._get
   DBAccess._post
   DBAccess.get_seq_annotations
'''


import requests
from collections import defaultdict
from logging import getLogger
from pkg_resources import parse_version

import numpy as np
import pandas as pd

from calour.dsfdr import dsfdr
from .term_pairs import get_enrichment_score

from . import __version__

logger = getLogger(__name__)


class DBAccess():
    '''Handle all dbBact REST-API facing functions, at various levels (with user auhtentication) including:
    low level: _get(), _post()
    medium level: get/update/add annotations
    high level: enrichment analysis

    When creating the DBAcess class, it stores the rest api webserver address and the username/password, which are then used
    in all function calls from this class.
    '''
    def __init__(self, dburl='http://api.dbbact.org', username=None, password=None):
        '''Create the DBAccess class

        Parameters:
        -----------
        dburl: str, optional
            Web address of the dbBact REST API server
        username: str, optional
            Username or None for annonimous user
        password: str, optional
            None if username is None
        '''
        self.dburl = dburl

        self.username = username
        self.password = password

        try:
            # check compatibility with dbBact server supported versions
            res = self._get('stats/get_supported_version', {'client': 'dbbact_calour'})
            if res is None:
                logger.warn('Could not validate dbbact_calour_version')
                return
            if 'min_version' in res:
                if parse_version(res['min_version']) > parse_version(self.version()):
                    logger.warn('Current dbbact_calour version (%s) is not supported. please update at least to %s' % (self.version(), res['min_version']))
                    return
            if 'current_version' in res:
                if parse_version(res['current_version']) > parse_version(self.version()):
                    logger.warn('A new dbbact_calour version (%s) is available. Update recommended' % res['current_version'])
                    logger.warn('(Current dbbact_calour version: %s)' % self.version())
        except Exception as err:
            logger.error('Connection to dbBact failed. Error: %s' % err)
            return

    def version(self):
        '''Get the dbbact version

        Returns
        -------
        version: str
        '''
        return __version__

    def _post(self, api, rdata):
        '''POST a request to dbBact using authorization parameters

        Parameters
        ----------
        api : str
            the REST API address to post the request to (without the heading /)
        rdata : dict
            parameters to pass to the dbBact REST API

        Returns
        -------
        res : request
            the result of the request
        '''
        rdata['user'] = self.username
        rdata['pwd'] = self.password
        res = requests.post(self.dburl + '/' + api, json=rdata)
        if res.status_code != 200:
            logger.warn('REST error %s enountered when accessing dbBact %s' % (res.reason, api))
        return res

    def _get(self, api, rdata, param_json=True):
        '''GET a request to dbBact using authorization parameters

        Parameters
        ----------
        api : str
            the REST API address to post the request to (without the heading /)
        rdata : dict
            parameters to pass to the dbBact REST API
        param_json : bool (optional)
            True (default) to use rdata as json, False to use rdata as parameter

        Returns
        -------
        res : request
            the result of the request
        '''
        rdata['user'] = self.username
        rdata['pwd'] = self.password
        if param_json:
            res = requests.get(self.dburl + '/' + api, json=rdata)
        else:
            res = requests.get(self.dburl + '/' + api, params=rdata)

        if res.status_code != 200:
            logger.warn('REST error %s encountered when accessing dbBact %s' % (res.reason, api))
        return res

    def get_primers(self):
        '''Get a list of all the primers available in the database

        Returns
        -------
        list of dict of {'primerid': int
                    dbbact internal id of the primer region (i.e. 1 for v4, etc.)
                'name': str,
                    name of the primer region (i.e. 'v4', 'its1', etc.)
                'fprimer': str
                'rprimer: str
                    name of the forward and reverse primers for the region (i.e. 515f, etc.)}
                'fprimerseq': str
                    the forward primer concensus sequence
        '''
        res = self._get('sequences/get_primers', rdata={})
        return res.json()['primers']

    def get_seq_annotations(self, sequence, max_id=None):
        '''Get the annotations for a sequence

        Parameters
        ----------
        sequence : str
            The DNA sequence to get the annotations for
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id

        Returns
        -------
        annotations : list of list of (annotation dict,list of [Type,Value] of annotation details)
            See dbBact sequences/get_annotations REST API documentation
        term_info : dict of {term: info}
            where key (str) is the ontology term, info is a dict of details containing:
                'total_annotations' : total annotations having this term in the database
                'total_sequences' : number of annotations with this term for the sequence
        taxonomy : str
            the dbbact taxonomy of the sequence
        '''
        rdata = {}
        rdata['sequence'] = sequence
        res = self._get('sequences/get_annotations', rdata)
        if res.status_code != 200:
            logger.warn('error getting annotations for sequence %s' % sequence)
            return []
        res = res.json()
        annotations = res.get('annotations')
        if max_id is not None:
            annotations = [x for x in annotations if x['annotationid'] <= max_id]

        term_info = res.get('term_info')
        taxonomy = res.get('taxonomy')
        logger.debug('Found %d annotations for sequence %s' % (len(annotations), sequence))
        return annotations, term_info, taxonomy

    def get_annotation_string(self, cann):
        '''Get nice string summaries of annotation

        Parameters
        ----------
        cann : dict
            annotations of the output of get_seq_annotations()

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
            cdesc += 'HIGH IN '
            for cval in chigh:
                cdesc += cval + ' '
            cdesc += 'COMPARED TO '
            for cval in clow:
                cdesc += cval + ' '
            cdesc += 'IN '
            for cval in call:
                cdesc += cval + ' '
        elif cann['annotationtype'] == 'positive association':
            chigh = []
            call = []
            for cdet in cann['details']:
                if cdet[0] == 'all':
                    call.append(cdet[1])
                    continue
                if cdet[0] == 'high':
                    chigh.append(cdet[1])
                    continue
            cdesc += 'POSITIVELY ASSOCIATED WITH '
            for cval in chigh:
                cdesc += cval + ' '
            cdesc += 'IN '
            for cval in call:
                cdesc += cval + ' '
        elif cann['annotationtype'] == 'negative association':
            clow = []
            call = []
            for cdet in cann['details']:
                if cdet[0] == 'all':
                    call.append(cdet[1])
                    continue
                if cdet[0] == 'low':
                    chigh.append(cdet[1])
                    continue
            cdesc += 'NEGATIVELY ASSOCIATED WITH '
            for cval in clow:
                cdesc += cval + ' '
            cdesc += 'IN '
            for cval in call:
                cdesc += cval + ' '
        elif cann['annotationtype'] == 'isa':
            cdesc += ' IS A'
            for cdet in cann['details']:
                cdesc = cdesc + ' ' + cdet[1] + ','
        elif cann['annotationtype'] == 'contamination':
            cdesc += 'CONTAMINATION IN'
            for cdet in cann['details']:
                cdesc = cdesc + ' ' + cdet[1] + ','
        else:
            cdesc += cann['annotationtype'].upper() + ' IN'
            for cdet in cann['details']:
                cdesc = cdesc + ' ' + cdet[1] + ','
        if cann['description']:
            cdesc += ')'
        return cdesc

    def get_seq_annotation_strings(self, sequence, max_id=None, get_summary=True, sort_order=('other', 'contamination', 'diffexp', 'dominant', 'common')):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        sequence : str
            the DNA sequence to query the annotation strings about
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        get_summary: bool, optional
            True (default) to get summary of all annotations in the first 3 results of the function (see output details)
            False to just get the annotations
        sort_order: list of str or None, optional
            The order of the annotations based on the annotation type. All types not specified in the list are returned last
            If None, do not sort.

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'annotationid' : int, the annotation id in the database
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short text summary of the annotation
            NOTE: if get_summary=True, the first 3 descriptions are a summary of all the annotations. they include:
            shortdesc[0] - taxonomy for the sequence
            shortdesc[1] - 5 highest f-score terms for the sequence
            shortdesc[2] - 5 highest precision
        '''
        shortdesc = []
        annotations, term_info, taxonomy = self.get_seq_annotations(sequence, max_id=max_id)
        if get_summary:
            if taxonomy is None:
                taxonomy = 'taxonomy not availble'
            elif len(taxonomy) == 0:
                taxonomy = 'Taxonomy not available'
            shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, taxonomy))
            if len(term_info) > 0:
                # terms = []
                # for cterm, cinfo in term_info.items():
                #     terms.append([cterm, 42, cinfo.get('total_annotations')])
                #     # terms.append([cterm, cinfo.get('total_sequences'), cinfo.get('total_annotations')])
                # terms = sorted(terms, key=lambda x: x[1])
                # terms = sorted(terms, key=lambda x: -x[2])

                # get the most f-score terms (i.e. most identifying this feature)
                summary = 'most: '
                most_terms, prec_terms = self.get_common_annotation_term(annotations, term_info, num_common=5)
                for cmterm in most_terms:
                    summary += cmterm + ', '
                shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, summary))
                summary = 'special: '
                for cmterm in prec_terms:
                    summary += cmterm + ', '
                shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, summary))

            # # get the most f-score terms (i.e. most identifying this feature)
            # summary = 'most: '
            # most_terms = self.get_common_annotation_term(annotations, term_info, num_common=5)
            # for cmterm in most_terms:
            #     summary += cmterm + ', '
            # shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, summary))

            # # get the special terms (surprising)
            # summary = 'special: '
            # terms = sorted(terms, key=lambda x: x[2])
            # for cterm in terms[:min(4, len(terms))]:
            #     summary += '%s, ' % cterm[0]
            # shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, summary))

        # and get all the annotations for the feature
        # first prepare the sort order
        if sort_order is None:
            sort_order = []
        sort_order_dict = {}
        for idx, ctype in enumerate(sort_order):
            sort_order_dict[ctype] = idx
        max_idx = idx
        # now iterate over annotations and create a list with the sort order (for later sorting)
        annotation_sort_list = []
        for cann in annotations:
            cdesc = self.get_annotation_string(cann)
            ctype = cann.get('annotationtype', 'NA')
            annotation_sort_list.append([sort_order_dict.get(ctype, max_idx+1),(cann, cdesc)])
        annotation_sort_list = [x for pos,x in sorted(annotation_sort_list, key=lambda y: y[0])]
        shortdesc.extend(annotation_sort_list)

        return shortdesc

    def _get_common_seqs_for_terms(self, terms):
        '''Get a list of annotation ids for all annotations of type 'common' that contain all the terms in terms
        Used in background enrichment analysis

        Parameters
        ----------
        terms: list of str
            the dbBact terms to search for (returned annotations must contain all the terms)

        Returns
        -------
        ids: list of int
            the annotation ids for all the annotations of type common that contain all the terms
        '''
        num_res = 0
        ids = set()
        for cterm in terms:
            logger.debug('getting annotations for term: %s' % cterm)
            res = self._get('ontology/get_annotations', rdata={'term': cterm}, param_json=False)
            # res=requests.get(server+'/ontology/get_annotations',params={'term':cterm})
            anno = res.json()['annotations']
            cids = [canno['id'] for canno in anno if canno['annotationtype'] == 'common']
            logger.debug('Found %d common annotations for term: %s' % (len(cids), cterm))
            if num_res == 0:
                ids.update(cids)
            else:
                ids = ids.intersection(set(cids))
            num_res += 1
        logger.debug('found %d common annotations total for all %d terms' % (len(ids), len(terms)))
        return ids

    def _get_sequences(self, annotation_ids, min_appearance=2):
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
            res = self.get_annotation_sequences(cid)
            for cseq in res:
                seqs_ids[cseq] += 1
        # keep only sequences appearing in enough annotations
        ok_seqs = [ck for ck, cv in seqs_ids.items() if cv >= min_appearance]
        res = self.get_seq_id_info(ok_seqs)
        seqs = [x.upper() for x in res]
        # create the inflated list of sequences (each sequence appears according to the number of annotations it is found in)
        for idx, x in enumerate(res):
            seqs_list.extend([x.upper()] * seqs_ids[ok_seqs[idx]])
        logger.debug('found %d sequences associated with annotations' % len(seqs))
        return seqs, seqs_list

    def _get_exp_annotations(self, annotations):
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
                exp_annotations[cexpid][cterm] += 1
        return exp_annotations

    def _get_all_term_counts(self, features, feature_annotations, annotations, ignore_exp=None, score_method='all_mean', use_term_pairs=False):
        '''Get the term scores for each feature in features

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
        Use_term_pairs: bool, optional
            True to get all term pair annotations (i.e. human+feces etc.)

        Returns
        -------
        dict of {feature:str, list of tuples of (term:str, score:float)}
            key is feature, value is list of (term, score)
        '''
        # set up the list of annotations per experiment if needed
        if score_method == 'all_mean':
            exp_annotations = self._get_exp_annotations(annotations)
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
            feature_terms[cfeature] = self.get_annotation_term_counts(annotation_list, exp_annotations=exp_annotations, score_method=score_method, use_term_pairs=use_term_pairs)
            # print(feature_terms)
        return feature_terms

    def get_annotation_term_counts(self, annotations, exp_annotations=None, score_method='all_mean', use_term_pairs=False):
        '''Get the annotation type corrected count (score) for all terms in annotations

        Parameters
        ----------
        annotations : list of dict
            list of annotations where the feature is present
        exp_annotations: dict of {expid:int : dict of {term:str : total:int}}
            from self._get_exp_annotations()
        score_method: str (optional)
                The method to use for score calculation:
                'all_mean' : score is the mean (per experiment) of the scoring for the term out of all annotations that have the term
                'sum' : score is sum of all annotations where the term appears (experiment is ignored)
        use_term_pairs: bool, optional
            True to get all term pair annotations (i.e. human+feces etc.)

        Returns
        -------
            list of tuples (term, count)
        '''
        term_count = defaultdict(int)

        # print('processing %d annotations' % len(annotations))
        for cannotation in annotations:
            details = cannotation['details']
            annotation_type = cannotation['annotationtype']
            if annotation_type == 'common':
                ascore = 1
            elif annotation_type == 'dominant':
                ascore = 2
            elif annotation_type == 'other':
                ascore = 0.5
            elif annotation_type == 'contamination':
                ascore = 1
                details = [('all', 'contamination')]
            elif annotation_type == 'diffexp':
                ascore = 1
            elif annotation_type == 'positive correlation':
                ascore = 1
            elif annotation_type == 'negative correlation':
                ascore = 1
            else:
                logger.warn('unknown annotation type %s encountered (annotationid %d). skipped' % (annotation_type, cannotation['annotationid']))
                continue

            # temp_details = deepcopy(details)
            temp_details = sorted(details, key=lambda x: x[1],)
            # print(temp_details)

            # if we need to create the pairs, let's add them to the details
            if use_term_pairs:
                if len(details) < 10:
                    new_details = []
                    for p1 in range(len(details)):
                        # print('now detail term idx %d' % p1)
                        for p2 in range(p1 + 1, len(details)):
                            det1 = details[p1]
                            det2 = details[p2]
                            term1 = det1[1]
                            term2 = det2[1]
                            type1 = det1[0]
                            type2 = det2[0]
                            if type1 == 'low':
                                term1 = '-' + term1
                            if type2 == 'low':
                                term2 = '-' + term2
                            cnew_type = 'all'
                            if type1 == type2:
                                cnew_type == type1
                            cnew_term = '%s+%s' % (term1, term2)
                            new_details.append((cnew_type, cnew_term))
                    temp_details.extend(new_details)
                    # print('new details: %d' % len(details))

            for cdetail in temp_details:
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
                        # print('term %s, exp %d, annotationid %d, score %d, scale %d' % (cterm, cannotation['expid'], cannotation['annotationid'], cscore, scale_factor))
                        scale_factor = 1
                else:
                    scale_factor = 1

                # fix for differential abundance
                term_count[cterm] += ascore * cscore / scale_factor

        res = []
        for k, v in term_count.items():
            # flip and add '-' to term if negative
            if v < 0:
                k = '-' + k
                v = -v
            res.append((k, v))
        return res

    def get_common_annotation_term(self, annotations, term_info, num_common=5):
        '''Get the most common term when combining all annotaions

        Parameters
        ----------
        annotations : list of dict
            The annotations to get the most common term from
        term_info : dict (key is term, value is dict)
            information about each term (total_sequences, total_annotations)
        num_common : int (optional)
            The number of common terms to return. if 0, return all

        Returns
        -------
        list of str
            list of the num_common most common terms
        '''
        # convert the annotations to the compact form. we need a dict of all annotations (key is id)
        annotations_dict = {}
        for cannotation in annotations:
            cid = cannotation['id']
            annotations_dict[str(cid)] = cannotation
        fscore, recall, precision, term_count, reduced_f = get_enrichment_score(annotations_dict, [(1, list(annotations_dict.keys()))], term_info=term_info)
        sorted_by_f = sorted(fscore, key=fscore.get, reverse=True)
        sorted_by_prec = sorted(precision, key=precision.get, reverse=True)
        if num_common > 0:
            sorted_by_f = sorted_by_f[:num_common]
            sorted_by_prec = sorted_by_prec[:num_common]
        return(sorted_by_f, sorted_by_prec)

    def find_experiment_id(self, datamd5='', mapmd5='', getall=False):
        """
        find the data id for the data/map md5 (which are calculated on load)
        note the md5s don't change following filtering/normalization/etc... - only the original data

        Parameters
        ----------
        datamd5 : str
            from Experiment.datamd5
        mapmd5 : str
            from Experiment.mapmd5
        getall : bool (optional)
            False (default) to get only 1st id, True to get a list of all

        Returns
        -------
        expids: int (if getall=False - default) or list of int (if getall=True)
            an id or a list of ids of matching dataID indices (or None if no match)
        """
        logger.debug('findexpid for datamd5 %s mapmd5 %s' % (datamd5, mapmd5))
        details = []
        if datamd5:
            details.append(['DataMD5', datamd5])
        if mapmd5:
            details.append(['MapMD5', mapmd5])
        if not details:
            logger.warn('Error. MapMD5 and DataMD5 both missing from find_experiment_id')
            return None

        rdata = {}
        rdata['details'] = details
        res = self._get('experiments/get_id', rdata)
        if res.status_code == 200:
            expids = res.json()['expId']
            if len(expids) == 0:
                logger.warn('No experiment found matching the details %s' % details)
                return None
            if not getall:
                if len(expids) > 1:
                    logger.warn('Problem. Found %d matches for data' % len(expids))
                logger.debug('Found experiment id %d' % expids[0])
                return expids[0]
            logger.debug("Found %d matches to data" % len(expids))
            return expids
        logger.error('Error getting expid from details')
        return None

    def get_experiment_info(self, expid):
        """
        get the information about a given experiment dataid

        Parameters
        ----------
        dataid : int
            The dataid on the experiment (DataID field)

        Returns
        -------
        info : list of (str,str)
            list of tuples for each entry in the experiment:
            (type,value) about dataid (i.e. ('PubMedID','100234'))
            empty if dataid not found
        """
        logger.debug('get experiment details for expid %d' % expid)
        rdata = {}
        rdata['expId'] = expid
        res = self._get('experiments/get_details', rdata)
        if res.status_code == 200:
            details = res.json()['details']
            logger.debug('Found %d details for experiment %d' % (len(details), expid))
            return details
        return []

    def get_experiment_annotations(self, expid):
        """Get the list of annotations for experiment expid

        Parameters
        ----------
        expid : int
            The dataid of the experiment

        Returns
        -------
        info: list of str
            the list of annotations for this experiment (1 item per annotation)
        """
        logger.debug('get experiment annotations for expid %d' % expid)
        rdata = {}
        rdata['expId'] = expid
        res = self._get('experiments/get_annotations', rdata)
        if res.status_code != 200:
            logger.warn('error getting annotations for experiment %d' % expid)
            return []
        annotations = res.json()['annotations']
        logger.debug('Found %d annotations for experiment %d' % (len(annotations), expid))
        # make it into a nice list of str
        info = []
        for cann in annotations:
            cstr = 'date:%s description:%s user:%s private:%s' % (cann['date'], cann['description'], cann['userid'], cann['private'])
            info.append(cstr)
        return info

    def add_experiment_data(self, data, expid=None):
        """add new data entries (for a new experiment)

        Parameters
        ----------
        data : list of tuples (Type:Value)
            a list of tuples of (Type,Value) to add to Data table (i.e. ("PUBMEDID","322455") etc)
        expid : list of int
            the ids in which this experiment appears (from find_experiment_id)

        Returns
        -------
        suid : int
            the value of DataID for the new experiment (from Data table)
        """
        # we need to get a new identifier for all entries in the experiment
        # there should be a more elegant way to do it
        logger.debug("add_experiment_data for %d enteries" % len(data))
        if expid is None:
            # add new experiment
            logger.debug("add_experiment_data for a new experiment")
        else:
            logger.debug('add_experiment_data for existing experiment %d' % expid)
        rdata = {}
        rdata['expId'] = expid
        rdata['details'] = data
        res = self._post('experiments/add_details', rdata)
        if res.status_code == 200:
            newid = res.json()['expId']
            logger.debug('experiment added. id is %d' % newid)
            return newid
        else:
            logger.error('error adding experiment. msg: %s' % res.content)
            return None

    def add_annotations(self, expid, sequences, annotationtype, annotations, submittername='NA',
                        description='', method='na', primerid='V4', agenttype='Calour', private='n'):
        """
        Add a new manual curation to the database

        Paramerers
        ----------
        expid : int
            the experiment ID - the value of ExpID from Experiments table
        sequences : list of str
            the list of DNA sequences to curate
        annotationtype : str
            the annotation type. can be:
            'COMMON' : the sequences are common in the samples (present in >0.5 of the samples)
            'DIFFEXP' : the sequences are different between two conditions in the samples
            'CONTAM' : the sequences are suspected experimental contaminants
            'DOMINANT' : the sequences are found in mean high frequency (>1%) in the samples
            'OTHER' : other general observations (i.e. the sequences are known pathogens)
            'POSITIVE CORRELATION': the phenotype is positively correlated with the sequence (i.e. high growth rate when bacteria is present)
            'NEGATIVE CORRELATION': the phenotype is negatively correlated with the sequence (i.e. low growth rate when bacteria is present)
        annotations : list of Type, Value
            The annotations to add to the AnnotationsList table (Type,Value).
            Value is the ontology term to add
            Type can be:
            'ALL' : the ontology term applies to all samples
            'LOW' : the sequences are present less in samples with the ontology term
            'HIGH' : the sequences are present more in samples with the ontology term
        submittername : str (optional)
            Name of the submitter (first,last) or NA (default)
        description : str (optional)
            text description of the annotation entry (i.e. "lower in whole wheat pita bread"). default is ''
        method : str (optional)
            text description of how the annotation was detected (mostly for 'DIFFEXP' annotations). default is ''
        primerid : str (optional)
            the primer region used for the sequence (usually 'v4')
        agenttype : str (optional)
            the program submitting the curation (defaul='Calour')
        private : str (optional)
            'n' (default) for public annotation, or 'y' for private annotation (visible only to the user submitting)

        Returns
        -------
        err: str
            if not empty (''), error enoucntered during adding annotation, annotation was not added
        annotationid : int or None
            the AnnotationID (in Annotations table) of the new annotation, or None if not added
        """
        logger.debug('addannotation - %d sequences' % len(sequences))
        if len(sequences) == 0:
            msg = 'No sequences to annotate!'
            logger.warn(msg)
            return msg, None
        if len(annotations) == 0:
            msg = 'No annotations to add!'
            logger.warn(msg)
            return msg, None

        # add the annotation
        rdata = {}
        rdata['expId'] = expid
        rdata['sequences'] = sequences
        rdata['region'] = primerid
        rdata['annotationType'] = annotationtype
        rdata['method'] = method
        rdata['agentType'] = agenttype
        rdata['description'] = description
        rdata['private'] = private
        rdata['annotationList'] = annotations

        res = self._post('annotations/add', rdata)
        if res.status_code == 200:
            newid = res.json()['annotationId']
            logger.debug('Finished adding experiment id %d annotationid %d' % (expid, newid))
            return '', newid
        logger.warn('problem adding annotations for experiment id %d' % expid)
        logger.warn(res.content)
        return res.content, None

    def get_seq_list_fast_annotations(self, sequences, get_taxonomy=False, get_parents=False, get_term_info=True, max_id=None):
        '''Get annotations for all sequences in list using compact format and with parent ontology terms

        Params
        ------
        sequences : list of str
            The list of DNA sequences to get the annotations for
        get_taxonomy: bool, optional
            True to get taxonomy info for each sequence using the dbbact assigned taxonomy (supplied in the 'taxonomy' dict)
        get_parents: bool, optional
            True to get all ontology parents for each annotation term.
        get_term_info: bool, optional
            True to get total annotation/experiments for each term appearing in the annotations
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id

        Returns
        -------
        sequence_terms : dict of {sequence: list of terms}
            key is sequence, value is list of ontology terms present in the bacteria.
        sequence_annotations : dict of {sequence: list of annotationIDs}
            key is sequence, value is list of annotationIDs present in the bacteria.
        annotations : dict of {annotationID : annotation_details}
            key is annotaitonID (int), value is the dict of annotation details.
        term_info : dict of {term, dict}:
            Information about each term which appears in the annotation parents. Key is the ontolgy term. the value dict is:
                'total_annotations' : int
                    total number of annotations where this term appears (as a parent)
                'total_sequences' : int
                    total number of sequences in annotations where this term appears (as a parent)
        taxonomy: dict of {sequence(str): taxonomy(str)}
            the dbbact assigned taxonomy for each sequence
        '''
        logger.info('Getting dbBact annotations for %d sequences, please wait...' % len(sequences))
        rdata = {}
        rdata['sequences'] = list(sequences)
        rdata['get_term_info'] = get_term_info
        rdata['get_taxonomy'] = get_taxonomy
        rdata['get_parents'] = get_parents
        rdata['get_all_exp_annotations'] = False
        res = self._get('sequences/get_fast_annotations', rdata)
        if res.status_code != 200:
            logger.warning('error getting fast annotations for sequence list. got status code %s' % res.status_code)
            raise ValueError('error getting fast annotations for sequence list. got status code %s' % res.status_code)
        res = res.json()
        logger.info('got %d annotations' % len(res['annotations']))

        sequence_terms = {}
        sequence_annotations = {}
        for cseq in sequences:
            sequence_terms[cseq] = []
            sequence_annotations[cseq] = []
        for cseqannotation in res['seqannotations']:
            cpos = cseqannotation[0]
            # need str since json dict is always string
            cseq = sequences[cpos]

            seq_annotation_ids = cseqannotation[1]
            # remove sequence annotations for annotations > max_id
            if max_id is not None:
                seq_annotation_ids = [x for x in seq_annotation_ids if int(x) <= max_id]
            sequence_annotations[cseq].extend(seq_annotation_ids)
            for cannotation in cseqannotation[1]:
                for k, v in res['annotations'][str(cannotation)]['parents'].items():
                    if k == 'high' or k == 'all':
                        for cterm in v:
                            sequence_terms[cseq].append(cterm)
                    elif k == 'low':
                        for cterm in v:
                            sequence_terms[cseq].append('-' + cterm)

        annotations = res['annotations']
        # replace the string in the key with an int (since in json key is always str)
        keys = list(annotations.keys())
        for cid in keys:
            annotations[int(cid)] = annotations.pop(cid)

        # remove the annotations with id > max_id
        if max_id is not None:
            ignored = 0
            keys = list(annotations.keys())
            for cid in keys:
                if cid > max_id:
                    annotations.pop(cid)
                    ignored += 1
            logger.warning('ignoring %d annotation with id > max_id %d' % (ignored, max_id))

        total_annotations = 0
        for cseq_annotations in sequence_annotations.values():
            total_annotations += len(cseq_annotations)
        logger.info('Got %d annotation-sequence pairs' % total_annotations)

        taxdict = {}
        if 'taxonomy' in res:
            if len(res['taxonomy']) > 0:
                for idx, cseq in enumerate(sequences):
                    taxdict[cseq] = res['taxonomy'][idx]
        if 'term_info' not in res:
            res['term_info'] = {}

        return sequence_terms, sequence_annotations, res['annotations'], res['term_info'], taxdict

    def delete_annotation(self, data):
        '''Delete the annotation from the database
        Tries to delete the annotation using the username/password from the config file.
        A user can only delete the annotations he created or anonymous annotations.
        Otherwise, an error will be returned

        Parameters
        ----------
        data : dict
            The annotationdetails to delete (using the 'annotationid' key)

        Returns
        -------
        str:
            empty if ok, non empty string if error encountered
        '''
        if 'annotationid' not in data:
            return 'No annotationID for selected annotation'
        annotationid = data['annotationid']
        rdata = {}
        rdata['annotationid'] = annotationid
        res = self._post('annotations/delete', rdata)
        if res.status_code != 200:
            msg = 'Deletion failed. error %s' % res.content
            logger.warn(msg)
            return msg
        logger.info('annotation %d deleted' % annotationid)
        return ''

    def get_user_id(self):
        '''Get the userid for the username/password
        Uses self.username, self.password

        Returns
        -------
        int or None
        the userid (0 is anonymous) or None if error encountered
        '''
        res = self._post('users/get_user_id', rdata={})
        if res.status_code != 200:
            logger.debug('get_user_id failed. User/pwd incorrect')
            return None
        res = res.json()
        return res.get('user')

    def register_user(self, user, pwd, email='', description='', publish='n', name=''):
        '''Register a new user in the database

        Parameters
        ----------
        user : str
            the username
        pwd : str
            the password for the user
        email : str (optional)
            the recovery and notification email
        description : str (optional)
            type of user (can be researcher/curious/etc)
        publish : str (optional)
            'n' (default) to hide the email from other users,
            'y' to enable other users to see the email
        name : str (optional)
            name of the user

        Returns
        -------
        err : str
            empty if ok or the error encountered
        '''
        rdata = {}
        rdata['email'] = email
        rdata['description'] = description
        rdata['name'] = name
        rdata['publish'] = publish
        self.username = user
        self.password = pwd

        res = self._post('users/register_user', rdata=rdata)
        if res.status_code != 200:
            msg = 'Error registering user: %s' % res.content
            logger.debug(msg)
            return msg, 0
        logger.debug('new user %s registered' % user)
        return ''

    def upadte_annotation(self, annotation, exp=None):
        '''Update an existing annotation

        Parameters
        ----------
        annotation : dict
            The annotation to update (keys/values are database specific)
        exp : ``Experiment`` (optional)
            The calour experiment from which the annotation is coming from
        Returns
        -------
        str
            empty if ok, otherwise the error encountered
        '''
        from . import dbannotation
        dbannotation.update_annotation_gui(self, annotation, exp)

    def send_update_annotation(self, annotationid, annotationtype=None, annotations=None, description=None, method=None, private=None):
        '''Update an annotation in the database

        None indicates field will not be updated

        Returns
        -------
        str
            empty if ok, error string if error encountered
        '''
        rdata = {}
        rdata['annotationId'] = annotationid
        rdata['annotationType'] = annotationtype
        rdata['annotationList'] = annotations
        rdata['description'] = description
        rdata['method'] = method
        rdata['private'] = private
        res = self._post('annotations/update', rdata)
        if res.status_code == 200:
            logger.debug('Finished updating annotationid %d' % annotationid)
            return ''
        msg = 'problem updating annotationid %d: %s' % (annotationid, res.content)
        logger.warn(msg)
        return msg

    def _get_all_annotation_string_counts(self, features, all_annotations, seq_annotations, ignore_exp=None):
        '''Get the annotations present for each feature

        Parameters
        ----------
        features: list of str
            the features to get the annotations for
        all_annotation: dict of {annotationid: annotationd(dict)}
            the details of each annotation, keyed by annotationid
        seq_annotations dict of {sequence:str, annotations_list(list of int)}
            the annotationids associated with each sequence
        ignore_exp: list of int or None (optional)
            the experimentIDs to ignore when getting the annotations

        Returns
        -------
        dict of {sequence(str): list of annotation details(str)}
            the annotations descriptions (text) associated with this sequence
        '''
        feature_annotations = {}
        for cseq, annotations_list in seq_annotations.items():
            if cseq not in features:
                continue
            newdesc = []
            for cannotation in annotations_list:
                if ignore_exp is not None:
                    if all_annotations[cannotation]['expid'] in ignore_exp:
                        continue
                cdesc = self.get_annotation_string(all_annotations[cannotation])
                newdesc.append((cdesc, 1))
            feature_annotations[cseq] = newdesc
        return feature_annotations

    def _filter_enriched_exps(self, orig_term_list, term_info, g1featurelist, g2featurelist, seq_annotations, all_annotations, ignore_exp, min_appearances, num_results_needed):
        '''get the terms enriched in enough experiments, at most getting num_results_needed experiments

        Parameters
        ----------
        orig_term_list:
        g1featurelist, g2featurelist:
            list of the features in the 2 groups
        seq_annotations:
        all_annotations:
        ignore_exp:
        min_appearances: int
            the minimal number of experiments the term needs to be enriched in
        num_results_needed: int
            the maximal number of enriched terms needed

        Returns
        -------
        term_idx: list of int
            the indices of the enriched experiments
        num_enriched_exps: list of int
            the number of enriched experiments for the term
        num_total_exps:
            the total number of experiments with the term
        '''
        term_idx = []
        num_enriched_exps = []
        num_total_exps = []
        num_found = 0

        for cidx, cterm in enumerate(orig_term_list):
            enriched_exps, enriched_annotations, total_exp_annotations = self.count_enriched_exps(cterm, g1featurelist, g2featurelist, seq_annotations, all_annotations, ignore_exp=ignore_exp)

            # if this term is enriched in enough experiments (>= min_experiments)
            if len(enriched_exps) >= min_appearances:
                num_found += 1
                term_idx.append(cidx)
                num_enriched_exps.append(len(enriched_exps))
                if term_info is not None:
                    if cterm in term_info:
                        num_total_exps.append(term_info[cterm]['total_experiments'])
                    else:
                        num_total_exps.append(0)
                else:
                    num_total_exps.append(-1)

                if num_found >= num_results_needed:
                    logger.debug('found enough')
                    break
        return term_idx, num_enriched_exps, num_total_exps

    def term_enrichment(self, g1_features, g2_features, all_annotations, seq_annotations, term_info=None, term_type='term', ignore_exp=None, min_appearances=0,
                        fdr_method='dsfdr', score_method='all_mean', random_seed=None, use_term_pairs=False, alpha=0.1,
                        method='meandiff', transform_type='rankdata', numperm=1000, min_exps=1, add_single_exp_warning=True, num_results_needed=0, focus_terms=None):
        '''Get the list of enriched terms in features compared to all features in exp.

        given uneven distribtion of number of terms per feature

        Parameters
        ----------
        g1_features : list of str
        g2_features : list of str
            The sequences to test for enrichmnt
        all_annotations :
            dict of {annotationid (int): annotation (dict)}
            all annotations for the features (from get_annotations_compact)
        seq_annotations: dict of {sequence (str): annnotationids (list of int)}
        term_info: dict of {term(str):info(dict)} or None, optional
            Not none to enable total experiment count for each term (returned in the 'num_total_exps' field in first returned pandas dataframe)
        term_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
            'parentterm' - ontology terms including parent terms associated with each feature.
            'annotation' - the full annotation strings associated with each feature
            'combined' - combine 'term' and 'annotation'
            'sigterm' - count number of significant experiments with the term
        ignore_exp: list of int or None or True(optional)
            List of experiments to ignore in the analysis
            None (default) to use annotations from all experiments including the current one
        min_appearances : int (optional)
            The minimal number of experiments a term is enriched in, in order to include in output.
            NOTE: can slow down computation! use with nun_results_needed to speed up.
        fdr_method : str (optional)
            The FDR method used for detecting enriched terms (permutation test). options are:
            'dsfdr' (default): use the discrete FDR correction
            'bhfdr': use Benjamini Hochbert FDR correction
        transform_type: str, optional
            the data trasformation to use before calculating the dsFDR. For options see dsfdr.dsfdr()
        score_method : str (optional)
            The score method used for calculating the term score. Options are:
            'all_mean' (default): mean over each experiment of all annotations containing the term
            'sum' : sum of all annotations (experiment not taken into account)
        method: str or None, optional
            'card_mean': use a null model keeping the number of annotations per each bacteria (useful for comparison to background sequences). NOTE: use transfrom_type=None for this method.
            otherwise, use as method for dsfdr.dsfdr() (i.e. 'meandiff')
        random_seed: int or None
            int to specify the random seed for numpy.random.
        use_term_pairs: bool, optional
            True to also test enrichment in pairs of terms (i.e. homo sapiens+feces, etc.)
        min_exps: int, optional
            the minimal number of experiments a term appears in in order to include in results (default = 1 so all are shown)
        add_single_exp_warning: bool, optional
            True to add a warning to terms present in one experiment, False to not add this warning
        num_results_needed: int, optional
            if min_appearances supplied, use this to speed up calculation by looking at terms in sorted (effect size) order until num_results_needed are met at each end
            0 to calculate for all
            NOTE: using this keeps only num_results_needed significant terms!
        focus_terms: list of str or None, optional
            if not None, use only annotations containing all the terms in focus_terms

        Returns
        -------
        pandas.DataFrame with  info about significantly enriched terms. columns:
            feature : str the feature
            pval : the p-value for the enrichment (float)
            qval: the FDR corrected qvalue for the enrichment (float)
            odif : the effect size (float)
            observed : the number of observations of this term in group1 (int)
            expected : the expected (based on all features) number of observations of this term in group1 (float)
            frac_group1 : fraction of total terms in group 1 which are the specific term (float)
            frac_group2 : fraction of total terms in group 2 which are the specific term (float)
            num_group1 : number of total terms in group 1 which are the specific term (float)
            num_group2 : number of total terms in group 2 which are the specific term (float)
            num_enriched_exps: number of experiments where the term is significantly enriched
            num_enriched_annotations: number of annotations where the term is significantly enriched
            num_total_exps: number of experiments with this annotation in the term annotations list
            description : the term (str)
        numpy.Array where rows are features (ordered like the dataframe), columns are terms and value is score
            for term in feature
        pandas.DataFrame with info about the features used. columns:
            group: int the group (1/2) to which the feature belongs
            sequence: str
        '''
        # set up the random seed if needed
        if random_seed is not None:
            np.random.seed(random_seed)

        # prepare arrays of the 2 features and the combined features
        # we sort the features to ensure replicability
        g1_features = np.array(sorted(list(g1_features)))
        g2_features = np.array(sorted(list(g2_features)))
        exp_features = np.hstack([g1_features, g2_features])

        # filter keeping only annotations containing focus_terms (if not None)
        if focus_terms is not None:
            logger.debug('filtering annotations for %d focus terms.' % len(focus_terms))
            seq_annotations = seq_annotations.copy()
            focus_terms = set(focus_terms)
            ok_annotations = set()
            for cid, cannotation in all_annotations.items():
                found_terms = set()
                for cdetail in cannotation['details']:
                    if cdetail[1] in focus_terms:
                        found_terms.add(cdetail[1])
                if len(found_terms) == len(focus_terms):
                    ok_annotations[cid] = cannotation
            logger.debug('Keeping %d annotations out of %d original' % (len(ok_annotations), len(all_annotations)))
            for k, v in seq_annotations.items():
                nv = [x for x in v if x in ok_annotations]
                seq_annotations[k] = nv

        # Get the feature_terms dict according to the scoring method
        # dict of {feature:str, list of tuples of (term:str, score:float)}, key is feature, value is list of (term, score)
        if term_type == 'term':
            # score for each term
            logger.debug('getting "terms" per feature')
            feature_terms = self._get_all_term_counts(exp_features, seq_annotations, all_annotations, ignore_exp=ignore_exp, score_method=score_method, use_term_pairs=use_term_pairs)
        elif term_type == 'parentterm':
            # score for each term including the parent terms (based on the ontology structure) for each term
            raise ValueError('parentterm not supported yet')
        elif term_type == 'annotation':
            # score for each annotation
            logger.debug('getting "annotations" per feature')
            feature_terms = self._get_all_annotation_string_counts(exp_features, all_annotations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
            # each annotation by definition is in one experiment, so no need to warn about single annotation :)
            add_single_exp_warning = False
        elif term_type == 'combined':
            # score is for each term or annotation
            logger.debug('getting combined "terms" and "annotations" per feature')
            feature_terms = self._get_all_term_counts(exp_features, seq_annotations, all_annotations, ignore_exp=ignore_exp, score_method=score_method, use_term_pairs=use_term_pairs)
            feature_annotations = self._get_all_annotation_string_counts(exp_features, all_annotations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
            for cfeature, cvals in feature_annotations.items():
                if cfeature not in feature_terms:
                    feature_terms[cfeature] = []
                feature_terms[cfeature].extend(cvals)
            # each annotation by definition is in one experiment
            add_single_exp_warning = False
        elif term_type == 'sigterm':
            # feature_terms = self._get_all_annotation_string_counts(exp_features, all_annotations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
            # get an str of each annotationid
            logger.debug('getting "sigterm" (number of experiment with the term) per feature')
            feature_terms = {}
            for cseq, cannotations in seq_annotations.items():
                feature_terms[cseq] = [(str(x), 1) for x in cannotations]
            # each annotation by definition is in one experiment
            add_single_exp_warning = False
        else:
            raise ValueError('term_type %s not supported for dbbact. possible values are: "term", "parentterm", "annotation", "sigterm', "combined")

        # arrays of terms (rows) x bacteria (cols)
        feature_array, term_list = self._get_term_features(g1_features, feature_terms)
        bg_array, term_list = self._get_term_features(g2_features, feature_terms)

        # and the array of terms (rows) x all bacteria (in both groups) (cols)
        all_feature_array = np.hstack([feature_array, bg_array])

        import hashlib
        # print(hashlib.sha256(keep.to_csv().encode()).hexdigest())
        print(hashlib.sha256(feature_array).hexdigest())
        print(hashlib.sha256(bg_array).hexdigest())
        print(hashlib.sha256(all_feature_array).hexdigest())

        # remove terms present in a small number of experiments (<min_exps)
        # get the number of experiments for each term
        term_exps = defaultdict(set)
        for cannotation in all_annotations.values():
            for cdetail in cannotation['details']:
                cterm = cdetail[1]
                term_exps[cterm].add(cannotation['expid'])
        # remove the ones appearing in not enough experiments
        remove_set = set()
        for cterm, cexps in term_exps.items():
            if len(cexps) < min_exps:
                remove_set.add(cterm)
                remove_set.add('-' + cterm)
        remove_pos = [pos for pos, term in enumerate(term_list) if term in remove_set]
        num_removed = len(remove_set)
        logger.debug('removing %d terms present in  <%d experiments' % (num_removed, min_exps))
        term_list = [x for x in term_list if x not in remove_set]
        all_feature_array = np.delete(all_feature_array, remove_pos, axis=0)

        # fix names ("LOWER IN" instead of "-") and add info about single experiment
        new_term_list = []
        orig_term_list = term_list.copy()
        for cterm in term_list:
            if cterm[0] == '-':
                ccterm = 'LOWER IN ' + cterm[1:]
                cterm = cterm[1:]
            else:
                ccterm = cterm
            if len(term_exps[cterm]) == 1:
                if add_single_exp_warning:
                    # cterm = '**%s**%s' % (list(term_exps[ccterm])[0], ccterm)
                    # ccterm = '%s {*single exp %s*}' % (ccterm, list(term_exps[cterm])[0])
                    ccterm = '%s *' % (ccterm)
            new_term_list.append(ccterm)
        term_list = new_term_list

        # set up the labels for the dsFDR test. 1s for features in the first group, 0 for the second
        labels = np.zeros(all_feature_array.shape[1])
        labels[:feature_array.shape[1]] = 1

        logger.debug('got %d labels=0, %d labels=1' % (np.sum(labels == 0), np.sum(labels == 1)))

        # find which terms are significantly enriched in either of the two feature groups
        if method == 'card_mean':
            if transform_type is not None:
                logger.warn('Using card_mean method for dbBact term enrichment but transform_type is not None is not recommended. Set transform_type to None')
            # count the total number of annotations per feature
            per_feature_annotations = np.zeros([len(exp_features)])
            for idx, cfeature in enumerate(exp_features):
                per_feature_annotations[idx] = len(seq_annotations[cfeature])

            # create the score function so that it normalizes each term by the total number of annotations in all features in the group
            def _mean_diff_toal_annotations(data, labels, per_feature_annotations=per_feature_annotations):
                mean0 = np.mean(data[:, labels == 0], axis=1) / np.mean(per_feature_annotations[labels == 0])
                mean1 = np.mean(data[:, labels == 1], axis=1) / np.mean(per_feature_annotations[labels == 1])
                tstat = mean1 - mean0
                return tstat

            # get the significantly enriched terms using the total annotations normalizing function
            logger.debug('card mean transfrom_type=%s' % transform_type)
            keep, odif, pvals, qvals = dsfdr(all_feature_array, labels, method=_mean_diff_toal_annotations, transform_type=transform_type, alpha=alpha, numperm=numperm, fdr_method=fdr_method, random_seed=random_seed)
        else:
            logger.debug('enrichment test method=%s transform_type=%s' % (method, transform_type))
            keep, odif, pvals, qvals = dsfdr(all_feature_array, labels, method=method, transform_type=transform_type, alpha=alpha, numperm=numperm, fdr_method=fdr_method, random_seed=random_seed)

        # keep only terms passing the FDR corrected significance threshold alpha (by using the keep field from dsFDR)
        keep = np.where(keep)[0]
        if len(keep) == 0:
            logger.info('no enriched terms found')
        term_list = np.array(term_list)[keep]
        orig_term_list = np.array(orig_term_list)[keep]
        all_feature_array = all_feature_array[keep, :].T
        all_feature_array = all_feature_array * 100
        odif = odif[keep]
        pvals = pvals[keep]
        qvals = qvals[keep]
        logger.debug('keeping %d significant terms' % len(pvals))

        # if sigterm, convert enriched annotations to the terms in them and count over terms
        if term_type == 'sigterm':
            new_term_list = []
            new_odif = []
            new_pvals = []
            exp_terms = defaultdict(set)
            # go over the id strings and count what terms are enriched in what experiments
            for cidstr in term_list:
                cann = all_annotations[int(cidstr)]
                cexp = cann['expid']
                for ctype, cterm in cann['details']:
                    if ctype == 'low':
                        cterm = '-' + cterm
                    exp_terms[cterm].add(cexp)
            for cterm, cexp_list in exp_terms.items():
                new_term_list.append(cterm)
                new_odif.append(len(cexp_list))
                new_pvals.append(0)
            new_all_feature_array = np.zeros([all_feature_array.shape[0], len(new_term_list)])
            # NEED TO ADD FILLING IT UP
            all_feature_array = new_all_feature_array
            term_list = np.array(new_term_list)
            odif = np.array(new_odif)
            pvals = np.array(new_pvals)

        # sort according to effect size
        si = np.argsort(odif)
        term_list = term_list[si]
        orig_term_list = orig_term_list[si]
        odif = odif[si]
        pvals = pvals[si]
        qvals = qvals[si]
        all_feature_array = all_feature_array[:, si]

        # get the per-term enriched experiments count if needed
        num_enriched_exps = np.zeros(len(term_list)) - 1
        num_total_exps = np.zeros(len(term_list)) - 1
        if min_appearances > 0:
            logger.debug('Getting per-experiment enrichment for min_appearances %d' % min_appearances)

            keep_min_exps = []
            g1featurelist = list(g1_features)
            g2featurelist = list(g2_features)
            num_enriched_exps = []
            num_total_exps = []
            # if num_results_needed == 0 then we want to scan all terms and not stop before
            if num_results_needed == 0:
                num_results_needed = len(orig_term_list)

            # find the negative terms
            npos = np.where(odif < 0)[0]
            ctermlist = orig_term_list[npos]
            term_idx, num_enriched_exps, num_total_exps = self._filter_enriched_exps(ctermlist, term_info, g1featurelist, g2featurelist, seq_annotations, all_annotations, ignore_exp, min_appearances, num_results_needed)
            logger.info('found %d negative' % len(term_idx))
            keep_min_exps = npos[term_idx]

            # find the positive terms.
            npos = np.where(odif > 0)[0]
            # reverse it since biggest effect size is at the end
            npos = npos[::-1]
            ctermlist = orig_term_list[npos]
            term_idx, cnum_enriched_exps, cnum_total_exps = self._filter_enriched_exps(ctermlist, term_info, g1featurelist, g2featurelist, seq_annotations, all_annotations, ignore_exp, min_appearances, num_results_needed)
            logger.info('found %d positive' % len(term_idx))
            keep_min_exps = np.concatenate((keep_min_exps, npos[term_idx]))
            num_total_exps = np.concatenate((num_total_exps, cnum_total_exps))
            num_enriched_exps = np.concatenate((num_enriched_exps, cnum_enriched_exps))

            # keep only the terms that are enriched in at least min_appearances experiments
            term_list = term_list[keep_min_exps]
            orig_term_list = orig_term_list[keep_min_exps]
            odif = odif[keep_min_exps]
            pvals = pvals[keep_min_exps]
            all_feature_array = all_feature_array[:, keep_min_exps]
            # and re-sorst since we went from edges to middle
            si = np.argsort(odif)
            term_list = term_list[si]
            orig_term_list = orig_term_list[si]
            odif = odif[si]
            pvals = pvals[si]
            all_feature_array = all_feature_array[:, si]
            num_enriched_exps = np.array(num_enriched_exps)[si]
            num_total_exps = np.array(num_total_exps)[si]

            # add the number of enriched experiments to each term string
            # for idx, cterm in enumerate(term_list):
            #     term_list[idx] = term_list[idx] + ' [%d/%d]' % (num_enriched_exps[idx], num_total_exps[idx])

        # for rank data, normalize the effect size to be in the [-1:1] range (0 for random, -1 / 1 for fully ordered)
        if transform_type == 'rankdata':
            n_g1 = len(g1_features)
            n_g2 = len(g2_features)
            odif = odif / ((((n_g1 + 1) / 2) + n_g2) - ((n_g2 + 1) / 2))

        # create dataframe for resulting terms.
        res = pd.DataFrame({'term': term_list, 'odif': odif, 'pvals': pvals, 'num_enriched_exps': num_enriched_exps, 'num_total_exps': num_total_exps}, index=term_list)
        # res = pd.DataFrame({'term': term_list, 'odif': odif, 'pvals': pvals, 'qvals': qvals, 'num_enriched_exps': num_enriched_exps, 'num_total_exps': num_total_exps}, index=term_list)
        features = list(g1_features)
        features.extend(g2_features)
        features = pd.DataFrame({'sequence': features, 'group': labels}, index=features)
        return res, all_feature_array, features

    def _get_term_enriched_annotations(self, g1_features, g2_features, all_annotations, seq_annotations, ignore_exp=None, **kwargs):
        '''DEPRACATED!!!!!
        Get a list of enriched annotations for each term
        Parameters
        ----------

        Returns
        -------
        enriched_annotations: dict of {annotationid(int): enriched(bool)}
        '''
        # get a dict with keys all experiment ids from the annotations
        all_exps = self._get_exp_annotations(all_annotations)

        # the number of experiments where each term was enriched (key is the term, value is the count of experiments)
        exp_term_count = defaultdict(int)
        for cexpid in all_exps.keys():
            # create the new ignore list for all experiments except the one we are processing
            exclude_exps = set(all_exps.keys())
            if cexpid in exclude_exps:
                exclude_exps.remove(cexpid)
            if ignore_exp is not None:
                exclude_exps = exclude_exps.union(set(ignore_exp))
            sterms, ftscore, sfeatures = self.term_enrichment(g1_features, g2_features, all_annotations, seq_annotations, ignore_exp=list(exclude_exps), **kwargs)
            for cterm in sterms.index.values:
                exp_term_count[cterm] += 1
        return exp_term_count

    def get_ontology_terms(self, min_term_id=None, ontologyid=None):
        '''Get all the terms in an ontology, starting with termid min_term_id
        By default it gets the terms for the dbbact ontology (i.e. all terms not in envo/gaz/etc.)

        Parameters
        ----------
        min_term_id: int or None, optional
            The minimal dbbact term id to get in the results (to make it faster)
        ontologyid: int or None, optional
            The ontologyid to get the results for (1 is dbbact ontology)
            None to get all ontologies

        Returns
        -------
        ontology:
            dict of {term(str): dbbact_id(int)}
                term: the ontology term (i.e. 'feces')
                dbbact_id: the internal dbbact id for the term
        ontology_term_ids:
            dict of {dbbact_id(int): ontology_term_id (str)}
                dbbact_id: the internal dbbact id for the term
                ontology_term_id: the id of the term in the ontology (i.e. 'UBERON:00004')
        '''
        rdata = {'min_term_id': min_term_id, 'ontologyid': ontologyid}
        res = self._get('ontology/get_all_terms', rdata)
        if res.status_code != 200:
            msg = 'Failed to get all ontology terms list. error %s' % res.content
            logger.warn(msg)
            return msg
        terms = res.json()
        logger.debug('got %d new terms' % len(terms['ontology']))
        return terms['ontology'], terms['ontology_term_ids']

    def _get_term_features(self, features, feature_terms):
        '''Get dict of number of appearances in each sequence keyed by term

        Parameters
        ----------
        features : list of str
            A list of DNA sequences
        feature_terms : dict of {feature: list of tuples of (term, amount)}
            The terms associated with each feature in exp
            feature (key) : str the feature (out of exp) to which the terms relate
            feature_terms (value) : list of tuples of (str or int the terms associated with this feature, count)

        Returns
        -------
        numpy array of T (terms) * F (features)
            total counts of each term (row) in each feature (column)
        list of str
            list of the terms corresponding to the numpy array rows
        '''
        # get all terms. store the index position for each term
        terms = {}
        cpos = 0
        for cfeature, ctermlist in sorted(feature_terms.items()):
            for cterm, ccount in ctermlist:
                if cterm not in terms:
                    terms[cterm] = cpos
                    cpos += 1

        res = np.zeros((len(terms), len(features)))
        for idx, cfeature in enumerate(features):
            for ctermlist in feature_terms[cfeature]:
                cterm = ctermlist[0]
                cscore = ctermlist[1]
                res[terms[cterm], idx] = cscore
        term_list = sorted(terms, key=terms.get)
        return res, term_list

    def get_term_annotations(self, term, features, feature_annotations, annotations):
        '''
        Get a matrix of annotations involving term for each feature

        Parameters
        ----------
        term: str
            the term to get annotations for
        features: list of str
            the features to get the annotations for
        feature_annotations: dict {str: list of int}
            per feature (key) list of annotation ids
        annotations : dict of {int: dict}
            key is annotationID, dict is the annotation details (see XXX)

        Returns
        -------
        numpy 2d array
            with annotatios (containing the term) (columns) x features (rows)
            value at each position is >0 if this feature has this annotation or 0 if not
            the number indicates the annotation type for the detail:
            4 - common
            8 - dominant
            5 - all (in diffexp annotations)
            2 - low (in diffexp annotations / negative correlation)
            16 - high (in diffexp annotations / positive correlation)
            32 - other annotation type
        pandas dataframe
            one row per annotation
        pandas dataframe
            one row per feature
        '''
        term_annotations = pd.DataFrame(columns=['annotation', 'annotationid', 'expid', 'annotation_type', 'detail_type'])
        for cannotation in annotations.values():
            details = cannotation['details']
            for cdetail in details:
                cterm = cdetail[1]
                ctype = cdetail[0]
                if cterm == term:
                    cdat = {'annotationid': cannotation['annotationid'], 'annotation_type': cannotation['annotationtype'], 'expid': cannotation['expid'], 'detail_type': ctype}
                    cdat['annotation'] = self.get_annotation_string(cannotation)
                    cdf = pd.DataFrame([cdat], index=[cannotation['annotationid']])
                    term_annotations = term_annotations.append(cdf)
                    break
        logger.debug('found %d annotations with term %s' % (len(term_annotations), term))
        term_mat = np.zeros([len(features), len(term_annotations)])
        for idx, cfeature in enumerate(features):
            for cannotation_id in feature_annotations[cfeature]:
                # cannotation_id = cannotation_data[1]
                if cannotation_id in term_annotations.index:
                    annotation_pos = term_annotations.index.get_loc(cannotation_id)
                    if term_annotations['annotation_type'][cannotation_id] == 'common':
                        # score = 4
                        score = 1
                    elif term_annotations['annotation_type'][cannotation_id] == 'dominant':
                        # score = 8
                        score = 2
                    elif term_annotations['annotation_type'][cannotation_id] == 'diffexp':
                        if term_annotations['detail_type'][cannotation_id] == 'all':
                            # score = 5
                            score = 3
                        elif term_annotations['detail_type'][cannotation_id] == 'low':
                            # score = 2
                            score = 4
                        elif term_annotations['detail_type'][cannotation_id] == 'high':
                            # score = 16
                            score = 5
                    else:
                        # score = 32
                        score = 6
                    term_mat[idx, annotation_pos] = score
        term_annotations.index = term_annotations.index.map(str)
        return term_mat, term_annotations, pd.DataFrame({'sequence': features}, index=features)

    def term_neighbors(self, term):
        '''Get the closest terms to a given term

        Based on the fraction of shared bacteria between the terms
        For database metanalysis.

        Parameters
        ----------
        term : str
            the term to get the neighbors for

        Returns
        -------

        '''
        # get the bacteria associated with the term
        term_features = self.get_db_term_features(term)

        # now go over all bacteria and find other terms associated with each
        terms_to_check = set()
        for cfeature in term_features.keys():
            annotations, term_info, taxonomy = self.get_seq_annotations(cfeature)
            for cannotation in annotations:
                for cdetail in cannotation['details']:
                    terms_to_check.add(cdetail[1])
        logger.debug('need to check %d details' % len(terms_to_check))
        term_dist = {}
        for cterm in terms_to_check:
            logger.debug(cterm)
            try:
                cterm_f = self.get_db_term_features(cterm)
                common = len(set(term_features.keys()).intersection(set(cterm_f.keys())))
                cscore = common / (len(term_features) + len(cterm_f))
                term_dist[cterm] = cscore
                logger.debug('oterm %d, cterm %d, common %d, score %f' % (len(term_features), len(cterm_f), common, cscore))
            except:
                logger.warn('failed for term %s' % cterm)
        return term_dist

    def get_db_term_features(self, terms, ignore_exp=(), max_id=None):
        '''Get all the features associated with a term in the database

        Parameters
        ----------
        terms : str or list of str
            The term to look for in the database. if list of str, look for annotations containing all terms
        ignore_exp: list of int, optional
            list of experiment ids to ignore when looking for the features

        Returns
        -------
        dict of {feature: num} {str: int}. Key is feature (sequence), value is number of annotations containing this feature for this term.
        '''
        rdata = {}
        # convert to list if single term
        if isinstance(terms, str):
            terms = [terms]

        rdata['term'] = terms
        res = self._get('ontology/get_annotations', rdata, param_json=False)
        if res.status_code != 200:
            logger.warn('error getting annotations for term %s' % terms)
            return []
        res = res.json()
        annotations = res.get('annotations')
        logger.info('found %d annotations with the term %s' % (len(annotations), terms))

        # get features only for annotations which don't contain "lower in " for one of the terms
        terms_set = set(terms)
        feature_num = defaultdict(int)
        num_ignored_annotations = 0
        for cannotation in annotations:
            # if annotation is from an experiment we are ignoring, skip it
            if cannotation['expid'] in ignore_exp:
                continue
            if max_id is not None:
                if cannotation['annotationid'] > max_id:
                    num_ignored_annotations += 1
                    continue
            foundit = set()
            for cdetail in cannotation['details']:
                if cdetail[1] not in terms_set:
                    continue
                if cdetail[0] == 'low':
                    logger.info('annotation %d is low' % cannotation['annotationid'])
                    continue
                foundit.add(cdetail[1])
                if len(foundit) >= len(terms):
                    break
            if len(foundit) < len(terms):
                continue
            # get the sequences for the annotation
            seqs = self.get_annotation_sequences(cannotation['annotationid'])
            for cseq in seqs:
                feature_num[cseq] += 1
        if max_id is not None:
            logger.info('ignored %d annotations' % num_ignored_annotations)
        return feature_num

    def get_sequences_ids(self, sequences, no_shorter=False, no_longer=False, use_sequence_translator=True):
        '''Get the dbbact IDs for a list of ACGT sequences

        Parameters
        ----------
        sequences: list of str
            the ACGT sequences to get ids for
        no_shorter: bool, optional
            True to not get sequences shorter than the query sequence (but match 100%)
        no_longer: bool, optional
            True to not get sequences that are longer than the query sequence (but match 100%)
        use_sequence_translator: bool, optional
            True to use also sequences from other regions (using SILVA sequence translator)
        Returns
        -------
        list of list of int
            List of dbbact sequence ids matching each query sequence
        '''
        rdata = {}
        rdata['no_shorter'] = no_shorter
        rdata['no_longer'] = no_longer
        rdata['use_sequence_translator'] = use_sequence_translator
        rdata['sequences'] = sequences
        res = self._get('sequences/getid_list', rdata)
        if res.status_code != 200:
            logger.warn('error getting sequenceids for %d sequences' % len(sequences))
            return []
        res = res.json()
        seqids = res.get('seqIds')
        return seqids

    def get_annotation_sequences(self, annotationid):
        '''Get all sequence ids associated with an annotation

        Parameters
        ----------
        annotationid : int
            the annotationid to get the sequences for

        Returns
        -------
        list of str. The sequences associated with the annotation
        '''
        rdata = {}
        rdata['annotationid'] = annotationid
        res = self._get('annotations/get_sequences', rdata)
        if res.status_code != 200:
            logger.warn('error getting sequences for annotationid %s' % annotationid)
            return []
        res = res.json()
        seqids = res.get('seqids')
        return seqids
        logger.debug('found %d sequences for annotationid %d' % (len(seqids), annotationid))
        seqs = self.get_seq_id_info(seqids)
        return seqs

    def get_seq_id_info(self, seqids):
        '''Get sequences for sequenceids

        Parameters
        ----------
        ids : list of int
            list of the sequenceids (from get_annotation_sequences)

        Returns
        -------
        list of str
            sequence of each id
        '''
        rdata = {}
        rdata['seqids'] = seqids
        res = self._get('sequences/get_info', rdata)
        if res.status_code != 200:
            logger.warn('error getting sequences for seqids %s' % seqids)
            return []
        res = res.json()
        seqdat = res.get('sequences')
        seqs = [x['seq'] for x in seqdat]
        return seqs

    def all_term_neighbors(self, limit=[]):
        '''Get all the teature_ids associated with each term in the database

        Parameters
        ----------
        limit : list of str (optional)
            Only use annotations that contain all these terms

        Returns
        -------
        dict of {str: {int: int}}
            dict of {term : {feature_id : count}}
            key is term, value is dict where key is the feature_id (int from dbBact), value is the number of annotations (containing the term) with the feature
        '''
        logger.info('getting all annotations list from dbBact')
        res = self._get('annotations/get_all_annotations', {}, param_json=True)
        res = res.json()
        annotations = res['annotations']
        logger.info('got %d annotations' % len(annotations))

        term_features = defaultdict(lambda: defaultdict(int))
        term_experiments = defaultdict(set)
        for cannotation in annotations:
            found_limits = True
            for climit_term in limit:
                climit_found = False
                for cdetail in cannotation['details']:
                    if cdetail[1] == climit_term:
                        climit_found = True
                        break
                if not climit_found:
                    found_limits = False
            if not found_limits:
                continue
            annotation_seqs = self.get_annotation_sequences(cannotation['annotationid'])
            for cdetail in cannotation['details']:
                if cdetail[0] == 'low':
                    cterm = '-' + cdetail[1]
                else:
                    cterm = cdetail[1]
                for cseq in annotation_seqs:
                    term_features[cterm][cseq] += 1
                    term_experiments[cterm].add(cannotation['expid'])
        return term_features, term_experiments

    def get_term_neighbors(self, term_features, term):
        term_dist = {}
        term_f = term_features[term]
        for cterm, cterm_f in term_features.items():
            common = len(set(term_f.keys()).intersection(set(cterm_f.keys())))
            cscore = common / (len(term_f) + len(cterm_f))
            term_dist[cterm] = cscore
        res = sorted(term_dist.items(), key=lambda x: x[1], reverse=True)
        return res

    def get_term_term_stats(self, term_features, term1, term2):
        term1_f = term_features[term1]
        term2_f = term_features[term2]
        common = len(set(term1_f.keys()).intersection(set(term2_f.keys())))
        cscore = common / (len(term1_f) + len(term2_f))
        logger.debug('term1 features %d, term2 features %d, common %d, score %f' % (len(term1_f), len(term2_f), common, cscore))

    def count_enriched_exps(self, term, g1features, g2features, seq_annotations, annotations, ignore_exp=None, **kwargs):
        '''Get experiments with enriched term annotations for a given term.

        Parameters
        ----------
        term: str
            the term to look for enriched annotations
        g1features, g2features: list of str
            list of bacterial sequences  for the two groups to test enrichment between
        seq_annotations: dict of {feature(str): list of annotationids (int)}
        annotations: dict of {annotationsid(int): }
        ignore_exp: list or None, optional
            experiment ids to ignore when counting
        **kwargs:
            parameters for dsfdr. can include:
            method=method, transform_type=transform_type, alpha=alpha, numperm=numperm, fdr_method=fdr_method

        Returns
        -------
        enriched_experiments: dict of {expid(int): num of enriched annotations(int)}
        enriched_annotations: list of [annotationids(int)]
        total_exp_annotations: dict of {expid(int): num of annotations with term(int)}
        '''
        # get rid of lower in, since it means higer in the other group????
        # TODO: FIX THIS!!!!!
        if term[0] == '-':
            lower = True
            term = term[1:]
        else:
            lower = False

        if ignore_exp is None:
            ignore_exp = []
        ignore_exp = set(ignore_exp)

        # get all annotations matrix for the term
        all_seqs = list(g1features) + list(g2features)
        tmat, tanno, tseqs = self.get_term_annotations(term, all_seqs, seq_annotations, annotations)

        labels = np.zeros(len(all_seqs))
        labels[:len(g1features)] = 1

        # find which terms are significantntly enriched in either of the two feature groups
        keep, odif, pvals, qvals = dsfdr(tmat.T, labels, **kwargs)

        # count enriched experiments
        enriched_experiments = defaultdict(int)
        total_exp_annotations = defaultdict(int)
        enriched_annotations = []
        for cpos, ckeep in enumerate(keep):
            cexp = tanno.iloc[cpos]['expid']
            if cexp in ignore_exp:
                continue
            canno = tanno.iloc[cpos]['annotationid']
            # if we look at LOWER IN annotations, only count annotations where it really appears LOWER IN
            if lower:
                if annotations[canno]['annotationtype'] != 'diffexp':
                    continue
                islower = False
                for cdet in annotations[canno]['details']:
                    if cdet[0] == 'low':
                        if cdet[1] == term:
                            islower = True
                            break
                if not islower:
                    continue
            total_exp_annotations[cexp] += 1
            if ckeep:
                enriched_experiments[cexp] += 1
                enriched_annotations.append(canno)

        keep = np.where(keep)[0]
        if len(keep) == 0:
            logger.debug('no enriched annotations found for term %s' % term)
        else:
            logger.debug('Found %d enriched annotations, %d enriched experiments for term %s' % (len(keep), len(enriched_experiments), term))

        return enriched_experiments, enriched_annotations, total_exp_annotations

    def get_sequences_primer(self, sequences):
        '''Get the primer name for the sequences

        Parameters
        ----------
        sequences: list of str
            the sequences to get the primer region for (from dbbact))

        Returns
        -------
        primer_name: str
            name of the primer region (i.e. 'v4')
        '''
        rdata = {}
        rdata['sequences'] = sequences
        res = self._get('sequences/guess_region', rdata)
        if res.status_code != 200:
            logger.warn('error getting sequences for sequences: %s' % res.content)
            return None
        res = res.json()
        if res['region'] == 'na':
            # force user to select
            return None
        return res['region']

    def remove_features_from_annotation(self, features, data):
        '''remove a set of features from a given annotation
        Tries to remove the features from the annotation using the username/password from the config file.
        A user can only change annotations he created or anonymous annotations.
        Otherwise, an error will be returned

        Parameters
        ----------
        features: list of str
            the features to remove from the annotation
        data : dict
            The annotationdetails to remove the features from (using the 'annotationid' key)

        Returns
        -------
        str:
            empty if ok, non empty string if error encountered
        '''
        if 'annotationid' not in data:
            return 'No annotationID for selected annotation'
        if len(features) == 0:
            return 'No sequences selected to delete from annotation'
        annotationid = data['annotationid']
        rdata = {}
        rdata['annotationid'] = annotationid
        rdata['sequences'] = features
        res = self._post('annotations/delete_sequences_from_annotation', rdata)
        if res.status_code != 200:
            msg = 'Delete sequeces from annotations failed. error %s' % res.content
            logger.warn(msg)
            return msg
        logger.info('%d sequences removed from annotation %d deleted' % (len(features), annotationid))
        return ''
