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
from logging import getLogger, NOTSET, basicConfig
from logging.config import fileConfig
from pkg_resources import resource_filename

import numpy as np
import pandas as pd

from .dsfdr import dsfdr

from . import __version_numeric__

logger = getLogger(__name__)

try:
    # get the logger config file location
    log = resource_filename(__name__, 'log.cfg')

    # set the logger output according to log.cfg
    # setting False allows other logger to print log.
    fileConfig(log, disable_existing_loggers=False)

    # set the log level to same value as calour log level
    clog = getLogger('calour')
    calour_log_level = clog.getEffectiveLevel()
    if calour_log_level != NOTSET:
        logger.setLevel(calour_log_level)
except:
    print('FAILED loading logging config file %s' % log)
    basicConfig(format='%(levelname)s:%(message)s')


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

    def version(self):
        return __version_numeric__

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
            logger.warn('REST error %s enountered when accessing dbBact %s' % (res.reason, api))
        return res

    def get_seq_annotations(self, sequence):
        '''Get the annotations for a sequence

        Parameters
        ----------
        sequence : str
            The DNA sequence to get the annotations for

        Returns
        -------
        annotations : list of list of (annotation dict,list of [Type,Value] of annotation details)
            See dbBact sequences/get_annotations REST API documentation
        term_info : dict of {term: info}
            where key (str) is the ontology term, info is a dict of details containing:
                'total_annotations' : total annotations having this term in the database
                'total_sequences' : number of annotations with this term for the sequence
        '''
        rdata = {}
        rdata['sequence'] = sequence
        res = self._get('sequences/get_annotations', rdata)
        if res.status_code != 200:
            logger.warn('error getting annotations for sequence %s' % sequence)
            return []
        res = res.json()
        annotations = res.get('annotations')
        term_info = res.get('term_info')
        logger.debug('Found %d annotations for sequence %s' % (len(annotations), sequence))
        return annotations, term_info

    def get_annotation_string(self, cann):
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

    def get_seq_annotation_strings(self, sequence):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        sequence : str
            the DNA sequence to query the annotation strings about

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'annotationid' : int, the annotation id in the database
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short summary of the annotation
        '''
        shortdesc = []
        annotations, term_info = self.get_seq_annotations(sequence)
        if len(term_info) > 0:
            terms = []
            for cterm, cinfo in term_info.items():
                terms.append([cterm, cinfo.get('total_sequences'), cinfo.get('total_annotations')])
            terms = sorted(terms, key=lambda x: x[1])
            terms = sorted(terms, key=lambda x: -x[2])
            summary = 'most: '
            most_terms = self.get_common_annotation_term(annotations, term_info, num_common=5)
            for cmterm in most_terms:
                summary += cmterm + ', '
            shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, summary))
            summary = 'special: '
            terms = sorted(terms, key=lambda x: x[2])
            for cterm in terms[:min(4, len(terms))]:
                summary += '%s, ' % cterm[0]
            shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, summary))
        for cann in annotations:
            cdesc = self.get_annotation_string(cann)
            shortdesc.append((cann, cdesc))
        return shortdesc

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
            elif annotation_type == 'highfreq':
                ascore = 2
            elif annotation_type == 'other':
                ascore = 0.5
            elif annotation_type == 'contamination':
                ascore = 1
                details = [('all', 'contamination')]
            elif annotation_type == 'diffexp':
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
            The number of common terms to return

        Returns
        -------
        list of str
            list of the num_common most common terms
        '''
        term_count = defaultdict(int)
        for cannotation in annotations:
            if cannotation['annotationtype'] == 'common':
                for cdesc in cannotation['details']:
                    term_count[cdesc[1]] += 1
                continue
            if cannotation['annotationtype'] == 'highfreq':
                for cdesc in cannotation['details']:
                    term_count[cdesc[1]] += 2
                continue
            if cannotation['annotationtype'] == 'other':
                for cdesc in cannotation['details']:
                    term_count[cdesc[1]] += 0.5
                continue
            if cannotation['annotationtype'] == 'contamination':
                term_count['contamination'] += 1
                continue
            if cannotation['annotationtype'] == 'diffexp':
                for cdesc in cannotation['details']:
                    if cdesc[0] == 'all':
                        term_count[cdesc[1]] += 1
                        continue
                    if cdesc[0] == 'high':
                        term_count[cdesc[1]] += 2
                        continue
                    if cdesc[0] == 'low':
                        term_count[cdesc[1]] -= 2
                        continue
                continue
            if cannotation['annotationtype'] == 'other':
                continue
            logger.warn('unknown annotation type %s encountered' % cannotation['annotationtype'])
        # now correct for common / uncommon terms
        for k, v in term_count.items():
            if k not in term_info:
                logger.debug('term %s not in term_info' % k)
                continue
            delta = 1 / term_info[k]['total_sequences']
            # flip and add '-' to term if negative
            if v < 0:
                k = '-' + k
                v = -v
            v += delta

        sorted_terms = sorted(term_count, key=term_count.get, reverse=True)
        return(sorted_terms[:num_common])

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
            logger.debug('error adding experiment. msg: %s' % res.content)
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
            'HIGHFREQ' : the sequences are found in mean high frequency (>1%) in the samples
            'PATHOGEN' : the sequences are known pathogens
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
        annotationid : int or None
            the AnnotationID (in Annotations table) of the new annotation, or None if not added
        """
        logger.debug('addannotation - %d sequences' % len(sequences))
        if len(sequences) == 0:
            logger.warn('No sequences to annotate!')
            return None
        if len(annotations) == 0:
            logger.warn('No annotations to add!')
            return None

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
            return newid
        logger.warn('problem adding annotations for experiment id %d' % expid)
        logger.debug(res.content)
        return None

    def get_seq_list_fast_annotations(self, sequences):
        '''Get annotations for all sequences in list using compact format and with parent ontology terms

        Params
        ------
        sequences : list of str
            The list of DNA sequences to get the annotations for

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
        '''
        logger.info('Getting dbBact annotations for %d sequences, please wait...' % len(sequences))
        rdata = {}
        rdata['sequences'] = list(sequences)
        res = self._get('sequences/get_fast_annotations', rdata)
        if res.status_code != 200:
            logger.warning('error getting fast annotations for sequence list. got status code %s' % res.status_code)
            raise ValueError('error getting fast annotations for sequence list. got status code %s' % res.status_code)
        res = res.json()

        sequence_terms = {}
        sequence_annotations = {}
        for cseq in sequences:
            sequence_terms[cseq] = []
            sequence_annotations[cseq] = []
        for cseqannotation in res['seqannotations']:
            cpos = cseqannotation[0]
            # need str since json dict is always string
            cseq = sequences[cpos]
            sequence_annotations[cseq].extend(cseqannotation[1])
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

        total_annotations = 0
        for cseq_annotations in sequence_annotations.values():
            total_annotations += len(cseq_annotations)
        logger.info('Got %d annotations' % total_annotations)

        return sequence_terms, sequence_annotations, res['annotations'], res['term_info']

    def get_feature_terms(self, features, exp=None, term_type=None, ignore_exp=None):
        '''Get list of terms per feature

        Parameters
        ----------
        features : list of str
            the features to get the terms for
        exp : calour.Experiment (optional)
            not None to store results in the exp (to save time for multiple queries)
        term_type : str or None (optional)
            The type of terms to return. optional values:
            None to use default ('seqstring')
            'annotation': the annotation string for each annotation (i.e. 'higher in fish compared to dogs...')
            'terms': the ontology terms without parent terms (with '-' attached to negative freq. terms)
            'parentterms': the ontology terms including all parent terms (with '-' attached to negative freq. terms)
            'contamination': get if bacteria is contaminant or not
        ignore_exp : list of int or None (optional)
            the list of experimentids to ignore (don't take info from them)

        Returns
        -------
        feature_terms : dict of list of str/int
            key is the feature, list contains all terms associated with the feature
        '''
        if term_type is None:
            term_type = 'terms'
        if exp is not None:
            if '__dbbact_sequence_terms' not in exp.exp_metadata:
                # if annotations not yet in experiment - add them
                self.add_all_annotations_to_exp(exp)
            # and filter only the ones relevant for features
            sequence_terms = exp.exp_metadata['__dbbact_sequence_terms']
            sequence_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
            annotations = exp.exp_metadata['__dbbact_annotations']
        else:
            sequence_terms, sequence_annotations, annotations, term_info = self.get_seq_list_fast_annotations(features)
        new_annotations = {}
        if term_type == 'annotation':
            for cseq, annotations_list in sequence_annotations.items():
                if cseq not in features:
                    continue
                newdesc = []
                for cannotation in annotations_list:
                    if ignore_exp is not None:
                        annotationexp = annotations[cannotation]['expid']
                        if annotationexp in ignore_exp:
                            continue
                    cdesc = self.get_annotation_string(annotations[cannotation])
                    newdesc.append(cdesc)
                new_annotations[cseq] = newdesc
        elif term_type == 'terms':
            for cseq, annotations_list in sequence_annotations.items():
                if cseq not in features:
                    continue
                newdesc = []
                for cannotation in annotations_list:
                    if ignore_exp is not None:
                        annotationexp = annotations[cannotation]['expid']
                        if annotationexp in ignore_exp:
                            continue
                    for cdesc in annotations[cannotation]['details']:
                        if cdesc[0] == 'all' or cdesc[0] == 'high':
                            cterm = cdesc[1]
                        else:
                            cterm = '-' + cdesc[1]
                        newdesc.append(cterm)
                new_annotations[cseq] = newdesc
        elif term_type == 'parentterms':
            for cseq, term_list in sequence_terms.items():
                if cseq not in features:
                    continue
                term_list = [x for x in term_list if x != 'na']
                new_annotations[cseq] = term_list
        elif term_type == 'contamination':
            for cseq, annotations_list in sequence_annotations.items():
                if cseq not in features:
                    continue
                newdesc = []
                is_contamination = 0
                for cannotation in annotations_list:
                    if ignore_exp is not None:
                        annotationexp = annotations[cannotation]['expid']
                        if annotationexp in ignore_exp:
                            continue
                    if annotations[cannotation]['annotationtype'] == 'contamination':
                        is_contamination += 1
                if is_contamination > 0:
                    new_annotations[cseq] = ['contamination']
                else:
                    new_annotations[cseq] = []
        else:
            raise ValueError('term_type %s not supported in get_feature_terms. Possible values are "annotation","terms","parentterms"' % term_type)
        return new_annotations
        # return sequence_annotations
        # return sequence_terms

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

    def term_enrichment(self, g1_features, g2_features, all_annotations, seq_annotations, term_type='term', ignore_exp=None, min_appearances=3, fdr_method='dsfdr', score_method='all_mean', random_seed=None, use_term_pairs=False, alpha=0.1, method='meandiff', transform_type='rankdata', numperm=1000):
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
        term_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
             'parentterm' - ontology terms including parent terms associated with each feature.
             'annotation' - the full annotation strings associated with each feature
             'combined' - combine 'term' and 'annotation'
        ignore_exp: list of int or None or True(optional)
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment
            None (default) to use annotations from all experiments including the current one
        min_appearances : int (optional)
            The minimal number of times a term appears in order to include in output list.
        fdr_method : str (optional)
            The FDR method used for detecting enriched terms (permutation test). options are:
            'dsfdr' (default): use the discrete FDR correction
            'bhfdr': use Benjamini Hochbert FDR correction
        transform_type: str, optional
            the data trasformation to use before calculating the dsFDR. For options see dsfdr doc
        score_method : str (optional)
            The score method used for calculating the term score. Options are:
            'all_mean' (default): mean over each experiment of all annotations containing the term
            'sum' : sum of all annotations (experiment not taken into account)
            'card_mean': use a null model keeping the number of annotations per each bacteria
        random_seed: int or None
            int to specify the random seed for numpy.random.
        use_term_pairs: bool, optional
            True to also test enrichment in pairs of terms (i.e. homo sapiens+feces, etc.)

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
        if random_seed is not None:
            np.random.seed(random_seed)

        g1_features = np.array(g1_features)
        g2_features = np.array(g2_features)
        exp_features = np.hstack([g1_features, g2_features])

        if term_type == 'term':
            feature_terms = self._get_all_term_counts(exp_features, seq_annotations, all_annotations, ignore_exp=ignore_exp, score_method=score_method, use_term_pairs=use_term_pairs)
        elif term_type == 'parentterm':
            pass
        elif term_type == 'annotation':
            feature_terms = self._get_all_annotation_string_counts(exp_features, all_annotations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
        elif term_type == 'combined':
            feature_terms = self._get_all_term_counts(exp_features, seq_annotations, all_annotations, ignore_exp=ignore_exp, score_method=score_method, use_term_pairs=use_term_pairs)
            feature_annotations = self._get_all_annotation_string_counts(exp_features, all_annotations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
            for cfeature, cvals in feature_annotations.items():
                if cfeature not in feature_terms:
                    feature_terms[cfeature] = []
                feature_terms[cfeature].extend(cvals)
        else:
            raise ValueError('term_type %s not supported for dbbact. possible values are: "term", "parentterm", "annotation"')

        # arrays of terms (rows) x bacteria (cols)
        feature_array, term_list = self._get_term_features(g1_features, feature_terms)
        bg_array, term_list = self._get_term_features(g2_features, feature_terms)

        # and the array of terms (rows) x all bacteria (in both groups) (cols)
        all_feature_array = np.hstack([feature_array, bg_array])

        # remove non-informative terms (present in not enough bacteria)
        non_informative = np.sum(all_feature_array > 0, 1) < min_appearances
        all_feature_array = np.delete(all_feature_array, np.where(non_informative)[0], axis=0)
        term_list = [x for idx, x in enumerate(term_list) if not non_informative[idx]]

        # remove terms present in one experiment
        term_exps = defaultdict(set)
        num_removed = 0
        for cannotation in all_annotations.values():
            for cdetail in cannotation['details']:
                cterm = cdetail[1]
                term_exps[cterm].add(cannotation['expid'])
        new_term_list = []
        for cterm in term_list:
            if cterm[0] == '-':
                ccterm = cterm[1:]
            else:
                ccterm = cterm
            if len(term_exps[ccterm]) == 1:
                cterm = '**%s**%s' % (list(term_exps[ccterm])[0], cterm)
                num_removed += 1
            new_term_list.append(cterm)
        term_list = new_term_list
        logger.info('removed %d terms' % num_removed)

        labels = np.zeros(all_feature_array.shape[1])
        labels[:feature_array.shape[1]] = 1

        keep, odif, pvals = dsfdr(all_feature_array, labels, method=method, transform_type=transform_type, alpha=alpha, numperm=numperm, fdr_method=fdr_method)
        keep = np.where(keep)[0]
        if len(keep) == 0:
            logger.info('no enriched terms found')
        term_list = np.array(term_list)[keep]
        all_feature_array = all_feature_array[keep, :].T
        all_feature_array = all_feature_array * 100
        odif = odif[keep]
        pvals = pvals[keep]
        si = np.argsort(odif)
        term_list = term_list[si]
        odif = odif[si]
        pvals = pvals[si]
        all_feature_array = all_feature_array[:, si]
        res = pd.DataFrame({'term': term_list, 'odif': odif, 'pvals': pvals}, index=term_list)
        features = list(g1_features)
        features.extend(g2_features)
        # tmat, tanno, tseqs = self.get_term_annotations("crohn's disease", features, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        # return tanno, tmat, tseqs
        features = pd.DataFrame({'sequence': features, 'group': labels}, index=features)
        return res, all_feature_array, features

    def _get_term_features2(self, features, feature_terms):
        '''Get dict of number of appearances in each sequence keyed by term.
        This function inflates each bacteria to the total number of annotations it has.
        So it can be used for the randomizing annotation count null hypothesis model

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
        # get all terms
        terms = {}
        cpos = 0
        for cfeature, ctermlist in feature_terms.items():
            for cterm, ccount in ctermlist:
                if cterm not in terms:
                    terms[cterm] = cpos
                    cpos += 1

        tot_features_inflated = 0
        feature_pos = {}
        for cfeature in features:
            ctermlist = feature_terms[cfeature]
            feature_pos[cfeature] = tot_features_inflated
            tot_features_inflated += len(ctermlist)

        res = np.zeros([len(terms), tot_features_inflated])

        for cfeature in features:
            for cterm, ctermcount in feature_terms[cfeature]:
                res[terms[cterm], feature_pos[cfeature]] += ctermcount
        term_list = sorted(terms, key=terms.get)
        return res, term_list

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
        for cfeature, ctermlist in feature_terms.items():
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
            value at each position is 1 if this feature has this annotation or 0 if not
        pandas dataframe
            one row per annotation
        pandas dataframe
            one row per feature
        '''
        term_annotations = pd.DataFrame()
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
        logger.info('found %d annotations with term' % len(term_annotations))
        term_mat = np.zeros([len(features), len(term_annotations)])
        for idx, cfeature in enumerate(features):
            for cannotation_id in feature_annotations[cfeature]:
                # cannotation_id = cannotation_data[1]
                if cannotation_id in term_annotations.index:
                    annotation_pos = term_annotations.index.get_loc(cannotation_id)
                    if term_annotations['annotation_type'][cannotation_id] == 'common':
                        score = 4
                    elif term_annotations['annotation_type'][cannotation_id] == 'highfreq':
                        score = 8
                    elif term_annotations['annotation_type'][cannotation_id] == 'diffexp':
                        if term_annotations['detail_type'][cannotation_id] == 'all':
                            score = 5
                        elif term_annotations['detail_type'][cannotation_id] == 'low':
                            score = 2
                        elif term_annotations['detail_type'][cannotation_id] == 'high':
                            score = 16
                    else:
                        score = 32
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
            annotations, term_info = self.get_seq_annotations(cfeature)
            for cannotation in annotations:
                for cdetail in cannotation['details']:
                    terms_to_check.add(cdetail[1])
        print('need to check %d details' % len(terms_to_check))
        term_dist = {}
        for cterm in terms_to_check:
            print(cterm)
            try:
                cterm_f = self.get_db_term_features(cterm)
                common = len(set(term_features.keys()).intersection(set(cterm_f.keys())))
                cscore = common / (len(term_features) + len(cterm_f))
                term_dist[cterm] = cscore
                print('oterm %d, cterm %d, common %d, score %f' % (len(term_features), len(cterm_f), common, cscore))
            except:
                logger.warn('failed for term %s' % cterm)
        return term_dist

    def get_db_term_features(self, term):
        '''Get all the features associated with a term in the database

        Parameters
        ----------
        term : str
            The term to look for in the database

        Returns
        -------
        dict of {feature: num} {str: int}. Key is feature (sequence), value is number of annotations containing this feature for this term.
        '''
        rdata = {}
        rdata['term'] = term
        res = self._get('ontology/get_annotations', rdata, param_json=False)
        if res.status_code != 200:
            logger.warn('error getting annotations for term %s' % term)
            return []
        res = res.json()
        annotations = res.get('annotations')
        logger.info('found %d annotations with the term %s' % (len(annotations), term))

        feature_num = defaultdict(int)
        for cannotation in annotations:
            foundit = False
            for cdetail in cannotation['details']:
                if cdetail[1] != term:
                    continue
                if cdetail[0] == 'low':
                    logger.info('annotation %d is low' % cannotation['annotationid'])
                    continue
                foundit = True
                break
            if not foundit:
                continue
            seqs = self.get_annotation_sequences(cannotation['annotationid'])
            for cseq in seqs:
                feature_num[cseq] += 1
        return feature_num

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
        print('term1 features %d, term2 features %d, common %d, score %f' % (len(term1_f), len(term2_f), common, cscore))
