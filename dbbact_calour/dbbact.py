import requests
import webbrowser
from collections import defaultdict
from logging import getLogger

import numpy as np
import pandas as pd

from calour.util import get_config_value
from calour.database import Database
from calour.dsfdr import dsfdr

logger = getLogger(__name__)


class DBBact(Database):
    def __init__(self, exp=None):
        super().__init__(database_name='dbBact', methods=['get', 'annotate', 'enrichment'])

        # Web address of the bact server
        self.dburl = 'http://dbbact.org/REST-API'
        # self.dburl = 'http://amnonim.webfactional.com/scdb_main'
        # self.dburl = 'http://amnonim.webfactional.com/scdb_develop'
        self.username = get_config_value('username', section='dbbact')
        self.password = get_config_value('password', section='dbbact')
        self.web_interface = 'http://dbbact.org'

    def _post(self, api, rdata):
        '''POST a request to dbBact using authorization parameters

        Parameters
        ----------
        api : str
            the REST API address to post the request to
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

    def _get(self, api, rdata):
        '''GET a request to dbBact using authorization parameters

        Parameters
        ----------
        api : str
            the REST API address to post the request to
        rdata : dict
            parameters to pass to the dbBact REST API

        Returns
        -------
        res : request
            the result of the request
        '''
        rdata['user'] = self.username
        rdata['pwd'] = self.password
        res = requests.get(self.dburl + '/' + api, json=rdata)
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

    def _get_all_annotation_string_counts(self, features, exp, ignore_exp=None):
        feature_annotations = {}
        sequence_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
        annotations = exp.exp_metadata['__dbbact_annotations']
        for cseq, annotations_list in sequence_annotations.items():
            if cseq not in features:
                continue
            newdesc = []
            for cannotation in annotations_list:
                if ignore_exp is not None:
                    if annotations[cannotation]['expid'] in ignore_exp:
                        continue
                cdesc = self.get_annotation_string(annotations[cannotation])
                newdesc.append((cdesc, 1))
            feature_annotations[cseq] = newdesc
        return feature_annotations

    def _get_all_term_counts(self, features, feature_annotations, annotations, ignore_exp=None):
        feature_terms = {}
        for cfeature in features:
            annotation_list = []
            for cannotation in feature_annotations[cfeature]:
                if ignore_exp is not None:
                    if annotations[cannotation]['expid'] in ignore_exp:
                        continue
                annotation_list.append(annotations[cannotation])
            # annotation_list = [annotations[x] for x in feature_annotations[cfeature]]
            feature_terms[cfeature] = self.get_annotation_term_counts(annotation_list)
        return feature_terms

    def get_annotation_term_counts(self, annotations):
        '''Get the annotation type corrected count for all terms in annotations

        Parameters
        ----------
        annotations : list of dict
            list of annotations

        Returns
        -------
            list of tuples (term, count)
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
                    logger.warn('unknown detail type %s encountered' % cdesc[0])
                continue
            if cannotation['annotationtype'] == 'other':
                continue
            logger.warn('unknown annotation type %s encountered' % cannotation['annotationtype'])
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

    def add_all_annotations_to_exp(self, exp):
        '''Get annotations for all sequences in experiment and store them in it

        Stores all the annotation details in exp.exp_metadata. Stored key/values are:
        '__dbbact_sequence_terms' : dict of {sequence: list of terms}
            key is sequence, value is list of ontology terms present in the bacteria.
        '__dbbact_sequence_annotations' : dict of {sequence: list of annotationIDs}
            key is sequence, value is list of annotationIDs present in the bacteria.
        '__dbbact_annotations':  dict of {annotationID : annotation_details}
            key is annotaitonID (int), value is the dict of annotation details.

        Parameters
        ----------
        exp : ``Experiment``
            The experiment to get the details for and store them in

        Returns:
        str
        '' if ok, otherwise error string
        '''
        logger.debug('Getting annotations for %d sequences' % len(exp.feature_metadata))
        sequence_terms, sequence_annotations, annotations = self.get_seq_list_fast_annotations(exp.feature_metadata.index.values)
        exp.exp_metadata['__dbbact_sequence_terms'] = sequence_terms
        exp.exp_metadata['__dbbact_sequence_annotations'] = sequence_annotations
        exp.exp_metadata['__dbbact_annotations'] = annotations
        logger.info('Added annotation data to experiment. Total %d annotations, %d terms' % (len(annotations), len(sequence_terms)))
        return ''

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
        '''
        logger.info('Getting dbBact annotations for %d sequences, please wait...' % len(sequences))
        rdata = {}
        rdata['sequences'] = list(sequences)
        res = self._get('sequences/get_fast_annotations', rdata)
        if res.status_code != 200:
            logger.warn('error getting fast annotations for sequence list')
            return None
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

        return sequence_terms, sequence_annotations, res['annotations']

    def show_annotation_info(self, annotation):
        '''Show details about the annotation

        Parameters
        ----------
        annotation : dict
            See dbBact REST API /annotations/get_annotation for keys / values
        '''
        # open in a new tab, if possible
        new = 2

        if 'annotationid' in annotation:
            address = '%s/annotation_info/%d' % (self.web_interface, annotation['annotationid'])
        elif 'sequence' in annotation:
            address = '%s/sequence_annotations/%s' % (self.web_interface, annotation['sequence'])
        else:
            logger.debug('Cannot show annotation info')
            return
        webbrowser.open(address, new=new)

    def add_annotation(self, features, exp):
        '''Open a GUI to Add a new annotation to the selected features

        Parameters
        ----------
        features : list of str
            The sequences to add annotation to the database
        exp : calour.Experiment
            The experiment for which the annotation is added

        Returns
        -------
        '''
        from . import dbannotation

        logger.debug('adding annotation for %d features' % len(features))
        res = dbannotation.annotate_bacteria_gui(self, features, exp)
        return res

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
            sequence_terms, sequence_annotations, annotations = self.get_seq_list_fast_annotations(features)
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
                            cterm = '-'+cdesc[1]
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

    def enrichment(self, exp, features, term_type='term', ignore_exp=None):
        '''Get the list of enriched terms in features compared to all features in exp.

        given uneven distribtion of number of terms per feature

        Parameters
        ----------
        exp : calour.Experiment
            The experiment to compare the features to
        features : list of str
            The features (from exp) to test for enrichmnt
        data_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
             'parentterm' - ontology terms including parent terms associated with each feature.
             'annotation' - the full annotation strings associated with each feature
        ignore_exp: list of int or None (optional)
            List of experiments to ignore in the analysis

        Returns
        -------
        pandas.DataFrame with  info about significantly enriched terms.
            columns:
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
        '''
        exp_features = set(exp.feature_metadata.index.values)
        bg_features = np.array(list(exp_features.difference(features)))

        # add all annotations to experiment if not already added
        if '__dbbact_sequence_terms' not in exp.exp_metadata:
            self.add_all_annotations_to_exp(exp)

        if term_type == 'term':
            feature_terms = self._get_all_term_counts(exp_features, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'], ignore_exp=ignore_exp)
        elif term_type == 'parentterm':
            pass
        elif term_type == 'annotation':
            feature_terms = self._get_all_annotation_string_counts(exp_features, exp=exp, ignore_exp=ignore_exp)
        else:
            raise ValueError('term_type %s not supported for dbbact. possible values are: "term", "parentterm", "annotation"')

        feature_array, term_list = self._get_term_features(features, feature_terms)
        bg_array, term_list = self._get_term_features(bg_features, feature_terms)

        all_feature_array = np.hstack([feature_array, bg_array])

        labels = np.zeros(all_feature_array.shape[1])
        labels[:feature_array.shape[1]] = 1

        keep, odif, pvals = dsfdr(all_feature_array, labels, method='meandiff', transform_type=None, alpha=0.1, numperm=1000, fdr_method='dsfdr')
        keep = np.where(keep)[0]
        if len(keep) == 0:
            logger.info('no enriched terms found')
        term_list = np.array(term_list)[keep]
        odif = odif[keep]
        pvals = pvals[keep]
        res = pd.DataFrame({'term': term_list, 'odif': odif, 'pvals': pvals})
        return res

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
