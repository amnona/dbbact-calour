import requests
import webbrowser
from logging import getLogger

from calour.util import get_config_value
from calour.database import Database

logger = getLogger(__name__)


class DBBact(Database):
    def __init__(self):
        super().__init__(database_name='dbBact', methods=['get', 'annotate', 'feature_terms'])

        # Web address of the bact server
        self.dburl = 'http://amnonim.webfactional.com/scdb_main'
        # self.dburl = 'http://amnonim.webfactional.com/scdb_develop'
        self.username = get_config_value('username', section='dbBact')
        self.password = get_config_value('password', section='dbBact')
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
            for cterm in terms[:min(4, len(terms))]:
                summary += '%s, ' % cterm[0]
            shortdesc.append(({'annotationtype': 'other'}, summary))
            summary = 'special: '
            terms = sorted(terms, key=lambda x: -x[2]/x[1])
            for cterm in terms[:min(4, len(terms))]:
                summary += '%s, ' % cterm[0]
            shortdesc.append(({'annotationtype': 'other'}, summary))
        for cann in annotations:
            cdesc = self.get_annotation_string(cann)
            shortdesc.append((cann, cdesc))
        return shortdesc

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
                        description='', method='', primerid=0, agenttype='Calour', private='n'):
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
        primerid : int (optional)
            the PrimerID from Primers table of the sequences (usually 1 - the V4 515F,806R)
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
        '''
        logger.debug('getting dbbact compact annottions for %d sequences' % len(sequences))
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

        logger.debug('got %d annotations' % (len(annotations)))
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

        address = '%s/annotation_info/%d' % (self.web_interface, annotation['annotationid'])
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

    def get_feature_terms(self, features, exp=None):
        '''Get list of terms per feature

        Parameters
        ----------
        features : list of str
            the features to get the terms for
        exp : calour.Experiment (optional)
            not None to store results inthe exp (to save time for multiple queries)

        Returns
        -------
        feature_terms : dict of list of str/int
            key is the feature, list contains all terms associated with the feature
        '''
        sequence_terms, sequence_annotations, annotations = self.get_seq_list_fast_annotations(features)
        new_annotations = {}
        for cseq, annotations_list in sequence_annotations.items():
            newdesc = []
            for cannotation in annotations_list:
                cdesc = self.get_annotation_string(annotations[cannotation])
                newdesc.append(cdesc)
            new_annotations[cseq] = newdesc
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
