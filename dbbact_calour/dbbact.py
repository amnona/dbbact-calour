'''
dbbact (:mod:`dbbact_calour.dbbact`)
====================================

.. currentmodule:: dbbact_calour.dbbact

Functions
^^^^^^^^^
.. autosummary::
   :toctree: generated

   DBBact.__init__
   DBBact.version
   DBBact.get_seq_annotation_strings
   DBBact.delete_annotation
   DBBact.remove_feature_from_annotation
   DBBact.get_annotation_website
   DBBact.show_annotation_info
   DBBact.add_all_annotations_to_exp
   DBBact.add_annotation
   DBBact.get_feature_terms
   DBBact.upadte_annotation
   DBBact.enrichment
   DBBact.get_terms_exp
   DBBact.show_term_details
   DBBact.plot_term_annotations
   DBBact.plot_term_annotations_venn
   DBBact.plot_term_venn_all
   DBBact.sample_enrichment
   DBBact.draw_wordcloud
   DBBact.get_enrichment_score
   DBBact.show_enrichment_qt5
'''

from collections import defaultdict
from logging import getLogger
from copy import deepcopy
import webbrowser

import numpy as np
import pandas as pd
import scipy.stats

from .db_access import DBAccess
from .term_pairs import get_enrichment_score, get_terms
from . import __version_numeric__
from calour.util import get_config_value
from calour.database import Database
from calour.experiment import Experiment
from calour.util import _to_list

logger = getLogger(__name__)


class DBBact(Database):
    '''Interface the dbBact database with calour Experiment

    The DBBact class stores the main functions linking a calour.Experiment to dbBact.
    the db attribute stores the DBAccess class with all the calour.Experiment independent dbbact access functions.

    Attributes:
    -----------
    db: DBAccess
        contains the calour independent dbbact access functions
    '''
    def __init__(self, exp=None, dburl='http://api.dbbact.org', web_interface='http://dbbact.org'):
        '''
        Parameters
        ----------
        exp : calour.Experiment, optional
            the experiment to link to
        dburl: str, optional
            the web address for the dbbact REST API
        web_interface: str, optional
            web address for the dbBact website (non REST API) - used for getting/showing dbBact information page
        '''
        super().__init__(database_name='dbBact', methods=['get', 'annotate', 'enrichment'])
        username = get_config_value('username', section='dbbact')
        password = get_config_value('password', section='dbbact')
        dburl = get_config_value('dburl', section='dbbact', fallback=dburl)
        web_interface = get_config_value('web_interface', section='dbbact', fallback=web_interface)
        # dburl = 'http://127.0.0.1:5001'
        # print('Using local database!!!!')
        self.db = DBAccess(dburl=dburl, username=username, password=password)
        self.web_interface = web_interface

    def version(self):
        '''Get the dbbact-calour interface version as a number

        Returns
        -------
        version: float
        '''
        return __version_numeric__

    def get_seq_annotation_strings(self, *kargs, **kwargs):
        '''Get a list of strings describing the sequence annotations, and the annotations details

        Parameters
        ----------
        sequence : str
            the DNA sequence to query the annotation strings about
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        get_summary: bool, optional
            True (default) to get summary of all annotations in the first 3 results of the function (see output details)
            False to just get the annotations

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
        return self.db.get_seq_annotation_strings(*kargs, **kwargs)

    def delete_annotation(self, *kargs, **kwargs):
        return self.db.delete_annotation(*kargs, **kwargs)

    def remove_feature_from_annotation(self, *kargs, **kwargs):
        return self.db.remove_feature_from_annotation(*kargs, **kwargs)

    def get_annotation_website(self, annotation):
        '''Get the database website address of information about the annotation.

        Parameters
        ----------
        annotation : dict
            keys/values are database specific.
            E.g. See dbBact REST API /annotations/get_annotation for keys / values


        Returns
        -------
        str or None
            The webaddress of the html page with details about the annotation,
            or None if not available
        '''
        if 'annotationid' in annotation:
            address = '%s/annotation_info/%d' % (self.web_interface, annotation['annotationid'])
        elif 'sequence' in annotation:
            address = '%s/sequence_annotations/%s' % (self.web_interface, annotation['sequence'])
        else:
            logger.warning('Cannot show annotation info since no annotationid/sequence')
            return None
        return address

    def show_annotation_info(self, annotation):
        '''Show details about the annotation

        Parameters
        ----------
        annotation : dict
            See dbBact REST API /annotations/get_annotation for keys / values
        '''
        # open in a new tab, if possible
        new = 2

        address = self.get_annotation_website(annotation)
        if address is None:
            logger.warning('Cannot find dbBact info address. Aborting')
            return
        webbrowser.open(address, new=new)

    def add_all_annotations_to_exp(self, exp, max_id=None, **kwargs):
        '''Get annotations for all sequences in experiment and store them in it

        Stores all the annotation details in exp.exp_metadata. Stored key/values are:
        '__dbbact_sequence_terms' : dict of {sequence: list of terms}
            key is sequence, value is list of ontology terms present in the bacteria.
        '__dbbact_sequence_annotations' : dict of {sequence: list of annotationIDs}
            key is sequence, value is list of annotationIDs present in the bacteria.
        '__dbbact_annotations':  dict of {annotationID : annotation_details}
            key is annotaitonID (int), value is the dict of annotation details.
        '__dbbact_term_info': dict of {term, {'total_annotations':XXX, 'total_sequences':YYY}}
            number of total annotations and sequences in the database having this term

        Parameters
        ----------
        exp : ``Experiment``
            The experiment to get the details for and store them in
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        **kwargs:
            extra parameters to pass to get_seq_list_fast_annotations()

        Returns
        -------
        str
        '' if ok, otherwise error string
        '''
        logger.debug('Getting annotations for %d sequences' % len(exp.feature_metadata))
        sequence_terms, sequence_annotations, annotations, term_info, taxonomy = self.db.get_seq_list_fast_annotations(exp.feature_metadata.index.values, max_id=max_id, **kwargs)
        exp.exp_metadata['__dbbact_sequence_terms'] = sequence_terms
        exp.exp_metadata['__dbbact_sequence_annotations'] = sequence_annotations
        exp.exp_metadata['__dbbact_annotations'] = annotations
        exp.exp_metadata['__dbbact_term_info'] = term_info
        exp.exp_metadata['__dbbact_taxonomy'] = taxonomy
        logger.info('Added annotation data to experiment. Total %d annotations, %d terms' % (len(annotations), len(sequence_terms)))
        return ''

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

    def get_feature_terms(self, features, exp=None, term_type=None, ignore_exp=None, term_method=('single'), max_id=None, **kwargs):
        '''Get dict of terms scores per feature

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
            'fscore': ...
        ignore_exp : list of int or None (optional)
            the list of experimentids to ignore (don't take info from them)
        term_method: list of str, optional
            the methods to get all the terms for each feature. can include:
                'singe': get the single terms per each feature (i.e. 'feces', '-control', etc.)
                'pairs': get the term pairs for each feature (i.e. 'feces+homo sapiens', etc.)
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        kwargs:
            Parameters to pass to db_access.get_seq_list_fast_annotations. can include:
            get_taxonomy=False, get_parents=False, get_term_info=True

        Returns
        -------
        feature_terms : dict of {feature(str: dict of {term(str): score(float)})}
            key is the feature, value is a dict of the score (value) for each term(key) for this feature
        '''
        if term_type is None:
            term_type = 'terms'
        if exp is not None:
            if '__dbbact_sequence_terms' not in exp.exp_metadata:
                # if annotations not yet in experiment - add them
                self.add_all_annotations_to_exp(exp, max_id=max_id, **kwargs)
            # and filter only the ones relevant for features
            sequence_terms = exp.exp_metadata['__dbbact_sequence_terms']
            sequence_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
            annotations = exp.exp_metadata['__dbbact_annotations']
            term_info = exp.exp_metadata['__dbbact_term_info']
            taxonomy = exp.exp_metadata['__dbbact_taxonomy']
        else:
            sequence_terms, sequence_annotations, annotations, term_info, taxonomy = self.db.get_seq_list_fast_annotations(features, max_id=max_id, **kwargs)

        # get the current experimentID to ignore if ignore_exp is True
        if ignore_exp is True:
            ignore_exp = self.db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['sample_metadata_md5'], getall=True)
            if ignore_exp is None:
                logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
            else:
                logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))

        new_annotations = {}
        term_scores = {}
        if term_type == 'annotation':
            for cseq, annotations_list in sequence_annotations.items():
                if cseq not in features:
                    continue
                term_scores[cseq] = defaultdict(float)
                newdesc = []
                for cannotation in annotations_list:
                    if ignore_exp is not None:
                        annotationexp = annotations[cannotation]['expid']
                        if annotationexp in ignore_exp:
                            continue
                    cdesc = self.db.get_annotation_string(annotations[cannotation])
                    term_scores[cseq][cdesc] += 1
                    newdesc.append(cdesc)
                new_annotations[cseq] = newdesc
        elif term_type == 'terms':
            for cseq, annotations_list in sequence_annotations.items():
                if cseq not in features:
                    continue
                term_scores[cseq] = defaultdict(float)
                newdesc = []
                annotation_terms = []
                for cannotation in annotations_list:
                    if ignore_exp is not None:
                        annotationexp = annotations[cannotation]['expid']
                        if annotationexp in ignore_exp:
                            continue
                    annotation_terms = get_terms(annotations[cannotation], term_types=term_method)
                    for cterm in annotation_terms:
                        term_scores[cseq][cterm] += 1
                    # for ctype, cterm in annotations[cannotation]['details']:
                    #     if ctype == 'low':
                    #         cterm = '-' + cterm
                    #     newdesc.append(cterm)
                    #     term_scores[cseq][cterm] += 1
                new_annotations[cseq] = annotation_terms

        # f-score for each term
        elif term_type == 'fscore':
            if ignore_exp is None:
                ignore_exp = []
            # we need to rekey the annotations with an str (old problem...)
            annotations = {str(k): v for k, v in annotations.items()}
            for cseq, annotations_list in sequence_annotations.items():
                if cseq not in features:
                    continue
                term_scores[cseq] = defaultdict(float)
                fscore, recall, precision, term_count, reduced_f = get_enrichment_score(annotations, [(cseq, annotations_list)], ignore_exp=ignore_exp, term_info=term_info, term_types=term_method)
                if len(fscore) == 0:
                    new_annotations[cseq] = ['NA']
                    continue
                sorted_fscore = sorted(fscore.items(), key=lambda x: x[1], reverse=True)
                new_annotations[cseq] = [sorted_fscore[0][0]]
                term_scores[cseq] = fscore
        elif term_type == 'parentterms':
            for cseq, term_list in sequence_terms.items():
                if cseq not in features:
                    continue
                term_scores = defaultdict(float)
                term_list = [x for x in term_list if x != 'na']
                for cterm in term_list:
                    term_scores[cseq][cterm] += 1
                new_annotations[cseq] = term_list
        elif term_type == 'contamination':
            for cseq, annotations_list in sequence_annotations.items():
                term_scores[cseq] = defaultdict(float)
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
                    term_scores[cseq]['contamination'] = is_contamination
                else:
                    new_annotations[cseq] = []
        else:
            raise ValueError('term_type %s not supported in get_feature_terms. Possible values are "annotation","terms","parentterms","fscore","contamination"' % term_type)
        return term_scores
        # return new_annotations
        # return sequence_annotations
        # return sequence_terms

    def filter_features_based_on_terms(self, exp, terms, filter_method='any', term_types=('single'), ignore_exp=None, negate=False, max_id=None):
        '''filter features based on how many experiments they appear in

        Parameters
        ----------
        exp: AmpliconExperiment
            The experiment to filter features from
        terms: dict of {term(str): min_experiments(int)}
            the minimal number of experiments (value) for each term (key) used in the filtering.
            keep the feature if we observe the tern in >= min_experiments
        filter_method: str, optional
            options are:
            'any': keep if any of the terms satisfy the minimal number of experiments (i.e. one term ok will keep the feature)
            'all': keep only if all terms satisfy the minimal number of experimetns (i.e. one term not ok will remove the features)
        term_types: list of str, optional
            'single': use single terms from the annotations (i.e. 'feces', '-control')
            'pairs': use term pairs from the annotations (i.e. 'feces+homo sapiens')
        ignore exp: list of int or None, optional
            list of xperiment IDs to ignore annotations from
        negate: bool, optional
            True to reverse the filtering (keep instead of remove)
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        '''
        if '__dbbact_sequence_terms' not in exp.exp_metadata:
                # if annotations not yet in experiment - add them
            self.add_all_annotations_to_exp(exp, max_id=max_id)
        # get the dbbact annotations
        sequence_terms = exp.exp_metadata['__dbbact_sequence_terms']
        sequence_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
        annotations = exp.exp_metadata['__dbbact_annotations']
        term_info = exp.exp_metadata['__dbbact_term_info']

        # filter
        keep_features = []
        for cseq, annotations_list in sequence_annotations.items():
            if cseq not in exp.feature_metadata.index:
                continue
            term_scores = defaultdict(float)
            for cannotation in annotations_list:
                if ignore_exp is not None:
                    annotationexp = annotations[cannotation]['expid']
                    if annotationexp in ignore_exp:
                        continue
                annotation_terms = get_terms(annotations[cannotation], term_types=term_types)
                for cterm in annotation_terms:
                    term_scores[cterm] += 1
            if filter_method == 'any':
                keep = False
                for cterm, cterm_min_count in terms.items():
                    if cterm not in term_scores:
                        continue
                    if term_scores[cterm] >= cterm_min_count:
                        keep = True
                        break
            elif filter_method == 'all':
                keep = True
                for cterm, cterm_min_count in terms.items():
                    if cterm not in term_scores:
                        keep = False
                        break
                    if term_scores[cterm] < cterm_min_count:
                        keep = False
                        break
            if keep:
                keep_features.append(cseq)
        newexp = exp.filter_ids(keep_features, axis='f', negate=negate)
        return newexp

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
        dbannotation.update_annotation_gui(self.db, annotation, exp)

    def enrichment(self, exp, features, max_id=None, **kwargs):
        '''Get the list of enriched terms in features compared to all other features in exp.

        given uneven distribtion of number of terms per feature

        Parameters
        ----------
        exp : calour.Experiment
            The experiment to compare the features to.
            NOTE: exp must contain the
        features : list of str
            The features (from exp) to test for enrichmnt (comapred to the other features in exp)
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id

        **kwargs: additional parameters supplied to db_access.term_enrichment(). These include:
        term_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
            'parentterm' - ontology terms including parent terms associated with each feature.
            'annotation' - the full annotation strings associated with each feature
            'combined' - combine 'term' and 'annotation'
        ignore_exp: list of int or None or True, optional
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment
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
        random_seed: int or None
            int to specify the random seed for numpy.random.
        use_term_pairs: bool, optional
            True to also test enrichment in pairs of terms (i.e. homo sapiens+feces, etc.)
        focus_terms: list of str or None, optional
            if not None, use only annotations containing all the terms in focus_terms

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
            num_enriched_exps: number of experiments where the term is significantly enriched
            num_enriched_annotations: number of annotations where the term is significantly enriched
            num_total_exps: number of experiments with this annotation in the term annotations list
            description : the term (str)
        numpy.Array where rows are features (ordered like the dataframe), columns are features and value is score
            for term in feature
        pandas.DataFrame with info about the features used. columns:
            group: int the group (1/2) to which the feature belongs
            sequence: str
        '''
        exp_features = set(exp.feature_metadata.index.values)
        bad_features = set(features).difference(exp_features)
        if len(bad_features) > 0:
            logger.warning('Some of the features for enrichment are not in the experiment. Ignoring %d features' % len(bad_features))
            features = list(set(features).intersection(exp_features))
            if len(features) == 0:
                raise ValueError("No features left after ignoring. Please make sure you test for enrichment with features from the experiment.")
        bg_features = np.array(list(exp_features.difference(features)))

        # add all annotations to experiment if not already added
        if '__dbbact_sequence_terms' not in exp.exp_metadata:
            self.add_all_annotations_to_exp(exp, max_id=max_id)

        ignore_exp = kwargs.get('ignore_exp')
        # if ignore exp is True, it means we should ignore the current experiment
        if ignore_exp is True:
            ignore_exp = self.db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['sample_metadata_md5'], getall=True)
            if ignore_exp is None:
                logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
            else:
                logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))
        kwargs['ignore_exp'] = ignore_exp

        res = self.db.term_enrichment(g1_features=features, g2_features=bg_features, all_annotations=exp.exp_metadata['__dbbact_annotations'], seq_annotations=exp.exp_metadata['__dbbact_sequence_annotations'], term_info=exp.exp_metadata.get('__dbbact_term_info'), **kwargs)
        return res

    def enrichmentcount(self, exp, features, **kwargs):
        '''DEPRACATED!!!!!!!
        Get the list of enriched terms in features compared to all features in exp.

        given uneven distribtion of number of terms per feature

        Parameters
        ----------
        exp : calour.Experiment
            The experiment to compare the features to
        features : list of str
            The features (from exp) to test for enrichmnt
        term_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
            'parentterm' - ontology terms including parent terms associated with each feature.
            'annotation' - the full annotation strings associated with each feature
            'combined' - combine 'term' and 'annotation'
        ignore_exp: list of int or None or True, optional
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment
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
        random_seed: int or None
            int to specify the random seed for numpy.random.
        use_term_pairs: bool, optional
            True to also test enrichment in pairs of terms (i.e. homo sapiens+feces, etc.)

        Returns
        -------
        count
        '''
        exp_features = set(exp.feature_metadata.index.values)
        bad_features = set(features).difference(exp_features)
        if len(bad_features) > 0:
            logger.warning('Some of the features for enrichment are not in the experiment. Ignoring %d features' % len(bad_features))
            features = list(set(features).intersection(exp_features))
            if len(features) == 0:
                raise ValueError("No features left after ignoring. Please make sure you test for enrichment with features from the experiment.")
        bg_features = np.array(list(exp_features.difference(features)))

        # add all annotations to experiment if not already added
        if '__dbbact_sequence_terms' not in exp.exp_metadata:
            self.add_all_annotations_to_exp(exp)

        ignore_exp = kwargs.get('ignore_exp')
        # if ignore exp is True, it means we should ignore the current experiment
        if ignore_exp is True:
            ignore_exp = self.db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['sample_metadata_md5'], getall=True)
            if ignore_exp is None:
                logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
            else:
                logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))
        kwargs['ignore_exp'] = ignore_exp

        res = self.db.get_term_enriched_annotations(g1_features=features, g2_features=bg_features, all_annotations=exp.exp_metadata['__dbbact_annotations'], seq_annotations=exp.exp_metadata['__dbbact_sequence_annotations'], **kwargs)
        return res

    def get_terms_exp(self, exp, term):
        '''Get an experiment with features (from exp) as columns, annotations cotaining terms as rows

        Parameters
        ----------
        exp: Experiment
        term: str
            The term to get annotations that include one of these terms

        Returns
        -------
        Experiment with features (from exp) as columns, annotations cotaining terms as rows
        '''
        tmat, tanno, tseqs = self.db.get_term_annotations(term, list(exp.feature_metadata.index.values), exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        newexp = Experiment(tmat, sample_metadata=tseqs, feature_metadata=tanno)
        newexp = newexp.cluster_features(1)
        newexp = newexp.sort_by_metadata(field='expid', axis='f')
        newexp.plot(feature_field='annotation', gui='qt5', yticklabel_kwargs={'rotation': 0}, yticklabel_len=35, cmap='tab20b', norm=None, bary_fields=['expid'], bary_label=False)

    def show_term_details(self, term, exp, features, group2_features, group1_name='group1', group2_name='group2', max_id=None, **kwargs):
        '''
        Plot a heatmap for all annotations containing term in experiment
        Rows are the annotations, columns are the sequences (sorted by features/group2_features)

        Parameters
        ----------
        term: str
            the term to plot for
        exp: calour.Experiment
            the experiment where the analysis should look at
        features, group2_features: list of str
            the list of sequences for the 2 groups to compare
        group1_name, group2_name: str, optional
            name for group1/ group2 (for the plot)
        **kwargs:
            passed to calour.plot() (i.e. gui='qt5', etc.)
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id

        Returns
        -------
        calour.Experiment
            with rows as annotations, columns as features from the 2 groups
        '''
        if term[0] == '-':
            term = term[1:]
        term = term.lstrip('LOWER IN ')
        term = term.rstrip(' *')
        all_seqs = list(features.copy())
        all_seqs.extend(list(group2_features))
        if '__dbbact_sequence_annotations' not in exp.exp_metadata:
            self.add_all_annotations_to_exp(exp, max_id=max_id)
        tmat, tanno, tseqs = self.db.get_term_annotations(term, all_seqs, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        seq_group = [str(group1_name)] * len(features)
        seq_group.extend([str(group2_name)] * len(group2_features))
        tseqs['group'] = seq_group
        newexp = Experiment(tmat.T, sample_metadata=tanno, feature_metadata=tseqs)
        # newexp = newexp.cluster_features(1)
        newexp = newexp.cluster_data(axis='s')
        # newexp = newexp.sort_by_metadata(field='expid', axis='s')
        newexp.plot(cmap='tab20b', norm=None, barx_fields=['expid'], barx_label=False, bary_fields=['group'], bary_label=True, **kwargs)
        return newexp

    def show_term_details_diff(self, term, exp, **kwargs):
        '''Plot all the annotations in a diff_abundance result exp for a given term, dividung according to the two groups

        Parameters
        ----------
        term: str
            the dbbact term to examine
        exp: calour.Experiment
            results of diff_abundance()
        **kwargs: passed to exp.plot() (i.e. gui='qt5', etc)

        Returns
        -------
        ca.Experiment
            with annotations as rows, features as columns
        '''
        names1 = exp.feature_metadata['_calour_direction'][exp.feature_metadata['_calour_stat'] < 0]
        if len(names1) > 0:
            group1_name = names1.values[0]
        else:
            group1_name = 'group1'
        names2 = exp.feature_metadata['_calour_direction'][exp.feature_metadata['_calour_stat'] > 0]
        if len(names2) > 0:
            group2_name = names2.values[0]
        else:
            group2_name = 'group2'

        negative = exp.feature_metadata['_calour_stat'] < 0
        group1_features = exp.feature_metadata.index.values[negative.values]
        positive = exp.feature_metadata['_calour_stat'] > 0
        group2_features = exp.feature_metadata.index.values[positive.values]

        newexp = self.show_term_details(term, exp, group1_features, group2_features, group1_name=group1_name, group2_name=group2_name, **kwargs)
        return newexp

    def plot_term_annotations(self, term, exp, features, group2_features, min_prevalence=0.01):
        '''Plot a nice graph summarizing all the annotations supporting the term in the 2 groups
        '''
        # from matplotlib import rc
        # rc('text', usetex=False)
        if term[0] == '-':
            term = term[1:]
        all_seqs = features.copy()
        all_seqs.extend(group2_features)
        tmat, tanno, tseqs = self.db.get_term_annotations(term, all_seqs, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        seq_group = np.ones(len(all_seqs))
        seq_group[:len(features)] = 0
        tseqs['group'] = seq_group
        newexp = Experiment(tmat, sample_metadata=tseqs, feature_metadata=tanno)
        newexp = newexp.filter_prevalence(min_prevalence)
        newexp = newexp.sort_by_metadata('expid', axis='f')
        # experiments = newexp.feature_metadata['expid'].unique()
        experiments = newexp.feature_metadata['expid']
        g1len = len(features)
        g2len = len(group2_features)

        import matplotlib.pyplot as plt
        colors = ['r', 'g', 'b', 'c', 'm', 'k']
        cdict = {}
        for idx, cexpid in enumerate(newexp.feature_metadata['expid'].unique()):
            cdict[cexpid] = colors[np.mod(idx, len(colors))]
        nrows = int(np.ceil(np.sqrt(len(experiments))))
        ncols = int(np.ceil(len(experiments) / nrows))
        plt.subplots(nrows=nrows, ncols=ncols, figsize=[15, 15])
        cexp = newexp
        # fig = plt.figure()
        # for idx, cexpid in enumerate(experiments):
        #     cexp = newexp.filter_by_metadata('expid',[cexpid],axis='f')
        #     plt.subplot(nrows, ncols, idx+1)
        #     plt.title(cexpid)
        for idx2, canno in enumerate(cexp.feature_metadata.iterrows()):
            canno = canno[1]
            plt.subplot(nrows, ncols, idx2 + 1)
            numg1 = (cexp.data[:g1len, idx2] > 0).sum()
            numg2 = (cexp.data[g1len:, idx2] > 0).sum()
            if canno['detail_type'] == 'low':
                mult = -1
            else:
                mult = 1
            # plt.pie([numg1, numg2, numother], colors=['r','g','b'])
            plt.pie([numg1 / g1len, 1 - numg1 / g1len], colors=['mediumblue', 'aliceblue'], center=[0.6, 0], radius=0.5)
            plt.pie([numg2 / g2len, 1 - numg2 / g2len], colors=['r', 'mistyrose'], center=[-0.6, 0], radius=0.5)
            # plt.bar(0,mult*numg1/g1len, width=0.1, color='r')
            # plt.bar(0.15,mult*numg2/g2len, width=0.1, color='b')
            # plt.barh(np.arange(2),[numg1/g1len, numg2/g2len])
            # plt.legend(['group1','group2','none'])
            ctitle = canno['annotation']
            tt = ctitle.split(' ')
            clen = 0
            otitle = ''
            for ttt in tt:
                otitle += ttt
                clen += len(ttt)
                if ttt == term:
                    ttt = r"$\bf{" + ttt + "}$"
                if clen > 20:
                    otitle += '\n'
                    clen = 0
                else:
                    otitle += ' '
                if len(otitle) > 100:
                    break
            plt.text(0, 0.5, otitle, fontsize=10, color=cdict[canno['expid']], horizontalalignment='center', verticalalignment='bottom')
            diff_title_high = []
            diff_title_low = []
            all_title = []
            print(canno)
            cannotation = exp.exp_metadata['__dbbact_annotations'][canno['annotationid']]
            for cdetailtype, cdetailterm in cannotation['details']:
                if cdetailtype == 'all':
                    all_title.append(cdetailterm)
                elif cdetailtype == 'high':
                    diff_title_high.append(cdetailterm)
                elif cdetailtype == 'low':
                    diff_title_low.append(cdetailterm)
            diff_title = ', '.join(diff_title_high) + ' > ' + ', '.join(diff_title_low)
            all_title = ', '.join(all_title)
            plt.text(0, -0.7, diff_title, fontsize=20, color='black', horizontalalignment='center', verticalalignment='bottom')
            plt.text(0, 0.7, all_title, fontsize=20, color='black', horizontalalignment='center', verticalalignment='bottom')
            # plt.title(otitle,fontsize=10, color=cdict[canno['expid']])
            # plt.ylim([-1,1])
            # plt.plot([-0.1,0.25],[0,0],':k')
            # plt.xlim([-0.1,0.25])
            # plt.xticks([], [])
        return plt.gcf()

    def plot_term_annotations2(self, term, exp, features, group2_features, min_prevalence=0.01):
        '''Plot a nice graph summarizing all the annotations supporting the term in the 2 groups
        '''
        # from matplotlib import rc
        # rc('text', usetex=False)
        if term[0] == '-':
            term = term[1:]
        all_seqs = features.copy()
        all_seqs.extend(group2_features)
        tmat, tanno, tseqs = self.db.get_term_annotations(term, all_seqs, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        seq_group = np.ones(len(all_seqs))
        seq_group[:len(features)] = 0
        tseqs['group'] = seq_group
        newexp = Experiment(tmat, sample_metadata=tseqs, feature_metadata=tanno)
        newexp = newexp.filter_prevalence(min_prevalence)
        newexp = newexp.sort_by_metadata('expid', axis='f')
        # experiments = newexp.feature_metadata['expid'].unique()
        experiments = newexp.feature_metadata['expid']
        g1len = len(features)
        g2len = len(group2_features)

        import matplotlib.pyplot as plt
        colors = ['r', 'g', 'b', 'c', 'm', 'k']
        cdict = {}
        for idx, cexpid in enumerate(newexp.feature_metadata['expid'].unique()):
            cdict[cexpid] = colors[np.mod(idx, len(colors))]
        nrows = int(np.ceil(np.sqrt(len(experiments))))
        ncols = int(np.ceil(len(experiments) / nrows))
        plt.subplots(nrows=nrows, ncols=ncols, figsize=[15, 15])
        cexp = newexp
        # fig = plt.figure()
        # for idx, cexpid in enumerate(experiments):
        #     cexp = newexp.filter_by_metadata('expid',[cexpid],axis='f')
        #     plt.subplot(nrows, ncols, idx+1)
        #     plt.title(cexpid)
        for idx2, canno in enumerate(cexp.feature_metadata.iterrows()):
            canno = canno[1]
            plt.subplot(nrows, ncols, idx2 + 1)
            numg1 = (cexp.data[:g1len, idx2] > 0).sum()
            numg2 = (cexp.data[g1len:, idx2] > 0).sum()
            if canno['detail_type'] == 'low':
                mult = -1
            else:
                mult = 1
            # plt.pie([numg1, numg2, numother], colors=['r','g','b'])
            plt.bar(0, mult * numg1 / g1len, width=0.1, color='r')
            plt.bar(0.15, mult * numg2 / g2len, width=0.1, color='b')
            # plt.barh(np.arange(2),[numg1/g1len, numg2/g2len])
            # plt.legend(['group1','group2','none'])
            ctitle = canno['annotation']
            tt = ctitle.split(' ')
            clen = 0
            otitle = ''
            for ttt in tt:
                otitle += ttt
                clen += len(ttt)
                if ttt == term:
                    ttt = r"$\bf{" + ttt + "}$"
                if clen > 20:
                    otitle += '\n'
                    clen = 0
                else:
                    otitle += ' '
                if len(otitle) > 100:
                    break
            plt.title(otitle, fontsize=10, color=cdict[canno['expid']])
            plt.ylim([-1, 1])
            plt.plot([-0.1, 0.25], [0, 0], ':k')
            plt.xlim([-0.1, 0.25])
            plt.xticks([], [])
        return plt.gcf()

    def plot_term_venn_all(self, terms, exp, bacteria_groups=None, set_colors=('red', 'green', 'mediumblue'), max_size=None, ignore_exp=[], max_id=None, use_exact=True):
        '''Plot a venn diagram for all sequences appearing in any annotation containing the term, and intersect with both groups of bacteria

        Parameters
        ----------
        terms: str or list of str
            the terms to test the overlap for. if more than one term supplied (list of str), look for all otus in the overlap
        exp: calour.Experiment or None
            the experiment containing the bacteria groups or None to use the bacteria_groups parameter
        bacteria_groups: (list of str, list of str) or None, optional
            if not None, use bacteria sequences from the two lists as the two groups.
            If None, use groups from exp
        set_colors: tuple of (str, str, str), optional
            The colors for group1, group2, term circles
        max_size: int or None, optional
            if not None, clip term circle size to max_size.
            Used to make cases where term has lots of sequences nicer.
            NOTE: it changes the circle size and number!
        ignore_exp: list of int or None or True, optional
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment
            None (default) to use annotations from all experiments including the current one
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        use_exact: bool, optional
            True (default) to search only for annotations exactly matching the query sequence (region and length)
            False to search including sequences from other regions (using SILVA based sequence translator) and also shorter/longer sequences.
            NOTE: if use_exact=False, the venn diagramm will include more sequences than the query sequence
        '''
        import matplotlib.pyplot as plt
        try:
            from matplotlib_venn import venn3
        except Exception as err:
            print(err)
            raise ValueError("Error importing matplotlib_venn. Is it installed? If not, install it using: pip install matplotlib-venn")

        # get the two bacteria groups
        if bacteria_groups is None:
            vals = exp.feature_metadata['_calour_direction'].unique()
            group1 = list(exp.feature_metadata.index[exp.feature_metadata['_calour_direction'] == vals[0]].values)
            group1_name = vals[0]
            group2 = list(exp.feature_metadata.index[exp.feature_metadata['_calour_direction'] == vals[1]].values)
            group2_name = vals[1]
        else:
            group1 = bacteria_groups[0]
            group2 = bacteria_groups[1]
            group1_name = 'group1'
            group2_name = 'group2'

        logger.debug('group1 (%s) :%d sequences, group2 (%s): %d sequences' % (group1_name, len(group1), group2_name, len(group2)))

        # get the dbbact sequenceIDs for the bacteria
        if use_exact:
            g1idst = self.db.get_sequences_ids(list(group1), no_shorter=True, no_longer=True, use_sequence_translator=False)
            g2idst = self.db.get_sequences_ids(list(group2), no_shorter=True, no_longer=True, use_sequence_translator=False)
        else:
            g1idst = self.db.get_sequences_ids(list(group1), no_shorter=False, no_longer=False, use_sequence_translator=True)
            g2idst = self.db.get_sequences_ids(list(group2), no_shorter=False, no_longer=False, use_sequence_translator=True)

        # convert from list of list to a set
        g1ids = set()
        for idx, cs in enumerate(g1idst):
            # if no sequence ids found for current sequence, keep it so venn diagram size will be correct
            # (it belongs to the non-intersect)
            if len(cs) == 0:
                g1ids.add('g1_%d' % idx)
                continue
            for cid in cs:
                g1ids.add(cid)

        g2ids = set()
        for idx, cs in enumerate(g2idst):
            # if no sequence ids found for current sequence, keep it so venn diagram size will be correct
            # (it belongs to the non-intersect)
            if len(cs) == 0:
                g2ids.add('g2_%d' % idx)
                continue
            for cid in cs:
                g2ids.add(cid)

        if not use_exact:
            logger.debug('after inexact match addition, group1 (%s): %d sequences, group2 (%s): %d sequences' % (group1_name, len(g1ids), group2_name, len(g2ids)))

        # if ignore exp is True, it means we should ignore the current experiment
        if ignore_exp is True:
            if exp is None:
                raise ValueError('Cannot ignore current experiment when exp=None. Please explicitly specify ignore_exp=[ID1,ID2...]')
            ignore_exp = self.db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['sample_metadata_md5'], getall=True)
            if ignore_exp is None:
                logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
                ignore_exp = []
            else:
                logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))

        # get the term sequence ids
        terms = _to_list(terms)
        termids = None

        # TODO: fix the lower in test
        # remove the "-" for the terms
        new_terms = []
        for cterm in terms:
            if cterm[0] == '-':
                cterm = cterm[1:]
            new_terms.append(cterm)
        terms = new_terms

        # get the sequence ids that have these terms
        termids = set(self.db.get_db_term_features(terms, ignore_exp=ignore_exp, max_id=max_id))

        og1 = len(termids.intersection(g1ids))
        og2 = len(termids.intersection(g2ids))
        ogg = 0
        oga = 0
        vvals = {}
        vvals['100'] = len(g1ids) - og1 - ogg
        vvals['010'] = len(g2ids) - og2 - ogg
        vvals['001'] = len(termids) - og1 - og2
        vvals['110'] = ogg
        vvals['101'] = og1
        vvals['011'] = og2
        vvals['111'] = oga

        if max_size is not None:
            if vvals['001'] > max_size:
                print('WARNING: clipped term circle size to %d. Real size (number of term seqs not overlapping) should be: %d' % (max_size, len(termids) - og1 - og2))
                vvals['001'] = max_size

        # for k, v in vvals.items():
        #     if v > 0:
        #         vvals[k] = np.log2(v)

        f = plt.figure()
        termstr = '+'.join(terms)
        # venn3([gg1, gg2, anno_group], set_labels=(group1_name, group2_name, 'annotation'), set_colors=['mediumblue', 'r', 'g'])
        venn3(vvals, set_labels=(group1_name, group2_name, termstr), set_colors=set_colors)
        plt.title('term overlaps for %s' % termstr)
        return f

    def plot_term_annotations_venn(self, term, exp, bacteria_groups=None, min_prevalence=0, annotation_types=None, set_colors=('red', 'green', 'mediumblue'), min_overlap=0, min_size=0):
        '''Plot a venn diagram for all annotations containing the term, showing overlap between the term and the two bacteria groups
        Parameters
        ----------
        term: str
            The term to get the annotations for
        exp: calour.Experiment
            The experiment containing the dbbact annotations
        bacteria_group: tuple of list of str or None, optional
            the list of bacterial sequences for the two groups to compare.
            None to get the bacteria from the diff. abundance results (assuming exp is after diff_abundance())
        min_prevalens: float, optional
            Only include annotations that have at least min_prevalence overlap with either of the groups
        annotation_types: list of str or None, optional
            Only include annotations of types in the list. can contain:
                'common', 'high_freq', 'diff', 'other'
            None to include all
            Not implemented yet
        set_colors: tuple of str, optional
            Colors to use for group1. group2, annotation
        min_overlap: float, optional
            The minimal fraction of annotation bacteria that overlap with either of both groups in order to show annotation
        min_size: int, optional
            minimal number of bacteria in the annotation in order to show it

        Returns
        -------
        list of figures
        '''
        try:
            from matplotlib_venn import venn3
        except Exception as err:
            print(err)
            raise ValueError("Error importing matplotlib_venn. Is it installed? If not, install it using: pip install matplotlib-venn")

        # get the two bacteria groups
        if bacteria_groups is None:
            vals = exp.feature_metadata['_calour_direction'].unique()
            group1 = list(exp.feature_metadata.index[exp.feature_metadata['_calour_direction'] == vals[0]].values)
            group1_name = vals[0]
            group2 = list(exp.feature_metadata.index[exp.feature_metadata['_calour_direction'] == vals[1]].values)
            group2_name = vals[1]
        else:
            group1 = bacteria_groups[0]
            group2 = bacteria_groups[1]
            group1_name = 'group1'
            group2_name = 'group2'

        if term[0] == '-':
            term = term[1:]
        all_seqs = group1.copy()
        all_seqs.extend(group2)
        tmat, tanno, tseqs = self.db.get_term_annotations(term, all_seqs, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        seq_group = np.ones(len(all_seqs))
        seq_group[:len(group1)] = 0
        tseqs['group'] = seq_group
        newexp = Experiment(tmat, sample_metadata=tseqs, feature_metadata=tanno)
        newexp = newexp.filter_prevalence(min_prevalence)
        newexp = newexp.sort_by_metadata('expid', axis='f')

        import matplotlib.pyplot as plt
        all_figures = []
        for idx2, canno in enumerate(newexp.feature_metadata.iterrows()):
            canno = canno[1]
            aseqs = self.db.get_annotation_sequences(int(canno['annotationid']))
            gg1 = set([exp.feature_metadata.index.get_loc(x) for x in group1])
            gg2 = set([exp.feature_metadata.index.get_loc(x) for x in group2])
            cdata = newexp.get_data(sparse=False)
            anno_group = set(np.where(cdata[:, idx2] > 0)[0])
            og1 = len(gg1.intersection(anno_group))
            og2 = len(gg2.intersection(anno_group))
            ogg = 0
            oga = 0
            vvals = {}
            vvals['100'] = len(gg1) - og1 - ogg
            vvals['010'] = len(gg2) - og2 - ogg
            vvals['001'] = len(aseqs) - og1 - og2
            vvals['110'] = ogg
            vvals['101'] = og1
            vvals['011'] = og2
            vvals['111'] = oga

            annotation_overlap = (og1 + og2) / len(aseqs)
            if annotation_overlap < min_overlap:
                continue
            if len(aseqs) < min_size:
                continue

            f = plt.figure()
            # venn3([gg1, gg2, anno_group], set_labels=(group1_name, group2_name, 'annotation'), set_colors=['mediumblue', 'r', 'g'])
            venn3(vvals, set_labels=(group1_name, group2_name, 'annotation'), set_colors=set_colors)
            plt.title('%s (expid %s)' % (canno['annotation'], canno['expid']))
            all_figures.append(f)
        return all_figures

    def sample_enrichment(self, exp, field, value1, value2=None, term_type='term', ignore_exp=None, min_appearances=3, fdr_method='dsfdr', score_method='all_mean', freq_weight='log', alpha=0.1, use_term_pairs=False, max_id=None):
        '''Get the list of enriched terms for all bacteria between two groups using frequencies from the Experiment.data table.

        It is equivalent to multiplying the (freq_weight transformed) feature X sample matrix by the database derived term X feature matrix
        (where the numbers are how strong is the term associated with the feature based on database annotations using score_method).
        A differntial abundance test (using dsFDR multiple hypothesis correction) is then performed on the resulting sample X term matrix.

        Parameters
        ----------
        exp : calour.Experiment
            The experiment to compare the features to
        field : str
            Name of the field to divide the samples by
        value1 : str or list of str
            Values (for selected field) for samples belonging to sample group 1
        value2 : str or list of str or None (optional)
            None (default) to select all samples not belonging to sample group 1
            str or list of str - Values (for selected field) for samples belonging to sample group 2.
        term_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
            'parentterm' - ontology terms including parent terms associated with each feature.
            'annotation' - the full annotation strings associated with each feature
            'combined' - combine 'term' and 'annotation'
        ignore_exp: list of int or None or True, optional
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment
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
                'rank' : use the rank of each feature across samples
        alpha : float (optional)
            the FDR level desired (0.1 means up to 10% of results can be due to null hypothesis)
        Use_term_pairs: bool, optional
            True to get all term pair annotations (i.e. human+feces etc.)
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id

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
        exp_features = list(exp.feature_metadata.index.values)

        # add all annotations to experiment if not already added
        if '__dbbact_sequence_terms' not in exp.exp_metadata:
            self.add_all_annotations_to_exp(exp, max_id=max_id)

        # if ignore exp is True, it means we should ignore the current experiment
        if ignore_exp is True:
            ignore_exp = self.db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['sample_metadata_md5'], getall=True)
            if ignore_exp is None:
                logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
            else:
                logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))

        all_annotations = exp.exp_metadata['__dbbact_annotations']
        seq_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
        if term_type == 'term':
            feature_terms = self.db._get_all_term_counts(exp_features, seq_annotations, all_annotations, ignore_exp=ignore_exp, score_method=score_method, use_term_pairs=use_term_pairs)
        elif term_type == 'parentterm':
            pass
        elif term_type == 'annotation':
            feature_terms = self.db._get_all_annotation_string_counts(exp_features, all_anntations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
        elif term_type == 'combined':
            feature_terms = self.db._get_all_term_counts(exp_features, seq_annotations, all_annotations, ignore_exp=ignore_exp, score_method=score_method, use_term_pairs=use_term_pairs)
            feature_annotations = self.db._get_all_annotation_string_counts(exp_features, all_annotations=all_annotations, seq_annotations=seq_annotations, ignore_exp=ignore_exp)
            for cfeature, cvals in feature_annotations.items():
                if cfeature not in feature_terms:
                    feature_terms[cfeature] = []
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
        fs_array = np.zeros([exp.data.shape[0], len(terms)])

        data = exp.get_data(sparse=False, copy=True)
        # how to weight the frequency of each bacteria
        if freq_weight == 'log':
            data[data < 1] = 1
            data = np.log2(data)
        elif freq_weight == 'binary':
            data = data > 0
        elif freq_weight == 'linear':
            pass
        elif freq_weight == 'rank':
            logger.debug('ranking the data')
            for ccol in range(np.shape(data)[1]):
                data[:, ccol] = scipy.stats.rankdata(data[:, ccol])

        else:
            raise ValueError('unknown freq_weight option %s. Can use ["log", "binary","linear"].')

        # iterate over all features and add to all terms associated with the feature
        for idx, cfeature in enumerate(exp_features):
            fterms = feature_terms[cfeature]
            for cterm, cval in fterms:
                fs_array[:, terms[cterm]] += cval * data[:, idx]

        # create the new experiment with samples x terms
        sm = deepcopy(exp.sample_metadata)
        sorted_term_list = sorted(terms, key=terms.get)
        fm = pd.DataFrame(data={'term': sorted_term_list}, index=sorted_term_list)
        fm['num_features'] = [term_features[d] for d in fm.index]
        newexp = Experiment(fs_array, sample_metadata=sm, feature_metadata=fm, description='Term scores')

        # get the differentially abundant terms between the two sample groups
        dd = newexp.diff_abundance(field, value1, value2, fdr_method=fdr_method, transform='log2data', alpha=alpha)
        return dd

    def draw_wordcloud(self, exp=None, features=None, term_type='fscore', ignore_exp=None, width=2000, height=1000, freq_weighted=False, relative_scaling=0.5, focus_terms=None, threshold=None, max_id=None):
        '''Draw a word_cloud for a given set of sequences

        Parameters
        ----------
        exp: calour.Experiment or None, optional
            The experiment containing the sequences (features) of interest.
            None to not ise the experiment (get the annotations for the features supplied from dbbact and use them)
        features: list of str or None, optional
            None to use the features from exp. Otherwise, a list of features ('ACGT' sequences) that is a subset of the features in exp (if exp is supplied)
            Note: if exp is None, must provide features.
        term_type: str, optional
            What score to use for the word cloud. Options are:
            'recall': sizes are based on the recall (fraction of dbbact term containing annotations that contain the sequences)
            'precision': sizes are based on the precision (fraction of sequence annotations of the experiment sequences that contain the term)
            'fscore': a combination of recall and precition (r*p / (r+p))
        ignore_exp: list of int or None or True, optional
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment (if exp is supplied)
            None (default) to use annotations from all experiments including the current one
        width, height: int, optional
            The width and heigh of the figure (high values are slower but nicer resolution for publication)
            If inside a jupyter notebook, use savefig(f, dpi=600)
        freq_weighted: bool, optional
            Only when supplying exp
            True to weight each bacteria by it's mean frequency in exp.data
            NOT IMPLEMENTED YET!!!
        relative_scaling: float, optional
            the effect of the score on the word size. 0.5 is a good compromise (passed to wordcloud.Wordcloud())
        focus_terms: list of str or None, optional
            if not None, use only annotations containing all terms in focus_terms list.
            for example, if focus_terms=['homo sapiens', 'feces'], will only draw the wordcloud for annotations of human feces
        threshold: float or None, optional
            if not None, show in word cloud only terms with p-value <= threshold (using the null model of binomial with term freq. as observed in all dbbact)
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id


        Returns
        -------
        matplotlib.figure
        '''
        import matplotlib.pyplot as plt
        try:
            from wordcloud import WordCloud
        except Exception as err:
            print(err)
            raise ValueError("Error importing wordcloud module. Is it installed? If not, install it using: pip install git+git://github.com/amueller/word_cloud.git")

        # get the annotations
        if exp is not None:
            if '__dbbact_sequence_terms' not in exp.exp_metadata:
                # if annotations not yet in experiment - add them
                self.add_all_annotations_to_exp(exp, max_id=max_id)
            # and filter only the ones relevant for features
            sequence_terms = exp.exp_metadata['__dbbact_sequence_terms']
            sequence_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
            annotations = exp.exp_metadata['__dbbact_annotations']
            term_info = exp.exp_metadata['__dbbact_term_info']
            if features is None:
                features = exp.feature_metadata.index.values
        else:
            if features is None:
                raise ValueError('Must supply experiment or features to use for wordcloud')
            sequence_terms, sequence_annotations, annotations, term_info, taxonomy = self.db.get_seq_list_fast_annotations(features)

        # filter based on focus_terms
        if focus_terms is not None:
            focus_terms = set(focus_terms)
            ok_annotations = {}
            for cid, cannotation in annotations.items():
                # check if an
                found_terms = set()
                for cdetail in cannotation['details']:
                    if cdetail[1] in focus_terms:
                        found_terms.add(cdetail[1])
                if len(found_terms) == len(focus_terms):
                    ok_annotations[cid] = cannotation
            logger.info('keeping %d out of %d annotations with all the terms (%s)' % (len(ok_annotations), len(annotations), focus_terms))
            print('keeping %d out of %d annotations with all the terms (%s)' % (len(ok_annotations), len(annotations), focus_terms))
            for k, v in sequence_annotations.items():
                nv = [x for x in v if x in ok_annotations]
                sequence_annotations[k] = nv

        # change the sequence annotations from dict to list of tuples
        sequence_annotations = [(k, v) for k, v in sequence_annotations.items()]

        # set the experiments to ignore in the wordcloud
        if ignore_exp is True:
            if exp is None:
                raise ValueError('Cannot use ignore_exp=True when exp is not supplied')
            ignore_exp = self.db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['sample_metadata_md5'], getall=True)
            if ignore_exp is None:
                logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
            else:
                logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))
        if ignore_exp is None:
            ignore_exp = []

        # we need to rekey the annotations with an str (old problem...)
        annotations = {str(k): v for k, v in annotations.items()}

        # calculate the recall, precision, fscore for each term
        fscores, recall, precision, term_count, reduced_f = get_enrichment_score(annotations, sequence_annotations, ignore_exp=ignore_exp, term_info=term_info, threshold=threshold)

        logger.debug('draw_cloud for %d words' % len(fscores))
        if len(fscores) == 0:
            logger.info('no words for wordcloud')
            return ''

        if term_type == 'fscore':
            score = fscores
        elif term_type == 'recall':
            score = recall
        elif term_type == 'precision':
            score = precision
        else:
            raise ValueError('term_type %s not recognized. options are: fscore, recall, precision')

        # weight by data frequency if needed
        if freq_weighted:
            new_score = defaultdict(float)
            for cseq in sequence_annotations:
                pass

        # normalize the fractions to a scale max=1
        new_scores = {}
        if score is not None:
            maxval = max(score.values())
            logger.debug('normalizing score. maxval is %f' % maxval)
            for ckey, cval in score.items():
                new_scores[ckey] = score[ckey] / maxval
        score = new_scores

        wc = WordCloud(width=width, height=height, background_color="white", relative_scaling=relative_scaling, stopwords=set(), color_func=lambda *x, **y: _get_color(*x, **y, fscore=score, recall=recall, precision=precision, term_count=term_count))
        # wc = WordCloud(width=400 * 3, height=200 * 3, background_color="white", relative_scaling=0.5, stopwords=set(), color_func=lambda *x, **y: _get_color(*x, **y, fscore=fscores, recall=recall, precision=precision, term_count=term_count))
        # else:
        #     wc = WordCloud(background_color="white", relative_scaling=0.5, stopwords=set(), color_func=lambda *x, **y: _get_color(*x, **y, fscore=fscores, recall=recall, precision=precision, term_count=term_count))

        if isinstance(score, str):
            logger.debug('generating from words list')
            wordcloud = wc.generate(score)
        elif isinstance(score, dict):
            logger.debug('generating from frequency dict')
            wordcloud = wc.generate_from_frequencies(score)
        else:
            logger.info('unknown type for generate_wordcloud!')

        fig = plt.figure()
        plt.imshow(wordcloud)
        plt.axis("off")
        fig.tight_layout()
        return fig

    def get_enrichment_score(self, *kargs, **kwargs):
        '''Get f score, recall and precision for set of annotations

        Parameters
        ----------
        annotations: dict of {annotationid (str): annotation(dict)}
        seqannotations: list of (seqid, [annotation ids])
        ingore_exp: list of str, optional
            list of experiment ids to ignore in the analysis
        term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
            The statistics about each term. if None, the function will contact dbbact to get the term_info
        term_types: list of str
            the terms to calculate enrichment scores for. options are:
            'single': each single term (i.e. 'feces')
            'pairs': all term pairs (i.e. 'feces+homo sapiens')
        threshold: float or None, optional
            if not None, return only terms that are significantly enriched in the annotations compared to complete database null with p-val <= threshold

        Returns
        -------
        fscore: dict of {term(str): fscore(float)}
        recall: dict of {term(str): recall(float)}
        precision: dict of {term(str): precision(float)}
        term_count: dict of {term(str): total experiments(float)}
            the number of experiments where annotations for each term appear
        reduced_f
        '''
        return get_enrichment_score(*kargs, **kwargs)

    def show_enrichment_qt5(self, group1, group2=None, exp=None, max_id=None, group1_name=None, group2_name=None, **kwargs):
        '''Show enriched terms between group1 and group2 using a qt5 GUI

        The gui shows a list of terms enriched in the 2 groups, and enables plotting per term heatmap and venn diagram.

        Parameters
        ----------
        cdb: DBBact.DBAccess
            the database interface for the analysis
        group1: list of str ('ACGT')
            the first group of sequences
        group2: list of str or None, optional
            the second group of sequences
            if None, group2 is taken from exp by using all sequences not in group1
        exp: calour.Experiment or None, optional
            the experiment on which the analysis is performed.
            If not None, annotations are taken from the experiment (if available)
            NOTE: if annotations are not available, they will be added to the experiment the first time the function is called
            If None, annotations are always queried from dbbact from each group and not stored for further queries.
        max_id: int or None, optional
            if not None, limit results to annotation ids <= max_id
        **kwargs: additional parameters supplied to db_access.term_enrichment(). These include:
        term_type : str or None (optional)
            The type of annotation data to test for enrichment
            options are:
            'term' - ontology terms associated with each feature.
            'parentterm' - ontology terms including parent terms associated with each feature.
            'annotation' - the full annotation strings associated with each feature
            'combined' - combine 'term' and 'annotation'
        ignore_exp: list of int or None or True, optional
            List of experiments to ignore in the analysis
            True to ignore annotations from the current experiment
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
        random_seed: int or None
            int to specify the random seed for numpy.random.
        use_term_pairs: bool, optional
            True to also test enrichment in pairs of terms (i.e. homo sapiens+feces, etc.)
        focus_terms: list of str or None, optional
            if not None, use only annotations containing all the terms in focus_terms
        '''
        from .enrichment_gui import show_enriched_terms_qt5
        show_enriched_terms_qt5(cdb=self, group1=group1, group2=group2, exp=exp, max_id=max_id, group1_name=group1_name, group2_name=group2_name, **kwargs)


def _get_color(word, font_size, position, orientation, font_path, random_state, fscore, recall, precision, term_count):
    '''Get the color for a wordcloud term based on the term_count and higher/lower

    If term starts with "-", it is lower in and colored red. otherwise colored blue
    If we have term_count, we use it to color from close to white(low count) to full color (>=10 experiments)

    Parameters
    ----------
    fscores: dict of {term(str): fscore(float)}
        between 0 and 1
    recall: dict of {term(str): recall score(float)}, optional
    precision: dict of {term(str): precision score(float)}, optional
    term_count: dict of {term(str): number of experiments with term(float)}, optional
        used to determine the color intensity


    Returns
    -------
    str: the color in hex "0#RRGGBB"
    '''
    import matplotlib as mpl
    if word in term_count:
        count = min(term_count[word], 10)
    else:
        count = 10

    if word[0] == '-':
        cmap = mpl.cm.get_cmap('Oranges')
        rgba = cmap(float(0.7 + count / 40), bytes=True)
    else:
        cmap = mpl.cm.get_cmap('Purples')
        rgba = cmap(float(0.6 + count / 40), bytes=True)

    red = format(rgba[0], '02x')
    green = format(rgba[1], '02x')
    blue = format(rgba[2], '02x')
    return '#%s%s%s' % (red, green, blue)
