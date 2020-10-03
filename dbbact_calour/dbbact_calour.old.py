'''
dbbact (:mod:`dbbact_calour.dbbact`)
====================================

.. currentmodule:: dbbact_calour.dbbact

Functions
^^^^^^^^^
.. autosummary::
   :toctree: generated

   DBBact.enrichment
   DBBact.add_annotations
   DBBact.add_all_annotations_to_exp
   DBBact.add_annotation
   DBBact.get_feature_terms
   DBBact.delete_annotation
   DBBact.upadte_annotation
   DBBact.show_term_details
   DBBact.plot_term_annotations
   DBBact.sample_enrichment
   DBBact
'''

from .db_access import DBAccess
from collections import defaultdict
from logging import getLogger, NOTSET, basicConfig
from logging.config import fileConfig
from copy import deepcopy
from pkg_resources import resource_filename

import numpy as np
import pandas as pd
import scipy.stats

from calour.util import get_config_value
from calour.database import Database
from calour.experiment import Experiment

logger = getLogger(__name__)

try:
    # get the logger config file location
    log = resource_filename(__name__, 'log.cfg')
    # log = path.join(path.dirname(path.abspath(__file__)), 'log.cfg')
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


class DBBact(Database):
    def __init__(self, exp=None):
        print('bobo')
        super().__init__(database_name='dbBact', methods=['get', 'annotate', 'enrichment'])
        username = get_config_value('username', section='dbbact')
        password = get_config_value('password', section='dbbact')
        self.db = DBAccess(username=username, password=password)

    def add_all_annotations_to_exp(self, exp, **kwargs):
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
        '__dbbact_taxonomy': dict of {sequence(str): taxonomy(str)}
            contains the dbbact assigned taxonomy for each feature

        Parameters
        ----------
        exp : ``Experiment``
            The experiment to get the details for and store them in
        **kwargs:
            extra parameters to pass to get_seq_list_fast_annotations()

        Returns:
        str
        '' if ok, otherwise error string
        '''
        print('bibi')
        logger.debug('Getting annotations for %d sequences' % len(exp.feature_metadata))
        sequence_terms, sequence_annotations, annotations, term_info, taxonomy = self.db.get_seq_list_fast_annotations(exp.feature_metadata.index.values, **kwargs)
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

    def get_feature_terms(self, features, exp=None, term_type=None, ignore_exp=None, **kwargs):
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
        **kwargs: extra parameters to pass to db_access.get_seq_list_fast_annotations()

        Returns
        -------
        feature_terms : dict of list of str/int
            key is the feature, list contains all terms associated with the feature
        '''
        print('baba')
        if term_type is None:
            term_type = 'terms'
        if exp is not None:
            if '__dbbact_sequence_terms' not in exp.exp_metadata:
                # if annotations not yet in experiment - add them
                self.add_all_annotations_to_exp(exp, **kwargs)
            # and filter only the ones relevant for features
            sequence_terms = exp.exp_metadata['__dbbact_sequence_terms']
            sequence_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
            annotations = exp.exp_metadata['__dbbact_annotations']
        else:
            sequence_terms, sequence_annotations, annotations, term_info, taxonomy = self.db.get_seq_list_fast_annotations(features, **kwargs)
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
                    cdesc = self.db.get_annotation_string(annotations[cannotation])
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

    def enrichment(self, exp, features, **kwargs):
        '''Get the list of enriched terms in features compared to all features in exp.

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
        exp_features = set(exp.feature_metadata.index.values)
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

        res = self.db.term_enrichment(g1_features=features, g2_features=bg_features, all_annotations=exp.exp_metadata['__dbbact_annotations'], seq_annotations=exp.exp_metadata['__dbbact_sequence_annotations'], **kwargs)
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

    def show_term_details(self, term, exp, features, group2_features, group1_name='group1', group2_name='group2', gui='qt5'):
        '''
        Plot a heatmap for all annotations containing term in experiment
        Rows are the annotations, columns are the sequences (sorted by features/group2_features)
        '''
        if term[0] == '-':
            term = term[1:]
        all_seqs = features.copy()
        all_seqs.extend(group2_features)
        tmat, tanno, tseqs = self.db.get_term_annotations(term, all_seqs, exp.exp_metadata['__dbbact_sequence_annotations'], exp.exp_metadata['__dbbact_annotations'])
        seq_group = [str(group1_name)] * len(features)
        seq_group.extend([str(group2_name)] * len(group2_features))
        tseqs['group'] = seq_group
        newexp = Experiment(tmat, sample_metadata=tseqs, feature_metadata=tanno)
        newexp = newexp.cluster_features(1)
        newexp = newexp.sort_by_metadata(field='expid', axis='f')
        newexp.plot(feature_field='annotation', gui=gui, yticklabel_kwargs={'rotation': 0}, yticklabel_len=35, cmap='tab20b', norm=None, bary_fields=['expid'], bary_label=True, barx_fields=['group'], barx_label=True)

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

    def plot_term_annotations_venn(self, term, exp, bacteria_groups=None, min_prevalence=0, annotation_types=None, set_colors=('mediumblue', 'r', 'g')):
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

        Returns
        -------
        list of figures
        '''
        from matplotlib_venn import venn3

        # get the two bacteria groups
        if bacteria_groups is None:
            vals = exp.feature_metadata['_calour_diff_abundance_group'].unique()
            group1 = list(exp.feature_metadata.index[exp.feature_metadata['_calour_diff_abundance_group'] == vals[0]].values)
            group1_name = vals[0]
            group2 = list(exp.feature_metadata.index[exp.feature_metadata['_calour_diff_abundance_group'] == vals[1]].values)
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
            f = plt.figure()
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
            # venn3([gg1, gg2, anno_group], set_labels=(group1_name, group2_name, 'annotation'), set_colors=['mediumblue', 'r', 'g'])
            venn3(vvals, set_labels=(group1_name, group2_name, 'annotation'), set_colors=set_colors)
            plt.title('%s (expid %s)' % (canno['annotation'], canno['expid']))
            all_figures.append(f)
        return all_figures

    def sample_enrichment(self, exp, field, value1, value2=None, term_type='term', ignore_exp=None, min_appearances=3, fdr_method='dsfdr', score_method='all_mean', freq_weight='log', alpha=0.1, use_term_pairs=False):
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
            self.add_all_annotations_to_exp(exp)

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
