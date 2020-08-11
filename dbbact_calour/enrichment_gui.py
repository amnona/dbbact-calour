from logging import getLogger

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QVBoxLayout, QListWidget, QDialogButtonBox, QHBoxLayout, QLabel, QPushButton
import numpy as np

from .dbannotation import init_qt5

logger = getLogger(__name__)


def show_enriched_terms_qt5(cdb, group1, group2=None, exp=None, max_id=None, group1_name=None, group2_name=None, **kwargs):
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
    if group1_name is None:
        group1_name = 'Group1'
    if exp is None:
        logger.debug('No experiment supplied - getting annotations')
        # need to write
        all_features = set(group1).union(set(group2))
        # get all the annotations....
        # all_annotations = exp.exp_metadata['__dbbact_annotations']
        # seq_annotations = exp.exp_metadata['__dbbact_sequence_annotations']
        # term_info = exp.exp_metadata.get('__dbbact_term_info')
        raise ValueError('exp=None not supported yet.')
    else:
        # add all annotations to experiment if not already added
        cdb.add_all_annotations_to_exp(exp, max_id=max_id, force=False)
        all_annotations = exp.databases['dbbact']['annotations']
        seq_annotations = exp.databases['dbbact']['sequence_annotations']
        term_info = exp.databases['dbbact'].get('term_info')

        # prepare group2 by taking all features from exp that are not in group1
        if group2 is None:
            logger.debug('Using experiment sequences for group2')
            exp_features = set(exp.feature_metadata.index.values)
            group2 = np.array(list(exp_features.difference(group1)))
            if group2_name is None:
                group2_name = 'NOT %s' % group1_name

        if group2_name is None:
            group2_name = 'Group2'

        # validate all features are in exp. otherwise ignore them
        exp_features = set(exp.feature_metadata.index.values)
        all_features = set(group1).union(set(group2))
        bad_features = all_features.difference(exp_features)
        if len(bad_features) > 0:
            logger.warning('Some of the features for enrichment are not in the experiment. Ignoring %d features' % len(bad_features))
            group1 = list(set(group1).intersection(exp_features))
            group2 = list(set(group2).intersection(exp_features))
            if len(group1) == 0:
                raise ValueError("No features left in group1 after ignoring. Please make sure you test for enrichment with features from the experiment.")
            if len(group2) == 0:
                raise ValueError("No features left in group2 after ignoring. Please make sure you test for enrichment with features from the experiment.")

    # if ignore exp is True, it means we should ignore the current experiment
    ignore_exp = kwargs.get('ignore_exp')
    if ignore_exp is True:
        ignore_exp = cdb.db.find_experiment_id(datamd5=exp.info['data_md5'], mapmd5=exp.info['sample_metadata_md5'], getall=True)
        if ignore_exp is None:
            logger.warn('No matching experiment found in dbBact. Not ignoring any experiments')
        else:
            logger.info('Found %d experiments (%s) matching current experiment - ignoring them.' % (len(ignore_exp), ignore_exp))
    kwargs['ignore_exp'] = ignore_exp

    # get the enriched terms (pandas dataframe)
    enriched, resmat, features = cdb.db.term_enrichment(g1_features=group1, g2_features=group2, all_annotations=all_annotations, seq_annotations=seq_annotations, term_info=term_info, **kwargs)
    logger.debug('Got %d enriched terms' % len(enriched))

    app, app_created = init_qt5()

    if len(enriched) == 0:
        QtWidgets.QMessageBox.information(None, "No enriched terms found",
                                          "No enriched annotations found when comparing\n%d selected sequences to %d "
                                          "other sequences" % (len(group1), len(group2)))
        return

    # sort by absolute value of effect size, so terms appear from high to low in both lists
    enriched['odif_abs'] = enriched['odif'].abs()
    enriched = enriched.sort_values('odif_abs', ascending=False)

    # add the terms to the appropriate list
    listwin = TermInfoListWindow(listname='enriched ontology terms', group1name=group1_name, group2name=group2_name)
    for idx, cres in enriched.iterrows():
        if cres['odif'] > 0:
            ccolor = 'blue'
            cgroup = 1
        else:
            ccolor = 'red'
            cgroup = 2
        cname = cres['term']
        # For each enriched term, double clicking will display a heatmap
        # where all annotations containing the term are the features,
        # and bacteria (from the two compared groups) are the samples.
        # This enables seeing where does the enrichment for this term come from.
        # i.e. which bacteria are present in each annotation containing this term.
        dblclick_data = {}
        dblclick_data['database'] = cdb
        dblclick_data['term'] = cname
        dblclick_data['exp'] = exp
        g1_seqs = set(group1)
        g2_seqs = set(group2)
        ordered_g1_seqs = [s for s in exp.feature_metadata.index.values[::-1] if s in g1_seqs]
        ordered_g2_seqs = [s for s in exp.feature_metadata.index.values[::-1] if s in g2_seqs]
        dblclick_data['features1'] = ordered_g1_seqs
        dblclick_data['features2'] = ordered_g2_seqs
        listwin.add_item('%s - effect %f, pval %f ' % (cname, cres['odif'], cres['pvals']), color=ccolor, dblclick_data=dblclick_data, group=cgroup)

    listwin.exec_()


class TermInfoListWindow(QtWidgets.QDialog):
    def __init__(self, group1data=[], group2data=[], listname=None, group1name=None, group2name=None):
        '''Create a list window with items in the list and the listname as specified

        Parameters
        ----------
        listdata: list of str, optional
            the data to show in the list
        listname: str, optional
            name to display above the list
        '''
        super().__init__()
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        if listname is not None:
            self.setWindowTitle(listname)

        self.layout = QVBoxLayout(self)

        listlayout = QHBoxLayout()
        g1layout = QVBoxLayout()
        g2layout = QVBoxLayout()
        self.w_list = QListWidget()
        self.w2_list = QListWidget()
        if group1name is None:
            group1name = 'group1'
        if group2name is None:
            group1name = 'group2'
        g1layout.addWidget(QLabel('higher in %s' % group1name))
        g2layout.addWidget(QLabel('higher in %s' % group2name))
        g1layout.addWidget(self.w_list)
        g2layout.addWidget(self.w2_list)

        self.group1name = group1name
        self.group2name = group2name

        listlayout.addLayout(g1layout)
        listlayout.addLayout(g2layout)

        self.layout.addLayout(listlayout)

        buttonlayout = QHBoxLayout()
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonlayout.addWidget(buttonBox)
        venn_button = QPushButton('Venn')
        venn_button.clicked.connect(self.venn)
        buttonlayout.addWidget(venn_button)
        heatmap_button = QPushButton('Term Heatmap')
        heatmap_button.clicked.connect(self.heatmap_clicked)
        buttonlayout.addWidget(heatmap_button)
        self.layout.addLayout(buttonlayout)
        self.layout.addWidget(buttonBox)

        for citem in group1data:
            self.w_list.addItem(citem)
        for citem in group2data:
            self.w2_list.addItem(citem)

        self.w_list.itemDoubleClicked.connect(self.list_double_click)
        self.w2_list.itemDoubleClicked.connect(self.list_double_click)

        self.w_list.currentItemChanged.connect(self.selection_change)
        self.w2_list.currentItemChanged.connect(self.selection_change)

        self.cselection = None

        self.show()
        self.adjustSize()

    def selection_change(self, current, previous):
        old_selection = self.cselection
        if old_selection is not None:
            old_selection.setSelected(False)
        self.cselection = current

    def venn(self):
        if self.cselection is None:
            logger.info('Must select term first')
            return
        data = self.cselection.data(QtCore.Qt.UserRole)
        cterm = data['term']
        if cterm.startswith('LOWER IN '):
            cterm = '-' + cterm[len('LOWER IN '):]
        f = data['database'].plot_term_venn_all(cterm, data['exp'], set_colors=('red', 'green', 'mediumblue'), max_size=500, ignore_exp=True)
        f.show()
        # print(data)
        # plot_term_venn_all(self, terms, exp, bacteria_groups=None, set_colors=('red', 'green', 'mediumblue'), max_size=None, ignore_exp=[]):

    def heatmap_clicked(self):
        if self.cselection is None:
            logger.info('Must select term first')
            return
        data = self.cselection.data(QtCore.Qt.UserRole)
        self.plot_heatmap(data)

    def plot_heatmap(self, data):
        data = self.cselection.data(QtCore.Qt.UserRole)
        data['database'].show_term_details(data['term'], data['exp'], data['features1'], data['features2'], gui='qt5', group1_name=self.group1name, group2_name=self.group2name, title='Annotations for term %s' % data['term'])

    def add_item(self, text, color='black', dblclick_data=None, group=1):
        '''Add an item to the list

        Parameters
        ----------
        text : str
            the string to add
        color : str, optional
            the color of the text to add
        dblclick_function : function or None
            the function to call when this item is double clicked (or None to ignore)
        '''
        item = QtWidgets.QListWidgetItem()
        item.setText(text)
        if color == 'black':
            ccolor = QtGui.QColor(0, 0, 0)
        elif color == 'red':
            ccolor = QtGui.QColor(155, 0, 0)
        elif color == 'blue':
            ccolor = QtGui.QColor(0, 0, 155)
        elif color == 'green':
            ccolor = QtGui.QColor(0, 155, 0)
        item.setForeground(ccolor)
        item.setData(QtCore.Qt.UserRole, dblclick_data)
        if group == 1:
            self.w_list.addItem(item)
        else:
            self.w2_list.addItem(item)

    def list_double_click(self, item):
        data = item.data(QtCore.Qt.UserRole)
        self.plot_heatmap(data)
