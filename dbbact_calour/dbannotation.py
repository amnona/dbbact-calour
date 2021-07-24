from logging import getLogger
from pkg_resources import resource_filename
import pickle
import sys
import os
import os.path

from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QListWidget, QDialogButtonBox

import pandas as pd
from calour.util import get_config_value, set_config_value, get_config_file

from . import dbbact

logger = getLogger(__name__)

# the global variable to store annotation history
# for autofill of 'ALL' details
# key is the md5 of the data, value is the annotations
history = {}


def init_qt5():
    '''Init the qt5 event loop

    Parameters
    ----------

    Returns
    -------
    app :
        QCoreApplication
    app_created : bool
        True if a new QApplication was created, False if using existing one
    '''
    app_created = False
    app = QtCore.QCoreApplication.instance()
    logger.debug('Qt app is %s' % app)
    if app is None:
        # app = QApplication(sys.argv)
        app = QApplication(sys.argv)
        app_created = True
        logger.debug('Qt app created')
    if not hasattr(app, 'references'):
        app.references = set()

    return app, app_created


def annotate_bacteria_gui(dbclass, seqs, exp):
    '''Create a dialog for annotating the bacteria into dbbact

    Parameters
    ----------
    dbclass : dbbact.DBBact
        the dbBact database to interact with
    seqs : list of sequences ('ACGT')
        the sequences to add to the database
    exp : Experiment
        Experiment containing the sample and experiment metadata (i.e. md5 etc.)

    Returns
    -------
    str
        empty if ok, error details str if error enoucntered
    '''
    app, app_created = init_qt5()

    # test if we have user/password set-up
    test_user_password(dbclass.db)

    # test if study already in database
    cdata = dbclass.db.find_experiment_id(datamd5=exp.info['data_md5'], mapmd5=exp.info['sample_metadata_md5'])
    if cdata is None:
        logger.info('Study does not exist in dbbact. Creating new study')
        cdata = study_data_ui(exp)
        if cdata is None:
            msg = 'no study information added. Please add at least one field. Annotation cancelled'
            logger.warn(msg)
            return msg

    show_and_ask('Please enter annotations for the study.\n'
                 'Choose the annotation type. Annotations can be:\n'
                 '"Differential abundance": High in one group compared to the other\n'
                 '"Common": present in majority of the samples (prevalence)\n'
                 '"Dominant": frequency of the bacteria is >1% (mean frequency)\n\n'
                 'Then choose the terms (preferably from ontology autocomplete) relevant\n'
                 'Note that for "Differential abundance", you need to choose the terms common to both groups (all),'
                 'the terms in the group where the bacteria is higher (high) and where it is lower (low).\n'
                 'For example, a fecal bacteria that is more common in females will be annotated as:\n'
                 '"all": feces, homo sapience\n'
                 '"low": male\n'
                 '"high": female\n'
                 'also, please supply as many terms as possible (host, material, country, disease, etc.)',
                 keyval='annotation_info')

    ui = DBAnnotateSave(seqs, exp, dbclass=dbclass, get_seq_region=True)
    res = ui.exec_()
    if res == QtWidgets.QDialog.Accepted:
        description = str(ui.bdescription.text())
        # TODO: need to get primer region!!!!
        primerid = ui.primerid
        # primerid = 'V4'
        method = str(ui.bmethod.text())
        if method == '':
            method = 'na'
        # TODO: fix submitter name
        submittername = 'Amnon Amir'
        annotations = []

        for citem in qtlistiteritems(ui.blistall):
            cdat = qtlistgetdata(citem)
            cval = cdat['value']
            cterm_id = cdat['term_id']
            ctype = cdat['type']
            # convert synonyms to original ontology terms
            if cval in DBAnnotateSave._ontology_from_id:
                cval = DBAnnotateSave._ontology_from_id[cval]
            else:
                logger.debug("item %s not found in ontologyfromid" % cval)

            # if the user selected a specific term_id, use it
            if cterm_id != '':
                annotations.append((ctype, cterm_id))
            # otherwise, use just the term description
            else:
                annotations.append((ctype, cval))
        # get annotation type
        if ui.bdiffpres.isChecked():
            annotation_type = 'DIFFEXP'
        elif ui.bisa.isChecked():
            curtypeval = ui.bisatype.currentText()
            if 'Common' in curtypeval:
                annotation_type = 'COMMON'
            elif 'Contam' in curtypeval:
                annotation_type = 'CONTAMINATION'
            elif 'Dominant' in curtypeval:
                annotation_type = 'DOMINANT'
            elif 'Positively associated' in curtypeval:
                annotation_type = 'POSITIVE ASSOCIATION'
            elif 'Negatively associated' in curtypeval:
                annotation_type = 'NEGATIVE ASSOCIATION'
            else:
                annotation_type = 'OTHER'
        else:
            msg = "No annotation type selected. Annotation not saved"
            logger.warn(msg)
            return msg

        logger.debug('Adding annotation to studyid %s' % cdata)
        err, annotation_id = dbclass.db.add_annotations(expid=cdata, sequences=seqs, annotationtype=annotation_type, annotations=annotations, submittername=submittername, description=description, method=method, primerid=primerid)
        if err:
            msg = 'Annotation not added. error: %s' % err
            logger.warn(msg)
            QtWidgets.QMessageBox.warning(None, 'Annotation not added', msg)
            return msg
        logger.debug('New annotation added. AnnotationId=%d' % annotation_id)
        err, annotation_res = dbclass.db.get_annotation(annotation_id)
        if err == '':
            # store the history of the annotation for the next annotation of the same experiment (so description, details etc. will auto fill)
            global history
            history[exp.info['data_md5']] = annotation_res
        else:
            logger.warn('Failed to get annotation for annotationID %d. please validate annotation was added to dbBact.' % annotation_id)
        return ''
    return 'Add annotation cancelled'


def update_annotation_gui(dbclass, annotation, exp):
    '''Update an existing annotation

    Parameters
    ----------
    db: dbbact.DBBact
        the interface to dbbact (for storing autocomplete)
    annotation: dict
        the annotation to update (dict including the key "annotationid")
    exp: Calour.Experiment
        the experiment on which to update the annotation
    '''
    app, app_created = init_qt5()

    # test if we have user/password set-up
    test_user_password(dbclass.db)

    annotationid = annotation['annotationid']
    primerid = annotation.get('primerid')

    ui = DBAnnotateSave([], exp, prefill_annotation=annotation, dbclass=dbclass)
    res = ui.exec_()
    if res == QtWidgets.QDialog.Accepted:
        primerid = ui.primerid
        description = str(ui.bdescription.text())
        method = str(ui.bmethod.text())
        if method == '':
            method = 'na'
        # TODO: fix submitter name
        annotations = []

        for citem in qtlistiteritems(ui.blistall):
            cdat = qtlistgetdata(citem)
            cval = cdat['value']
            ctype = cdat['type']
            # convert synonyms to original ontology terms
            if cval in DBAnnotateSave._ontology_from_id:
                cval = DBAnnotateSave._ontology_from_id[cval]
            else:
                logger.debug("item %s not found in ontologyfromid" % cval)
            annotations.append((ctype, cval))
        # get annotation type
        if ui.bdiffpres.isChecked():
            annotation_type = 'DIFFEXP'
        elif ui.bisa.isChecked():
            curtypeval = ui.bisatype.currentText()
            if 'Common' in curtypeval:
                annotation_type = 'COMMON'
            elif 'Contam' in curtypeval:
                annotation_type = 'CONTAMINATION'
            elif 'Dominant' in curtypeval:
                annotation_type = 'DOMINANT'
            elif 'Positively associated' in curtypeval:
                annotation_type = 'POSITIVE ASSOCIATION'
            elif 'Negatively associated' in curtypeval:
                annotation_type = 'NEGATIVE ASSOCIATION'
            else:
                annotation_type = 'OTHER'
        else:
            msg = "No annotation type selected. Annotation not saved"
            logger.warn(msg)
            return msg
        logger.debug('Updating annotations for annotationid %d' % annotationid)
        res = dbclass.db.send_update_annotation(annotationid=annotationid, annotationtype=annotation_type, annotations=annotations, description=description, method=method)
        if res is None:
            msg = 'Annotation not updated.'
            logger.warn(msg)
            return msg
        logger.debug('Annotation updated.')

        err, annotation_res = dbclass.db.get_annotation(annotationid)
        if err == '':
            # store the history of the annotation for the next annotation of the same experiment (so description, details etc. will auto fill)
            global history
            history[exp.info['data_md5']] = annotation_res
        else:
            logger.warn('Failed to get annotation for annotationID %d. please validate annotation was added to dbBact.' % annotationid)
        return ''
    return 'Update annotation cancelled'


def test_user_password(db):
    '''
    Test if the config file has user/password and if not ask for one.
    Also don't ask if the 'ask again' checkbox is checked

    Parameters
    ----------
    db : DBBact.DBAccess
        the database interface class.
    '''
    logger.debug('Testing if user/pwd in config file')
    username = get_config_value('username', section='dbbact')
    if username is not None:
        logger.debug('found user %s' % username)
        return
    if get_config_value('show_user_request', section='dbbact') is not None:
        logger.debug('user/password not set, but show_user_request flag in config file is set, so ignoring')
        return
    logger.debug('no username in config file')
    if QtWidgets.QMessageBox.warning(None, 'Register/login for better annotation', 'You can add annotations as an anonymous user, '
                                           'or register/login in order to create non-anonymous annotations.\n'
                                           'NOTE: annonymous annotations can be modified/deleted by any user, '
                                           'whereas non-annonymous annotations can be modified/deleted only by the user who created them.\n\n'
                                           'Click "Yes" to register/login or "No" to continue as annonymous user',
                                     QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.No:
        return
    ui = UserPasswordDlg()
    res = ui.exec_()
    if res == QtWidgets.QDialog.Accepted:
        username = str(ui.username.text())
        password = str(ui.password.text())
        db.username = username
        db.password = password
        userid = db.get_user_id()
        if userid is None:
            logger.debug('userid/password not found. Registering')
            email = str(ui.email.text())
            description = str(ui.interest.currentText())
            err = db.register_user(username, password, email=email, description=description, publish='n', name='')
            if err:
                logger.warn(err)
                QtWidgets.QMessageBox.warning(None, 'Login failed', 'login for user %s failed.\n'
                                                    'You are now logged in as anonymous user.' % username)
                return
            QtWidgets.QMessageBox.information(None, 'Registered new user', 'New user %s registered in dbBact.\n'
                                                    'You can change the Calour dbBact user/password in the config file:\n'
                                                    '%s' % (username, get_config_file()))
        else:
            userid = db.get_user_id()
            if userid is not None:
                QtWidgets.QMessageBox.information(None, 'Logged in existing user', 'user %s logged into dbBact.\n'
                                                        'You can change the Calour dbBact user/password in the config file:\n'
                                                        '%s' % (username, get_config_file()))
            else:
                QtWidgets.QMessageBox.warning(None, 'Login failed', 'login for user %s failed.\n'
                                                    'You are now logged in as anonymous user.' % username)
                return

        logger.info('storing username %s in config file' % username)
        set_config_value('username', username, section='dbbact')
        set_config_value('password', password, section='dbbact')


def study_data_ui(cexp):
    """
    open the study info window and show/get new references for the study data

    Parameters
    ----------
    cexp : Experiment
        the experiment for which to show the data (uses the datamd5 and mapmd5)

    Returns
    -------
    dataid : int or None
        studyid if study  has data, None if no data
    """
    bdb = dbbact.DBBact()
    cid = bdb.db.find_experiment_id(datamd5=cexp.info['data_md5'], mapmd5=cexp.info['sample_metadata_md5'])

    # add study details we have
    study_details = []
    if cid is None:
        study_details.append(('DataMD5', cexp.info['data_md5']))
        study_details.append(('MapMD5', cexp.info['sample_metadata_md5']))
        interesting_columns = {'sra_study_s': 'sra', 'project_name_s': 'name', 'experiment_title': 'name', 'experiment_design_description': 'name', 'BioProject_s': 'sra', 'run_date': 'date'}
        for ccol in cexp.sample_metadata.columns:
            if ccol.lower() in interesting_columns:
                if len(pd.unique(cexp.sample_metadata[ccol].ravel())) == 1:
                    study_details.append((interesting_columns[ccol.lower()], cexp.sample_metadata[ccol][0]))

    show_and_ask('Study is not in database.\n'
                 '(based on data or mapping file md5)\n'
                 'Please supply details that identify the source of the data\n'
                 'Preferably a dataID (such as sra/qiita/etc)'
                 ', the name of the study (name), and the DOID\n',
                 keyval='study_info')

    dbsi = DBStudyInfo(cid, study_details=study_details)
    res = dbsi.exec_()
    if res == QtWidgets.QDialog.Accepted:
        newstudydata = []
        allstudydata = []
        for citem in qtlistiteritems(dbsi.blist):
            cdata = qtlistgetdata(citem)
            allstudydata.append((cdata['type'], cdata['value']))
            if cdata['fromdb'] is False:
                newstudydata.append((cdata['type'], cdata['value']))

        if len(allstudydata) <= 2:
            logger.warn('not enough details. not saving anything')
            return None
        if len(newstudydata) == 0:
            logger.warn('No new items. not saving anything')
            return cid
        dataid = bdb.db.add_experiment_data(newstudydata, expid=cid)
        if dataid is None:
            logger.error('Error enountered - experiment not added')
            return None
        logger.debug('Study data saved to id %d' % dataid)
        return dataid
    return None


class DBStudyInfo(QtWidgets.QDialog):
    '''Show/get details about a study
    '''
    def __init__(self, experimentid=None, study_details=[]):
        '''Init the values in the study info dialog

        Parameters
        ----------
        experimentid : int
            the id pf the experiment in the database or None if new experiment
        exp_details : list of ('str','str')
            A list of tuples ('detail name','detail value') to add to the study
            usually init to the data_d5 and sample_metadata_md5
        '''
        super(DBStudyInfo, self).__init__()
        uic.loadUi(resource_filename(__package__, 'ui/studyinfo.ui'), self)
        self.dataid = 0
        bdb = dbbact.DBBact()

        # add experiment details from the database
        if experimentid is not None:
            info = bdb.get_experiment_info(experimentid)
            for cinfo in info:
                qtlistadd(self.blist, cinfo[0] + ':' + cinfo[1], {'fromdb': True, 'type': cinfo[0], 'value': cinfo[1]}, color='grey')
        else:
            # if not in database, get details from map file
            # the supplied details (supplied to the init function via exp_details)
            for cinfo in study_details:
                qtlistadd(self.blist, cinfo[0] + ':' + str(cinfo[1]), {'fromdb': False, 'type': cinfo[0], 'value': str(cinfo[1])}, color='black')

        self.bplus.clicked.connect(self.plus)
        self.bvalue.returnPressed.connect(self.plus)
        self.bminus.clicked.connect(self.minus)
        self.bannotations.clicked.connect(self.annotations)
        # self.bvalue.returnPressed.connect(self.plus)
        # self.cexp=expdat
        # self.setWindowTitle(self.cexp.studyname)
        self.bvalue.setFocus()

    def keyPressEvent(self, e):
        '''
        override the enter event so will not close dialog
        '''
        e.ignore()

    def addentry(self, fromdb, ctype, value, color='black'):
        if len(ctype) > 0 and len(value) > 0:
            newentry = '%s:%s' % (ctype, value)
            for citem in getqtlistitems(self.blist):
                if citem == newentry:
                    logger.debug('item already in list %s' % newentry)
                    return
            qtlistadd(self.blist, newentry, {'fromdb': False, 'type': ctype, 'value': value}, color="black")

    def plus(self):
        ctype = str(self.btype.currentText())
        cval = str(self.bvalue.text())
        self.addentry(fromdb=False, ctype=ctype, value=cval, color='black')
        self.bvalue.setText('')

    def minus(self):
        items = self.blist.selectedItems()
        for citem in items:
            cdata = qtlistgetdata(citem)
            if cdata['fromdb']:
                print('delete from db')
            self.blist.takeItem(self.blist.row(citem))

    def annotations(self):
        pass

    def inputCheck(self):
        '''Validate the new study information is good enough
        We link here from the accept slot

        Returns
        -------
        str
        '' (empty string) if details valid, non-empty error string if problem encountered
        '''
        # get all the details about the study
        newstudydata = []
        allstudydata = []
        name_found = False
        for citem in qtlistiteritems(self.blist):
            cdata = qtlistgetdata(citem)
            if cdata['type'].lower() == 'name':
                name_found = True
            allstudydata.append((cdata['type'], cdata['value']))
            if cdata['fromdb'] is False:
                newstudydata.append((cdata['type'], cdata['value']))
        # validate we have something else than the md5 for data/map
        if len(allstudydata) <= 2:
            QtWidgets.QMessageBox.warning(self, 'Missing study information', 'Please enter an identifier for the study (name recommended)')
            return 'Please enter an identifier for the study (name recommended)'
        if not name_found:
            if QtWidgets.QMessageBox.warning(self, 'Study information missing', '"name" field not supplied Do you want to continue?',
                                             QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.No:
                return 'Missing name field'
        return ''

    def accept(self):
        err = self.inputCheck()
        if err == '':
            self.done(1)
        else:
            logger.warn(err)


class DBAnnotateSave(QtWidgets.QDialog):
    '''The gui class for getting the annotations about a set of sequences
    '''

    # Attributes:
    # the ontologies preloaded for the autocomplete
    _ontologies_loaded = set()
    # used to store the ontology values for the autocomplete
    _ontology_dict = None
    # used to store the dbbact (user) ontology terms
    _ontology_dbbact = None
    # the maximal ontology term id present local (for syncing with the dbbact server for new terms)
    _ontology_dbbact_max_id = 0
    # the ontologies autocomplete sorted list
    _ontology_sorted_list = []

    def __init__(self, selected_features, expdat, prefill_annotation=None, dbclass=None, get_seq_region=False):
        '''Create the manual annotation window

        Parameters
        ----------
        selected_features : list of sequences ('ACGT')
            The sequences to annotate
        expdat : Experiment
            The experiment being annotated
        prefill_annotation : dict or None (optional)
            None (default) to prefill from history
            dict to pre-fill using dict fields instead
        dbclass: dbbact or None, optional
            the dbBact database to interact with (used to get the primer regions possible)
        get_seq_region: bool, optional
            if True, query dbbact for the sequences region if needed (i.e. no history) and prefill the region (i.e. 'v4')

        '''
        super(DBAnnotateSave, self).__init__()

        # create the gui
        uic.loadUi(resource_filename(__package__, 'ui/manualdata.ui'), self)
        self.bplus.clicked.connect(self.plus)
        self.bminus.clicked.connect(self.minus)
        self.bontoinput.returnPressed.connect(self.plus)

        # this is used to limit the autocompleter to only when enough chars are entered
        self.bontoinput.textEdited.connect(self.ontochanged)

        self.bstudyinfo.clicked.connect(self.studyinfo)
        self.bisa.toggled.connect(self.radiotoggle)
        self.bdiffpres.toggled.connect(self.radiotoggle)
        self.bisatype.currentIndexChanged.connect(self.isatypechanged)

        # set the primer region button
        if dbclass:
            self.bhistory.clicked.connect(lambda: self.select_primers(dbclass))

        self.cexp = expdat
        self.lnumbact.setText(str(len(selected_features)))

        # create the autocomplete ontology list
        completer = QtWidgets.QCompleter()
        self.bontoinput.setCompleter(completer)
        model = QtCore.QStringListModel()
        completer.setModel(model)
        # completer.setCompletionMode(QCompleter.InlineCompletion)
        completer.maxVisibleItems = 10
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        # make the completer selection also erase the text edit
        # completer.activated.connect(self.cleartext, type=Qt.QueuedConnection)
        # on selection from the completer, add to the ontology list
        completer.activated.connect(self.plus, type=Qt.QueuedConnection)
        # init the ontology values
        self._load_ontologies(dbclass)
        model.setStringList(self._ontology_sorted_list)
        # store the completer
        self._ontology_completer = completer

        self.setWindowTitle(expdat.description)

        self.primerid = None

        if prefill_annotation is None:
            # try to autofill from history if experiment annotated
            global history
            expmd5 = expdat.info['data_md5']
            if expmd5 in history:
                logger.debug('Filling annotation details from history')
                self.fill_from_annotation(history[expmd5], onlyall=True)
        else:
            logger.debug('Filling annotation details from supplied annotation')
            self.fill_from_annotation(prefill_annotation, onlyall=False)

        # set the default primer id for the annotation sequences
        if get_seq_region:
            if self.primerid is None:
                if dbclass is not None:
                    self.primerid = dbclass.db.get_sequences_primer(selected_features)

        # set the primer region button label to the region
        if self.primerid is not None:
            self.bhistory.setText(self.primerid)
        else:
            self.bhistory.setText('Select Region')

        # self.prefillinfo()
        self.bontoinput.setFocus()
        self.show()

    def _load_ontologies(self, dbclass=None):
        '''Load the ontology term files from the local computer (for autocomplete).
        Also check dbbact server if new ontology terms have been added and update the local files

        The data is stored locally in DBAnnotateSave._ontology_dict and DBAnnotateSave._ontology_from_id

        The needed files are created by a script from ontology files and are:
        data/ontology.pickle:
        dict of {name(str): ontologyid(str)}
            name:
                contains the full term/sysnonim name + "(+"ONTOLOGY NAME+"original term + ")". This is the string displayed to the user
            ontologyid:
                contains a unique id for this term that appears in the data/ontologyfromid.pickle file (loaded to DBAnnotateSave._ontology_from_id).

        data/ontologyfromid.pickle:
        dict of {ontologyid(str): term(str)}
            ontologyid:
                contains a unique id for each of the terms (linked from data/ontologies.pickle or DBAnnotateSave._ontology_dict)
            term:
                the dbbact term name

        For example for the term "united states of america" we have in DBAnnotateSave._ontology_dict key "U.S.A. :GAZ(United States of America)" with value GAZ:00002459
        and in DBAnnotateSave._ontology_from_id we have key "GAZ:00002459" with value "United States of America"

        Parameters
        ----------
        dbclass: dbbact.DBBact
            the class where terms are stored, so we can just ask for new dbbact terms
        '''
        _changed = False

        # if ontology files have not been loaded yet, load them
        if len(DBAnnotateSave._ontologies_loaded) == 0:
            _changed = True
            DBAnnotateSave._ontology_from_id = {}
            DBAnnotateSave._ontology_dict = {}
            logger.info('Loading ontology autocomplete')
            ontology_dir = resource_filename(__package__, 'data')
            for cfile in os.listdir(ontology_dir):
                if not cfile.endswith('.ontology.pickle'):
                    continue
                # load the ontology file (key=full name (with ENVO:XXXX), value=dbbact_id)
                logger.debug('Loading ontology %s' % cfile)
                ontology = pickle.load(open(os.path.join(ontology_dir, cfile), 'rb'))
                DBAnnotateSave._ontology_dict.update(ontology)
                # load the ontology ids file (key=dbbact_id, value=name (without ENVO:XXXX))
                logger.debug('Loading ontology ids %s' % cfile)
                ontology_from_id_file = os.path.join(ontology_dir, '.'.join(cfile.split('.')[:-1]) + '.ids.pickle')
                ids_dict = pickle.load(open(ontology_from_id_file, 'rb'))
                DBAnnotateSave._ontology_from_id.update(ids_dict)
                DBAnnotateSave._ontologies_loaded.add(cfile)
                logger.debug('Loaded')
            DBAnnotateSave._ontology_dbbact_max_id = DBAnnotateSave._ontology_dict[max(DBAnnotateSave._ontology_dict, key=DBAnnotateSave._ontology_dict.get)]
            logger.debug('Finished loading ontologies for autocomplete')

        # get the new terms added to dbbact since we got the ontologies files
        if dbclass is not None:
            try:
                logger.debug('current ontology term max_id is %d' % DBAnnotateSave._ontology_dbbact_max_id)
                logger.debug('getting updates to ontology terms')
                new_terms_names, new_terms_ids = dbclass.db.get_ontology_terms(min_term_id=DBAnnotateSave._ontology_dbbact_max_id)
                # we have new terms, let's add them to our autocomplete
                if len(new_terms_names) > 0:
                    logger.debug('found %d new terms' % len(new_terms_names))
                    _changed = True
                    # create the nice text formatting for the new terms
                    # it is of the form:
                    # term (term_ontology_id)
                    # i.e.:
                    # feces (UBERON:000004)
                    new_term_strings = {}
                    new_ontology_from_id = {}
                    for cterm, cid in new_terms_names.items():
                        cterm_str = '%s (%s)' % (cterm, new_terms_ids[str(cid)])
                        new_term_strings[cterm_str] = cid
                        new_ontology_from_id[cid] = cterm
                    DBAnnotateSave._ontology_dict.update(new_term_strings)
                    DBAnnotateSave._ontology_from_id.update(new_ontology_from_id)
                    DBAnnotateSave._ontology_dbbact_max_id = DBAnnotateSave._ontology_dict[max(DBAnnotateSave._ontology_dict, key=DBAnnotateSave._ontology_dict.get)]
                    logger.debug('new ontology dbbact_max_id: %d' % DBAnnotateSave._ontology_dbbact_max_id)
                # if more than 10 new terms, write them to the updated_terms pickle files
                if len(new_terms_names) > 10:
                    logger.info('found %d new ontology terms. updating local file' % len(new_terms_names))
                    # load the existing, update and save for updated_terms.ontology.pickle
                    # it is {autocomplete_term_name(str): id(int)}
                    user_ontology_file = resource_filename(__package__, 'data/updated_terms.ontology.pickle')
                    if os.path.isfile(user_ontology_file):
                        logger.debug('loading user ontology file')
                        contology = pickle.load(open(user_ontology_file, 'rb'))
                        logger.debug('loaded %d terms' % len(contology))
                    else:
                        contology = {}
                    contology.update(new_term_strings)
                    with open(user_ontology_file, 'wb') as ofl:
                            pickle.dump(contology, ofl)
                            logger.debug('Created user ontology file %s' % user_ontology_file)
                    # load the existing, update and save for updated_terms.ontology.pickle
                    # it is {autocomplete_term_name(str): id(int)}
                    user_ontology_file = resource_filename(__package__, 'data/updated_terms.ontology.ids.pickle')
                    if os.path.isfile(user_ontology_file):
                        logger.debug('loading user ontology id file %s' % user_ontology_file)
                        contology_ids = pickle.load(open(user_ontology_file, 'rb'))
                        logger.debug('loaded %d terms' % len(contology))
                    else:
                        contology_ids = {}
                    contology_ids.update(new_ontology_from_id)
                    with open(user_ontology_file, 'wb') as ofl:
                            pickle.dump(contology_ids, ofl)
                            logger.debug('Created user ontology-ids file %s' % user_ontology_file)
            except AttributeError as e:
                logger.warn('error encountered when updating user ontology files: %s' % e)

        if len(DBAnnotateSave._ontology_sorted_list) == 0:
            _changed = True

        # if we changed the terms, re-sort
        if _changed:
            logger.debug('Sorting ontologies')
            DBAnnotateSave._ontology_sorted_list = sorted(list(DBAnnotateSave._ontology_dict.keys()), key=lambda s: s.lower())

    # def history(self):
    #     curtext = []
    #     for cur in hs.lastcurations:
    #         ct = ''
    #         for dat in cur:
    #             ct += dat[0] + '-' + dat[1] + ','
    #         curtext.append(ct)
    #     slistwin = SListWindow(curtext, 'select curation from history')
    #     res = slistwin.exec_()
    #     if res:
    #         items = slistwin.lList.selectedItems()
    #         for citem in items:
    #             print(citem)
    #             spos = slistwin.lList.row(citem)
    #             print(spos)
    #             self.fillfromcuration(hs.lastcurations[spos], onlyall=False)

    def fill_from_annotation(self, annotation, onlyall=True, clearit=True):
        """Fill gui list from annotation list

        Parameters
        ----------
        annotation : dict containing:
            'annotations': list of (annotation type, annotation)
            'description' : str
            'method' : str
            'annotation_type' : str
            'primerid' : str
        onlyall : bool
            True to show only annotations which have ALL, False to show also HIGH/LOW
        clearit : bool
            True to remove previous annotations from list, False to keep
        """
        if clearit:
            self.blistall.clear()
        if 'details' in annotation:
            for cdat in annotation['details']:
                if onlyall:
                    if cdat[0].lower() != 'all':
                        continue
                if len(cdat) >= 3:
                    cterm_id = cdat[2]
                else:
                    cterm_id = None
                self.addtolist(cdat[0], cdat[1], cterm_id)
        if 'description' in annotation:
            self.bdescription.setText(annotation['description'])
        if 'method' in annotation:
            self.bmethod.setText(annotation['method'])
        if 'primer' in annotation:
            self.primerid = annotation['primer']
        # activate the appropriate annotation type buttons
        if 'annotationtype' in annotation:
            atype = annotation['annotationtype'].lower()
            if atype == 'common':
                ctypeidx = self.bisatype.findText('Common', Qt.MatchContains)
                self.bisa.setChecked(True)
                self.bisatype.setCurrentIndex(ctypeidx)
            elif atype == 'dominant':
                ctypeidx = self.bisatype.findText('Dominant', Qt.MatchContains)
                self.bisa.setChecked(True)
                self.bisatype.setCurrentIndex(ctypeidx)
            elif atype == 'other':
                ctypeidx = self.bisatype.findText('Other', Qt.MatchContains)
                self.bisa.setChecked(True)
                self.bisatype.setCurrentIndex(ctypeidx)
            elif atype == 'positive association':
                ctypeidx = self.bisatype.findText('Positively associated', Qt.MatchContains)
                self.bisa.setChecked(True)
                self.bisatype.setCurrentIndex(ctypeidx)
            elif atype == 'negative association':
                ctypeidx = self.bisatype.findText('Negatively associated', Qt.MatchContains)
                self.bisa.setChecked(True)
                self.bisatype.setCurrentIndex(ctypeidx)
            elif atype == 'diffexp':
                self.bdiffpres.setChecked(True)
                pass
            elif atype == 'contamination':
                ctypeidx = self.bisatype.findText('Contam', Qt.MatchContains)
                self.bisa.setChecked(True)
                self.bisatype.setCurrentIndex(ctypeidx)
            self.radiotoggle()

    def radiotoggle(self):
        if self.bisa.isChecked():
            if 'negatively associated' in self.bisatype.currentText().lower():
                self.blow.setDisabled(False)
            else:
                self.blow.setDisabled(True)
            if 'positively associated' in self.bisatype.currentText().lower():
                self.bhigh.setDisabled(False)
            else:
                self.bhigh.setDisabled(True)
        if self.bdiffpres.isChecked():
            self.blow.setEnabled(True)
            self.bhigh.setEnabled(True)

    def isatypechanged(self):
        """
        changed the selection of isatype combobox so need to activate the isa radio button
        """
        self.bisa.setChecked(True)
        self.radiotoggle()

    def select_primers(self, dbclass=None):
        primer_info = dbclass.db.get_primers()
        primers = ['%s (%s:%s - %s)' % (cpi['name'], cpi['fprimer'], cpi['fprimerseq'], cpi['rprimer']) for cpi in primer_info]
        listwin = SListWindow(primers, listname='Select region')
        res = listwin.exec_()
        if res == QtWidgets.QDialog.Accepted:
            # pos = listwin.get_selected()
            pos = listwin.selected_pos
            cprime = primer_info[pos]['name']
            # set the primer region name on the button
            self.bhistory.setText(cprime)
            self.primerid = cprime
        return None

    def studyinfo(self):
        study_data_ui(self.cexp)

    def keyPressEvent(self, e):
        """
        override the enter event so will not close dialog
        """
        e.ignore()

    def minus(self):
        """
        delete selected item from current list
        """
        items = self.blistall.selectedItems()
        for citem in items:
            self.blistall.takeItem(self.blistall.row(citem))

    def cleartext(self):
        self.bontoinput.setText('')

    def plus(self):
        conto = str(self.bontoinput.text())
        cgroup = self.getontogroup()
        self.addtolist(cgroup, conto)
        self.cleartext()

    def ontochanged(self):
        '''called when the ontology text box is changed by the user
        we use it to enable/disable the autocompleter based on the text length in the ontology textbox in order to speed up the response when < 4 chars are entered (too many options)
        '''
        conto = str(self.bontoinput.text())
        if len(conto) > 3:
            self.bontoinput.setCompleter(self._ontology_completer)
        else:
            self.bontoinput.setCompleter(None)

    def addtolist(self, cgroup, conto, conto_term_id=None):
        """
        add an ontology term to the list

        input:
        cgroup : str
            the group (i.e. 'low/high/all')
        conto : str
            the ontology term to add
        conto_term_id: str or None, optional
            the ontologytable term_id for the given ontoloty term (i.e. GAZ:00003).
            If none, deduce the term_id from the conto string (by taking the last part of the string "XXX (YYY GAZ:0004)"
        """
        if conto == '':
            logger.debug('no string to add to list')
            return
        logger.debug('addtolist %s %s (%s)' % (cgroup, conto, conto_term_id))

        # is it a term from an existing ontology?
        if conto in self._ontology_dict:
            # this is a hack to get the ontology term_id (i.e. "GAZ:00001" from the description string)
            # we could do it more ordered but will take more memory - so let's hack in the meanwhile
            if conto_term_id is None:
                conto_term_id = conto.split('(')[-1].split(' - ')[-1].split(')')[0]
                if ':' not in conto_term_id:
                    logger.warn('No ontology term_id in current string %s' % conto)
            conto = self._ontology_from_id[self._ontology_dict[conto]]
        else:
            if conto_term_id is None:
                logger.debug('Term %s not in ontology.' % conto)
                res = show_and_ask('Term %s does not appear in any ontology.\n'
                                   'Do you want to add it as a new term to dbBact?' % conto,
                                   keyval='new_ontology_term', show_cancel=True)
                if not res:
                    logger.debug('cancel all term to list')
                    return
                conto_term_id = ''

        # TODO: check if there are more than 1 ontology terms with the same description (i.e. 'feces'), let the user select which ones to add

        # check if item already in list
        for citem in qtlistiteritems(self.blistall):
            cdata = qtlistgetdata(citem)
            if cdata['value'] == conto:
                # if already in list and it is a new term, do nothing
                if conto_term_id == '':
                    logger.debug('item already in list')
                    return
                # if specific term_id already in the list of term_ids, do nothing
                if conto_term_id == cdata['term_id']:
                    logger.debug('term_id %s for term %s already in list' % (conto_term_id, conto))
                    return

        if cgroup.lower() == 'low':
            ctext = "LOW:%s (%s)" % (conto, conto_term_id)
            qtlistadd(self.blistall, ctext, {'type': 'LOW', 'value': conto, 'term_id': conto_term_id}, color='red')
        if cgroup.lower() == 'high':
            ctext = "HIGH:%s (%s)" % (conto, conto_term_id)
            qtlistadd(self.blistall, ctext, {'type': 'HIGH', 'value': conto, 'term_id': conto_term_id}, color='blue')
        if cgroup.lower() == 'all':
            ctext = "ALL:%s (%s)" % (conto, conto_term_id)
            qtlistadd(self.blistall, ctext, {'type': 'ALL', 'value': conto, 'term_id': conto_term_id}, color='black')

    def getontogroup(self):
        if self.ball.isChecked():
            return('ALL')
        if self.blow.isChecked():
            return('LOW')
        if self.bhigh.isChecked():
            return('HIGH')

    def accept(self):
        err = self.inputCheck()
        if err == '':
            self.done(1)
        else:
            logger.warn(err)

    def inputCheck(self):
        '''Validate we have enough information in the annotation
        We link here from the accept slot

        Returns
        -------
        str
        '' (empty string) if details valid, non-empty error string if problem encountered
        '''
        # check we have details
        if len(list(qtlistiteritems(self.blistall))) == 0:
            QtWidgets.QMessageBox.warning(self, 'Missing information', 'No entries in annotation')
            return ('No entries in annotation')

        # if small amount of details, verify
        if len(list(qtlistiteritems(self.blistall))) < 3:
            if QtWidgets.QMessageBox.warning(self, 'Submit annotation?', 'Less than three entries for the annotation. Do you want to submit?',
                                             QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.No:
                return 'Fill more entries'

        # get the annotations
        types = set()
        for citem in qtlistiteritems(self.blistall):
            cdat = qtlistgetdata(citem)
            ctype = cdat['type'].lower()
            types.add(ctype)

        # check we chose a primer region
        if self.primerid is None:
            msg = 'Please select primer region for annotation'
            return msg

        # if differential expression, check there is high and low
        if self.bdiffpres.isChecked():
            msg = None
            if 'high' not in types:
                msg = 'Missing "high" entries for differential abundance'
            if 'low' not in types:
                msg = 'Missing "low" entries value for differential abundance'
            if msg is not None:
                QtWidgets.QMessageBox.warning(self, 'Missing information', msg)
                return msg
        return ''

    def prefillinfo(self):
        """
        prefill "ALL" data fields based on mapping file
        if all samples have same info
        """
        logger.debug('prefill info')
#       ontologyfromid=self.ontologyfromid
# #     fl=open('/Users/amnon/Python/git/heatsequer/db/ncbitaxontofromid.pickle','rb')
#       fl=open(os.path.join(hs.heatsequerdir,'db/ncbitaxontofromid.pickle'),'rb')
#       ncbitax=pickle.load(fl)
#       fl.close()

#       cexp=self.cexp
#       for cfield in cexp.fields:
#           uvals=[]
#           if cfield in cexp.fields:
#               uvals=hs.getfieldvals(cexp,cfield,ounique=True)
#           # if we have 1 value
#           if len(uvals)==1:
#               cval=uvals[0]
#               hs.Debug(1,'found 1 value %s' % cval)
#               if cfield=='HOST_TAXID' or cfield=='host_taxid':
#                   hs.Debug(2,'%s field has 1 value %s' % (cfield,cval))
#                   # if ncbi taxonomy (field used differently)
#                   cval='NCBITaxon:'+cval
#                   if cval in ncbitax:
#                       hs.Debug(2,'found in ncbitax %s' % cval)
#                       cval=ncbitax[cval]
#               else:
#                   # get the XXX from ENVO:XXX value
#                   uvalspl=cval.split(':',1)
#                   if len(uvalspl)>1:
#                       cval=uvalspl[1]
#                       cval=uvalspl[1]+' :'+uvalspl[0]
#               if cval in self.ontology:
#                   cval=ontologyfromid[self.ontology[cval]]
#                   hs.Debug(2,'term %s found in ontologyfromid' % cval)
#                   conto=cval
#                   hs.Debug(1,'add prefill %s' % conto)
#                   self.addtolist('ALL',conto)
#               else:
#                   hs.Debug(3,'term %s NOT found in ontologyfromid' % uvals[0])

#           else:
#               hs.Debug(1,'found %d values' % len(uvals))


def getqtlistitems(qtlist):
    """
    get a list of strings of the qtlist
    input:
    qtlist : QTListWidget

    output:
    item : list of str
    """
    items = []
    for index in range(qtlist.count()):
        items.append(str(qtlist.item(index).text()))
    return items


def qtlistadd(qtlist, text, data, color="black"):
    """
    Add an entry (text) to qtlist and associaxte metadata data
    input:
    qtlist : QTListWidget
    text : str
        string to add to list
    data : arbitrary python var
        the data to associate with the item (get it by qtlistgetdata)
    color : (R,G,B)
        the color of the text in the list
    """
    item = QtWidgets.QListWidgetItem()
    item.setText(text)
    ccol = QtGui.QColor()
    ccol.setNamedColor(color)
    item.setForeground(ccol)
    item.setData(Qt.UserRole, data)
    qtlist.addItem(item)


def qtlistgetdata(item):
    """
    Get the metadata associated with item as position pos
    input:
    qtlist : QtListWidget
    index : QtListWidgetItem
        the item to get the info about

    output:
    data : arbitrary
        the data associated with the item (using qtlistadd)
    """
    data = item.data(Qt.UserRole)
    return data


def qtlistiteritems(qtlist):
    """
    iterate all items in a list
    input:
    qtlist : QtListWidget
    """
    for i in range(qtlist.count()):
        yield qtlist.item(i)


class UserPasswordDlg(QtWidgets.QDialog):
    '''Show/get details about a study
    '''
    def __init__(self):
        '''Init the values in the study info dialog

        '''
        super().__init__()
        uic.loadUi(resource_filename(__package__, 'ui/user_password.ui'), self)
        self.login.clicked.connect(self.login_func)
        self.anonymous.clicked.connect(self.anonymous_func)

    def login_func(self):
        if self.email.text() == '':
            if QtWidgets.QMessageBox.warning(self, 'email missing', 'No email supplied, you will not be able to recover your password.\nContinue?',
                                             QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.No:
                return
        self.accept()

    def anonymous_func(self):
        self.close()


def show_and_ask(msg, keyval, show_cancel=False):
    '''Show the message unless the user has checked the "Don't show this again" box

    Parameters
    ----------
    msg: str
        the message to display
    keyval: str
        the config file key to indicate user has asked to not show again the message
    show_cancel: bool, optional
        if True, also show a cancel button.
        If False (default), only have an ok button

    Returns
    -------
    bool:
        True if ok was pressed, False if Cancel was pressed (only if show_cancel is True)
    '''
    keyval = 'skip_msg_' + keyval
    res = get_config_value(keyval, section='dbbact', fallback='no')
    if res.lower() == 'yes':
        return
    a = QtWidgets.QMessageBox()
    a.setText(msg)
    a.setWindowTitle('dbBact-Calour')
    a.setIcon(QtWidgets.QMessageBox.Information)
    if not show_cancel:
        a.setStandardButtons(QtWidgets.QMessageBox.Ok)
    else:
        a.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)
    cb = QtWidgets.QCheckBox(text="Don't show this again")
    a.setCheckBox(cb)
    res = a.exec_()
    if cb.isChecked():
        set_config_value(keyval, 'yes', section='dbbact')
    if res == QtWidgets.QMessageBox.Ok:
        return True
    return False


class SListWindow(QtWidgets.QDialog):
    def __init__(self, listdata=[], listname=None):
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

        self.w_list = QListWidget()
        self.layout.addWidget(self.w_list)
        self.selected_pos = None

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept_get_pos)
        # buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        self.layout.addWidget(buttonBox)

        for citem in listdata:
            self.w_list.addItem(citem)

        self.w_list.itemDoubleClicked.connect(self.list_double_click)

        self.show()
        self.adjustSize()

    def accept_get_pos(self):
        self.selected_pos = self.get_selected()
        self.accept()

    def add_item(self, text, color='black', dblclick_data=None):
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
        self.w_list.addItem(item)

    def list_double_click(self, item):
        data = item.data(QtCore.Qt.UserRole)
        if data is not None:
            data['database'].show_term_details(data['term'], data['exp'], data['features1'], data['features2'], gui='qt5')

    def get_selected(self):
        items = self.w_list.selectedItems()
        for citem in items:
            spos = self.w_list.row(citem)
        return spos
