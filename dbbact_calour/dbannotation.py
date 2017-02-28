from logging import getLogger
from pkg_resources import resource_filename
import pickle
import sys

from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication
import pandas as pd
from calour.util import get_config_value, set_config_value

from . import dbbact


logger = getLogger(__name__)


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


def annotate_bacteria_gui(db, seqs, exp):
    '''Create a dialog for annotating the bacteria into dbbact

    Parameters
    ----------
    db : dbbact.DBBact
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
    test_user_password(db)

    # test if study already in database
    cdata = db.find_experiment_id(datamd5=exp.exp_metadata['data_md5'], mapmd5=exp.exp_metadata['map_md5'])
    if cdata is None:
        logger.info('Study does not exist in dbbact. Creating new study')
        cdata = study_data_ui(exp)
        if cdata is None:
            msg = 'no study information added. Please add at least one field. Annotation cancelled'
            logger.warn(msg)
            return msg

    ui = DBAnnotateSave(seqs, exp)
    res = ui.exec_()
    if res == QtWidgets.QDialog.Accepted:
        description = str(ui.bdescription.text())
        # TODO: need to get primer region!!!!
        primerid = 'V4'
        method = str(ui.bmethod.text())
        if method == '':
            method = 'na'
        # TODO: fix submitter name
        submittername = 'Amnon Amir'
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
            elif 'High' in curtypeval:
                annotation_type = 'HIGHFREQ'
            else:
                annotation_type = 'OTHER'
        else:
            msg = "No annotation type selected. Annotation not saved"
            logger.warn(msg)
            return msg

        logger.debug('Adding annotation to studyid %s' % cdata)
        res = db.add_annotations(expid=cdata, sequences=seqs, annotationtype=annotation_type, annotations=annotations, submittername=submittername, description=description, method=method, primerid=primerid)
        if res is None:
            msg = 'Annotation not added.'
            logger.warn(msg)
            return msg
        logger.debug('New annotation added. AnnotationId=%d' % res)
        # # store the history
        # try:
        #     hs.lastcurations.append(curations)
        # except:
        #     hs.lastcurations=[curations]
        # hs.lastdatamd5=self.cexp.datamd5
        return ''
    return 'Add annotation cancelled'


def test_user_password(db):
    '''
    Test if the config file has user/password and if not ask for one.
    Also don't ask if the 'ask again' checkbox is checked

    Parameters
    ----------
    db : DBBact
        the database interface class
    '''
    username = get_config_value('username', section='dbBact')
    if username is not None:
        return
    if get_config_value('show_user_request', section='dbBact') is not None:
        logger.debug('user/password not set, but show_user_request flag in config file is set, so ignoring')
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
                return

        logger.info('storing username %s in config file' % username)
        set_config_value('username', username, section='dbBact')
        set_config_value('password', password, section='dbBact')


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
    cid = bdb.find_experiment_id(datamd5=cexp.exp_metadata['data_md5'], mapmd5=cexp.exp_metadata['map_md5'])

    # add study details we have
    study_details = []
    if cid is None:
        study_details.append(('DataMD5', cexp.exp_metadata['data_md5']))
        study_details.append(('MapMD5', cexp.exp_metadata['map_md5']))
        interesting_columns = {'sra_study_s': 'sra', 'project_name_s': 'name', 'experiment_title': 'name', 'experiment_design_description': 'name', 'BioProject_s': 'sra', 'run_date': 'date'}
        for ccol in cexp.sample_metadata.columns:
            if ccol.lower() in interesting_columns:
                if len(pd.unique(cexp.sample_metadata[ccol].ravel())) == 1:
                    study_details.append((interesting_columns[ccol.lower()], cexp.sample_metadata[ccol][0]))

    dbsi = DBStudyInfo(cid, study_details=study_details)
    res = dbsi.exec_()
    if res == QtWidgets.QDialog.Accepted:
        newstudydata = []
        allstudydata = []
        for citem in qtlistiteritems(dbsi.blist):
            cdata = qtlistgetdata(citem)
            allstudydata.append((cdata['type'], cdata['value']))
            if cdata['fromdb'] == False:
                newstudydata.append((cdata['type'], cdata['value']))

        if len(allstudydata) <= 2:
            logger.warn('not enough details. not saving anything')
            return None
        if len(newstudydata) == 0:
            logger.warn('No new items. not saving anything')
            return cid
        dataid = bdb.add_experiment_data(newstudydata, expid=cid)
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
            usually init to the data_d5 and map_md5
        '''
        super(DBStudyInfo, self).__init__()
        uic.loadUi(resource_filename(__package__, 'ui/studyinfo.ui'), self)
        self.dataid = 0
        bdb = dbbact.DBBact()

        # add experiment details from the database
        if experimentid is not None:
            info = bdb.get_experiment_info(experimentid)
            for cinfo in info:
                qtlistadd(self.blist, cinfo[0]+':'+cinfo[1], {'fromdb': True, 'type': cinfo[0], 'value': cinfo[1]}, color='grey')

        # and the supplied details (supplied to the init function via exp_details)
        for cinfo in study_details:
            qtlistadd(self.blist, cinfo[0] + ':' + cinfo[1], {'fromdb': False, 'type': cinfo[0], 'value': cinfo[1]}, color='black')

        self.bplus.clicked.connect(self.plus)
        self.bvalue.returnPressed.connect(self.plus)
        self.bminus.clicked.connect(self.minus)
        self.bannotations.clicked.connect(self.annotations)
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
            if cdata['fromdb'] == False:
                newstudydata.append((cdata['type'], cdata['value']))
        # validate we have something else than the md5 for data/map
        if len(allstudydata) <= 2:
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
    # used to store the ontology values for the autocomplete
    _ontology_dict = None

    def __init__(self, selected_features, expdat):
        '''Create the manual annotation window

        Parameters
        ----------
        selected_features : list of sequences ('ACGT')
            The sequences to annotate
        expdat : Experiment
            The experiment being annotated
        '''
        super(DBAnnotateSave, self).__init__()

        # create the gui
        uic.loadUi(resource_filename(__package__, 'ui/manualdata.ui'), self)
        self.bplus.clicked.connect(self.plus)
        self.bminus.clicked.connect(self.minus)
        self.bontoinput.returnPressed.connect(self.plus)
        self.bstudyinfo.clicked.connect(self.studyinfo)
        self.bisa.toggled.connect(self.radiotoggle)
        self.bdiffpres.toggled.connect(self.radiotoggle)
        self.bisatype.currentIndexChanged.connect(self.isatypechanged)
        self.bhistory.clicked.connect(self.history)
        self.cexp = expdat
        self.lnumbact.setText(str(len(selected_features)))
        completer = QtWidgets.QCompleter()
        self.bontoinput.setCompleter(completer)

        # create the autocomplete ontology list
        model = QtCore.QStringListModel()
        completer.setModel(model)
        # completer.setCompletionMode(QCompleter.InlineCompletion)
        completer.maxVisibleItems = 10
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        # make the completer selection also erase the text edit
        completer.activated.connect(self.cleartext, type=Qt.QueuedConnection)
        # init the ontology values
        self._load_ontologies()
        model.setStringList(self._ontology_sorted_list)

        self.setWindowTitle(expdat.description)
        # try:
        #   tt=hs.lastdatamd5
        # except:
        #   hs.lastdatamd5=''
        # if self.cexp.datamd5==hs.lastdatamd5:
        #   self.fillfromcuration(hs.lastcurations[-1],onlyall=True)

        self.prefillinfo()
        self.bontoinput.setFocus()
        self.show()

    def _load_ontologies(self):
        if DBAnnotateSave._ontology_dict is None:
            print('loading')
            ontology_file = resource_filename(__package__, 'data/ontology.pickle')
            ontology = pickle.load(open(ontology_file, 'rb'))
            DBAnnotateSave._ontology_dict = ontology
            print('sorting')
            olist = list(DBAnnotateSave._ontology_dict.keys())
            DBAnnotateSave._ontology_sorted_list = sorted(olist, key=lambda s: s.lower())

            ontology_from_id_file = resource_filename(__package__, 'data/ontologyfromid.pickle')
            DBAnnotateSave._ontology_from_id = pickle.load(open(ontology_from_id_file, 'rb'))

    def history(self):
        curtext = []
        for cur in hs.lastcurations:
            ct = ''
            for dat in cur:
                ct += dat[0] + '-' + dat[1] + ','
            curtext.append(ct)
        slistwin = SListWindow(curtext, 'select curation from history')
        res = slistwin.exec_()
        if res:
            items = slistwin.lList.selectedItems()
            for citem in items:
                print(citem)
                spos = slistwin.lList.row(citem)
                print(spos)
                self.fillfromcuration(hs.lastcurations[spos], onlyall=False)

    def fillfromcuration(self, curation, onlyall=True, clearit=True):
        """
        fill gui list from curation
        input:
        curation : from hs.lastcurations
        onlyall : bool
            True to show only curations which have ALL, False to show also HIGH/LOW
        clearit : bool
            True to remove previous curations from list, False to keep
        """
        if clearit:
            self.blistall.clear()
        for cdat in curation:
            if onlyall:
                if cdat[0] != 'ALL':
                    continue
            self.addtolist(cdat[0], cdat[1])

    def radiotoggle(self):
        if self.bisa.isChecked():
            self.blow.setDisabled(True)
            self.bhigh.setDisabled(True)
        if self.bdiffpres.isChecked():
            self.blow.setEnabled(True)
            self.bhigh.setEnabled(True)

    def isatypechanged(self):
        """
        changed the selection of isatype combobox so need to activate the isa radio button
        """
        self.bisa.setChecked(True)

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

    def addtolist(self, cgroup, conto):
        """
        add an ontology term to the list

        input:
        cgroup : str
            the group (i.e. 'low/high/all')
        conto : str
            the ontology term to add
        """
        if conto == '':
            logger.debug('no string to add to list')
            return
        print('addtolist %s %s' % (cgroup, conto))
        if conto in self._ontology_dict:
            conto = self._ontology_from_id[self._ontology_dict[conto]]
        else:
            logger.debug('Not in ontology!!!')
            # TODO: add are you sure... not in ontology list....

        # if item already in list, don't do anything
        for citem in qtlistiteritems(self.blistall):
            cdata = qtlistgetdata(citem)
            if cdata['value'] == conto:
                logger.debug('item already in list')
                return

        if cgroup == 'LOW':
            ctext = "LOW:%s" % conto
            qtlistadd(self.blistall, ctext, {'type': 'LOW', 'value': conto}, color='red')
        if cgroup == 'HIGH':
            ctext = "HIGH:%s" % conto
            qtlistadd(self.blistall, ctext, {'type': 'HIGH', 'value': conto}, color='blue')
        if cgroup == 'ALL':
            ctext = "ALL:%s" % conto
            qtlistadd(self.blistall, ctext, {'type': 'ALL', 'value': conto}, color='black')

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

        # if differential expression, check there is high and low
        if self.bdiffpres.isChecked():
            if 'high' not in types:
                return 'Missing "high" entries for differential abundance'
            if 'low' not in types:
                return 'Missing "low" entries value for differential abundance'
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
