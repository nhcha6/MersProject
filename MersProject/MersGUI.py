# pyinstaller MersGUI --> this command from the relevant file location creates executable file
# Need to fix no parameter
import sys
import subprocess
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit, QInputDialog, QGroupBox, QFormLayout
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
import matplotlib as mpl
import queue
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from Mers import *
from MGFMain import *
import functools
from functools import partial
from datetime import datetime
from pathlib import Path
import signal
import platform


class WorkerSignals(QObject):
    """
    This class declares a series of custom signals that QRunnable objects can emit to call multiple functions
    from within a thread.
    """

    finished = pyqtSignal()
    updateProgBar = pyqtSignal()
    disableButtons = pyqtSignal()
    plot = pyqtSignal(list, list)


class ProgressGenerator(QRunnable):
    """
    A QRunnable class that enables a progress bar to be created and updated while the output is being run. This must
    be done in a separate thread to the main GUI functionality to ensure the rest of the interface doesn't freeze.
    Used to create an object in self.outputPreStep().
    """

    def __init__(self, *args, **kwargs):
        super(ProgressGenerator, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.flag = True

    def changeFlag(self):
        """
        Called from within the MyTableWidget.outputFinished to change self.flag to False and quit the thread.
        :return:
        """
        self.flag = False

    @pyqtSlot()
    def run(self):
        """
        This function is run when a ProgressGenerator QRunnable is added to the threadpool. It simply updates the
        progress bar every 10ms via the updateProgBar signal until self.flag is changed to False.
        :return:
        """
        self.signals.disableButtons.emit()
        while self.flag:
            self.signals.updateProgBar.emit()
            time.sleep(0.01)
        self.signals.finished.emit()


class OutputGenerator(QRunnable):
    """
    A QRunnable class which enables self.output() to run in a different thread to the main GUI and the progress bar.
    This ensures that the GUI doesn't freeze while the spliced peptides are being created. This class is used to create
    an object in self.outputPreStep().
    """

    def __init__(self, fn, *args, **kwargs):
        super(OutputGenerator, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        """
        This function runs when an OutputGenerator class is added to the threadpool. It simply runs the worker function
        self.output() and then emits the finished signal which runs self.outputFinished().
        :return:
        """
        self.fn(*self.args)
        self.signals.finished.emit()


class MGFImporter(QRunnable):
    """
    A QRunnable class which enables the MGF file data to be extracted in an alternative thread to the main interface.
    This class is used to create an object in self.returnPath().
    """

    def __init__(self, fn, *args, **kwargs):
        super(MGFImporter, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        """
        Runs the worker function self.importMGF() before emitting the the finished signal, which run self.importedMGf().
        :return:
        """
        start = time.time()
        self.fn(*self.args)
        self.signals.finished.emit()
        end = time.time()
        print("Uploading mgf took: " + str(end - start))


class MGFPlotter(QRunnable):
    """
    A QRunnable class which is used to create an alternative thread to plot the maximum intensities of each spectra
    in the input MGF file. This class is used to create an object in self.uploadMgfPreStep().
    """

    def __init__(self, fn, *args, **kwargs):
        super(MGFPlotter, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        """
        This function runs when an MGFPlotter object is added to the threadpool. It runs the worker function plotData()
        from MGFMain.py which returns the x and y data to be plotted. It then emits the plot signal which is connected
        to self.onlyImportMGF() and ultimately creats the plot. It then emits the finished signal which is connected to
        self.intensityPlotFin() cleanly break the thread and enable all widgets.
        :return:
        """
        start = time.time()
        ms2Thresh, intensityPoints = self.fn(*self.args)

        self.signals.plot.emit(ms2Thresh, intensityPoints)
        end = time.time()
        print("Getting plot data took: " + str(end - start))
        self.signals.finished.emit()


class App(QMainWindow):
    """
    App serves as the parent class for the embedded MyTableWidget.
    Methods: center(), closeEvent()
    Class Variables:
        self.title: the title of the App, present window no matter which tab the user navigates to.
        self.statusBar: a status bar that can be updated to allow the offer information to the user.
        self.table_widget: the widget inbedded within the main App window. This widget stores most of the detail
        of the GUI.
    """

    def __init__(self):
        """
         Initialisation of main window class

        """

        super().__init__()
        self.title = 'Peptide Splicer'
        self.statusbar = self.statusBar()

        self.center()

        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)
        self.setWindowTitle(self.title)

        self.show()

    def center(self):
        """
        Called by self.__init__(), this function centres the App window in the middle of the users screen.
        :return:

        """

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def closeEvent(self, event):
        """
        Called automatically on click of the exit button. This function is configured to test if the system is being
        run on windows or mac, and uses the appropriate command to close the window and kill all processes which were
        generated by the program.

        :param event: details that close event is taking place.
        :return:
        """
        print('closed')
        # windows close command
        if platform.system() == 'Windows':
            os.system('taskkill /f /fi "WINDOWTITLE eq Peptide Splicer" /t')
        # mac close command
        else:
            os.system("ps aux |grep MersGUI | grep -v 'pattern_of_process_you_dont_want_to_kill' | awk '{print $2}' |xargs kill")


class MyTableWidget(QWidget):
    """
    MyTableWidget class is the child of the App class. It holds the tabs where most of the GUI functionality occurs.
    """

    def __init__(self, parent):

        """
        Initialisation of table child
        """

        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        # self.fasta will store the fasta class once initialised
        self.fasta = None
        self.mgf = None
        self.mgfPath = None
        self.outputPath = None
        self.settingString = None

        # Init threading
        self.threadpool = QThreadPool()

        # Default values for the input parameters
        self.minDefault = '8'
        self.maxDefault = '12'
        self.maxDistDefault = '25'
        # self.minByIonDefault = '50'
        # self.byIonAccDefault = '0.1'

        # Initialisation of two tabs
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tabs.resize(300, 200)

        # Add tabs to table class (self)
        self.tabs.addTab(self.tab1, "Select File and Path")
        self.tabs.addTab(self.tab2, "Input Parameters")
        # self.tabs.setTabEnabled(1, False)

        # Creation of tab layout and widgets within tab
        self.tab1.layout = QGridLayout(self)
        self.tab1.layout.setSpacing(10)
        self.createTab1ParameterWidgets()
        self.addTab1ParameterWidgets()

        self.tab1.setLayout(self.tab1.layout)

        # counter for ensuring additional mods can be added
        self.addModsFlag = False

        # Create second tab layout and widgets within each tab
        self.tab2.layout = QGridLayout(self)
        self.tab2.layout.setSpacing(10)

        self.createTab2ParameterWidgets()
        self.addTab2ParameterWidgets()

        # set layout
        self.tab2.setLayout(self.tab2.layout)

        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    # First set of functions are those which are involved with creating and editing the widgets on each tab.

    def createTab1ParameterWidgets(self):
        """
        Called when MyTableWidget is initialised. Creates all the wigets which are to be added to the first tab of the GUI.
        :return:
        """
        self.pushButton1 = QPushButton("Select Fasta File")
        self.pushButton1.clicked.connect(self.uploadFasta)
        self.addMultipleFasta = QPushButton("Add Another Fasta")
        self.addMultipleFasta.clicked.connect(self.uploadFasta)
        self.addMultipleFasta.setEnabled(False)
        self.mgfButton = QPushButton("Select MGF File")
        self.mgfButton.clicked.connect(self.uploadMgfPreStep)
        self.mgfPlotFlag = QCheckBox('Produce Intensity Plot')
        self.nextTab = QPushButton("Next Tab")
        self.nextTab.clicked.connect(self.nextTabFunc)
        self.nextTab.setEnabled(False)
        self.tabs.setTabEnabled(1, False)
        self.mgfFlag = QCheckBox("No MGF Comparison")
        self.mgfFlag.stateChanged.connect(self.controlMGFInput)

        self.tab1.ppmLabel = QLabel('PPM (0.1 - 1000): ')
        self.tab1.ppmText = QLineEdit(self)
        self.tab1.ppmText.setEnabled(False)
        self.tab1.ppmLabel.setEnabled(False)
        self.tab1.ppmText.textChanged[str].connect(self.textBoxChanged)
        self.tab1.ppmText.textChanged[str].connect(self.enableControl)
        self.tab1.ppmStatus = QLabel("")

        self.tab1.toleranceLabel = QLabel('Intensity Threshold: ')
        self.tab1.toleranceText = QLineEdit(self)
        self.tab1.toleranceText.setEnabled(False)
        self.tab1.toleranceLabel.setEnabled(False)
        self.tab1.toleranceText.textChanged[str].connect(self.textBoxChanged)
        self.tab1.toleranceText.textChanged[str].connect(self.enableControl)
        self.tab1.toleranceStatus = QLabel("")

        self.tab1.minByIonLabel = QLabel('Minimum b/y Ion Matches(%): ')
        self.tab1.minByIonText = QLineEdit(self)
        self.tab1.minByIonStatus = QLabel("")
        self.tab1.minByIonText.setEnabled(False)
        self.tab1.minByIonLabel.setEnabled(False)
        self.tab1.minByIonText.textChanged[str].connect(self.textBoxChanged)
        self.tab1.minByIonText.textChanged[str].connect(self.enableControl)

        self.tab1.byIonAccLabel = QLabel('b/y Ion Accuracy (0.01 - 0.2): ')
        self.tab1.byIonAccText = QLineEdit(self)
        self.tab1.byIonAccText.setEnabled(False)
        self.tab1.byIonAccLabel.setEnabled(False)
        self.tab1.byIonAccText.textChanged[str].connect(self.textBoxChanged)
        self.tab1.byIonAccText.textChanged[str].connect(self.enableControl)
        self.tab1.byIonAccStatus = QLabel("")

        self.tab1.byIonFlag = QCheckBox('Apply b/y Ion Comparison: ')
        self.tab1.byIonFlag.setEnabled(False)
        # self.tab1.byIonFlag.stateChanged.connect(self.disableByInputs)
        self.tab1.byIonFlag.stateChanged.connect(self.enableControl)

    def addTab1ParameterWidgets(self):
        """
        Called when MyTableWidget is initialised. Adds all the widgets that have been created for tab 1 to the grid layout of tab 1.
        :return:
        """
        self.tab1.layout.setColumnStretch(0, 1)
        self.tab1.layout.setColumnStretch(5, 1)
        self.tab1.layout.setRowStretch(0, 1)
        self.tab1.layout.setRowStretch(5, 1)
        self.tab1.layout.addWidget(self.pushButton1, 1, 2)
        self.tab1.layout.addWidget(self.addMultipleFasta, 3, 2)
        self.tab1.layout.addWidget(self.mgfButton, 2, 2)
        self.tab1.layout.addWidget(self.mgfFlag, 1, 3)
        self.tab1.layout.addWidget(self.mgfPlotFlag, 2, 3)
        self.tab1.layout.addWidget(self.tab1.ppmLabel, 4, 2)
        self.tab1.layout.addWidget(self.tab1.ppmText, 4, 3)
        self.tab1.layout.addWidget(self.tab1.ppmStatus, 4, 4)
        self.tab1.layout.addWidget(self.tab1.toleranceLabel, 5, 2)
        self.tab1.layout.addWidget(self.tab1.toleranceText, 5, 3)
        self.tab1.layout.addWidget(self.tab1.toleranceStatus, 5, 4)

        self.tab1.layout.addWidget(self.tab1.minByIonLabel, 6, 2)
        self.tab1.layout.addWidget(self.tab1.minByIonText, 6, 3)
        self.tab1.layout.addWidget(self.tab1.minByIonStatus, 6, 4)
        self.tab1.layout.addWidget(self.tab1.byIonAccLabel, 7, 2)
        self.tab1.layout.addWidget(self.tab1.byIonAccText, 7, 3)
        self.tab1.layout.addWidget(self.tab1.byIonAccStatus, 7, 4)
        self.tab1.layout.addWidget(self.tab1.byIonFlag, 8, 2)
        self.tab1.layout.addWidget(self.nextTab, 8, 3)

    def createTab2ParameterWidgets(self):
        """
        Called when MyTableWidget is initialised. Calls a series of functions which create all the widgets which are to be added to the second tab.
        :return:
        """
        self.addMinMaxAndDist()
        self.addModifications()
        self.addFlagChecks()
        self.addChargeStates()

        # AN EXTRA ADD "WRITE TO CSV FUNCTION CHECKBOX"
        self.tab2.csv = QCheckBox('Write To Csv')
        self.tab2.pepToProt = QCheckBox('Pep to Prot.csv')
        self.tab2.protToPep = QCheckBox('Prot to Pep.csv')

        # create generate output push button
        self.tab2.output = QPushButton('Generate Output!', self)
        #self.tab2.stop = QPushButton('Stop Process', self)
        self.tab2.output.clicked.connect(self.confirmationFunction)
        #self.tab2.stop.clicked.connect(self.stopFunction)
        #self.tab2.stop.setEnabled(False)

        self.setDefaultParameters()

    def addTab2ParameterWidgets(self):

        """
        Called when MyTableWidget is initialised. This function adds all the widgets which have been created for the
        second tab to the grid layout corresponding to this tab.
        """

        # All the labels added to grid layout of tab2
        self.tab2.layout.addWidget(self.tab2.minimum, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximum, 2, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistance, 3, 3)
        self.tab2.layout.addWidget(self.tab2.mod1, 4, 3)
        self.tab2.layout.addWidget(self.tab2.mod2, 5, 3)
        self.tab2.layout.addWidget(self.tab2.mod3, 6, 3)
        if self.addModsFlag:
            self.tab2.layout.addWidget(self.tab2.removeMods, 10, 3)
        else:
            self.tab2.layout.addWidget(self.tab2.addMods, 10, 3)
        self.tab2.layout.addWidget(self.tab2.maxMod, 11, 3)
        self.tab2.layout.addWidget(self.tab2.linear, 12, 3)
        self.tab2.layout.addWidget(self.tab2.cis, 13, 3)
        self.tab2.layout.addWidget(self.tab2.trans, 14, 3)
        self.tab2.layout.addWidget(self.tab2.overlap, 15, 3)
        self.tab2.layout.addWidget(self.tab2.concat, 15, 4)

        self.tab2.layout.addWidget(self.tab2.csv, 12, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.pepToProt, 13, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.protToPep, 14, 4, 1, 3)

        self.tab2.layout.addWidget(self.tab2.chargeLabel, 16, 3)

        # all dynamic elements added to the grid layout of tab 2
        self.tab2.layout.addWidget(self.tab2.minimumCombo, 1, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximumCombo, 2, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistCombo, 3, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod1Combo, 4, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod2Combo, 5, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod3Combo, 6, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maxModCombo, 11, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.plusOne, 17, 4)
        self.tab2.layout.addWidget(self.tab2.plusTwo, 17, 5)
        self.tab2.layout.addWidget(self.tab2.plusThree, 17, 6)
        self.tab2.layout.addWidget(self.tab2.plusFour, 18, 4)
        self.tab2.layout.addWidget(self.tab2.plusFive, 18, 5)
        self.tab2.layout.addWidget(self.tab2.output, 19, 5, 1, 2)
        #self.tab2.layout.addWidget(self.tab2.stop, 15, 3)

    def addMinMaxAndDist(self):
        """
        Called by self.createTab2ParameterWidgets() to create the min, max and max distance combo boxes.
        :return:
        """
        # Minimum/maximum Combo boxes and connector functions
        self.tab2.minimum = QLabel('Minimum Peptide Length : ')
        self.tab2.minimumCombo = QComboBox(self)
        self.tab2.minimumCombo.activated[str].connect(self.minMaxChanged)
        self.tab2.maximum = QLabel('Maximum Peptide Length : ')
        self.tab2.maximumCombo = QComboBox(self)
        self.tab2.maximumCombo.activated[str].connect(self.minMaxChanged)

        # Max distance combo box
        self.tab2.maxDistance = QLabel('Maximum Distance : ')
        self.tab2.maxDistCombo = QComboBox(self)

        # Adding values to the max/min/maxDist combos
        self.tab2.maxDistCombo.addItem('None')
        for i in range(int(self.minDefault), 26):
            self.tab2.maximumCombo.addItem(str(i))
        for i in range(2, int(self.maxDefault) + 1):
            self.tab2.minimumCombo.addItem(str(i))
        for i in range(2, 26):
            self.tab2.maxDistCombo.addItem(str(i))

    def addFlagChecks(self):
        """
        Called by self.createTab2ParameterWidgets to create the concatenation, overlep, cis, linear and trans
        checkboxes.
        :return:
        """
        # initialise overlap, trans, cis and linear check boxes
        self.tab2.overlap = QCheckBox('Cis Overlap Off', self)
        self.tab2.concat = QCheckBox('Concat Output', self)
        self.tab2.trans = QCheckBox('Trans', self)
        self.tab2.trans.stateChanged.connect(self.enableControl)
        self.tab2.cis = QCheckBox('Cis', self)
        self.tab2.cis.stateChanged.connect(self.enableControl)
        self.tab2.cis.stateChanged.connect(self.disableMaxDist)
        self.tab2.linear = QCheckBox('Linear', self)
        self.tab2.linear.stateChanged.connect(self.enableControl)
        self.tab2.linear.stateChanged.connect(self.disableMaxDist)

    def addChargeStates(self):
        """
        Called by self.createTab2ParameterWidgets to create the five charge state checkboxes for charge states of 1+
        through to 5+.
        """
        self.tab2.chargeLabel = QLabel('Charge states (z): ')
        self.tab2.plusOne = QCheckBox('+1', self)
        self.tab2.plusOne.stateChanged.connect(self.enableControl)
        self.tab2.plusTwo = QCheckBox('+2', self)
        self.tab2.plusTwo.stateChanged.connect(self.enableControl)
        self.tab2.plusThree = QCheckBox('+3', self)
        self.tab2.plusThree.stateChanged.connect(self.enableControl)
        self.tab2.plusFour = QCheckBox('+4', self)
        self.tab2.plusFour.stateChanged.connect(self.enableControl)
        self.tab2.plusFive = QCheckBox('+5', self)
        self.tab2.plusFive.stateChanged.connect(self.enableControl)

    def addModifications(self):
        """
        Called by self.createTab2ParameterWidgets to create the initial three modification combo boxes. An additional
        three combos are added later if the user clicks the button to do so.
        :return:
        """
        # Modifications combo boxes and labels
        self.tab2.mod1 = QLabel('Modification 1 : ')
        self.tab2.mod2 = QLabel('Modification 2 : ')
        self.tab2.mod3 = QLabel('Modification 3 : ')
        self.tab2.maxMod = QLabel('Max Mods Per Pep: ')
        self.tab2.mod1Combo = QComboBox(self)
        self.tab2.mod1Combo.activated[str].connect(self.modSelected)
        self.tab2.mod2Combo = QComboBox(self)
        self.tab2.mod2Combo.activated[str].connect(self.modSelected)
        self.tab2.mod3Combo = QComboBox(self)
        self.tab2.mod3Combo.activated[str].connect(self.modSelected)
        self.tab2.maxModCombo = QComboBox(self)
        self.tab2.addMods = QPushButton("Add Additional Mods")
        self.tab2.addMods.clicked.connect(self.addMods)

        # Adding values to modification combo boxes
        self.tab2.mod1Combo.addItem("None")
        self.tab2.mod2Combo.addItem("None")
        self.tab2.mod3Combo.addItem("None")
        for modification in modTable.keys():
            self.tab2.mod1Combo.addItem(modification)
            self.tab2.mod2Combo.addItem(modification)
            self.tab2.mod3Combo.addItem(modification)
        self.tab2.mod1Combo.addItem('Custom Modification')
        self.tab2.mod2Combo.addItem('Custom Modification')
        self.tab2.mod3Combo.addItem('Custom Modification')

        # Add values to maxMod combo box:
        self.tab2.maxModCombo.addItem('None')
        for i in range(1,6):
            self.tab2.maxModCombo.addItem(str(i))

    def addMods(self):
        """
        Called on click of the self.tab2.addMods button. This function adds an additional 3 modification combo boxes
        and corresponding labels to tab 2. It also deletes self.tab2.addMods and adds a different button,
        self.tab2.removeMods, to the GUI in its place.
        :return:
        """
        self.addModsFlag = True
        self.tab2.removeMods = QPushButton("Remove Additional Mods")
        self.tab2.removeMods.clicked.connect(self.removeMods)

        self.tab2.layout.removeWidget(self.tab2.addMods)
        self.tab2.addMods.deleteLater()

        self.addTab2ParameterWidgets()

        # create three new modifications and add them to the GUI
        self.createNewMods()

    def createNewMods(self):
        """
        Called by self.addMods(), this function creates 3 new mods, fills their values and adds them to tab2 of the
        GUI.
        :return:
        """
        # Modifications combo boxes and labels
        self.tab2.mod4 = QLabel('Modification 4 : ')
        self.tab2.mod5 = QLabel('Modification 5 : ')
        self.tab2.mod6 = QLabel('Modification 6 : ')
        self.tab2.mod4Combo = QComboBox(self)
        self.tab2.mod4Combo.activated[str].connect(self.modSelected)
        self.tab2.mod5Combo = QComboBox(self)
        self.tab2.mod5Combo.activated[str].connect(self.modSelected)
        self.tab2.mod6Combo = QComboBox(self)
        self.tab2.mod6Combo.activated[str].connect(self.modSelected)

        # Store the values of the other modification boxes so that they are not added to the list.
        currentMods = [self.tab2.mod1Combo.currentText(), self.tab2.mod2Combo.currentText(), self.tab2.mod3Combo.currentText()]

        # Adding values to modification combo boxes
        self.tab2.mod4Combo.addItem("None")
        self.tab2.mod5Combo.addItem("None")
        self.tab2.mod6Combo.addItem("None")
        for modification in modTable.keys():
            if modification in currentMods:
                continue
            self.tab2.mod4Combo.addItem(modification)
            self.tab2.mod5Combo.addItem(modification)
            self.tab2.mod6Combo.addItem(modification)
        self.tab2.mod4Combo.addItem('Custom Modification')
        self.tab2.mod5Combo.addItem('Custom Modification')
        self.tab2.mod6Combo.addItem('Custom Modification')

        # Add the new mods to the layout
        self.tab2.layout.addWidget(self.tab2.mod4, 7, 3)
        self.tab2.layout.addWidget(self.tab2.mod4Combo, 7, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod5, 8, 3)
        self.tab2.layout.addWidget(self.tab2.mod5Combo, 8, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod6, 9, 3)
        self.tab2.layout.addWidget(self.tab2.mod6Combo, 9, 4, 1, 3)

    def removeMods(self):
        """
        Called on click of self.removeMods, this function deletes the additional three modifications from the GUI. It
        also deletes the self.removeMods button from the GUI and replaces it with the self.addMods button.
        :return:
        """
        # change addModsFlag to False and reset the window
        self.addModsFlag = False
        self.tab2.layout.removeWidget(self.tab2.removeMods)
        self.tab2.removeMods.deleteLater()

        self.tab2.addMods = QPushButton("Add Additional Mods")
        self.tab2.addMods.clicked.connect(self.addMods)

        self.addTab2ParameterWidgets()

        widgetList = [self.tab2.mod4, self.tab2.mod4Combo, self.tab2.mod5, self.tab2.mod5Combo, self.tab2.mod6, \
                      self.tab2.mod6Combo]

        # delete all widgets
        for widget in widgetList:
            self.tab2.layout.removeWidget(widget)
            widget.deleteLater()
            widget = None

    def setDefaultParameters(self):

        """
        Called when MyTableWidget is initialised. This sets the default values of the program to the most common input
        required by users.
        """
        maxDistIndex = self.tab2.maxDistCombo.findText(str(self.maxDistDefault))
        self.tab2.maxDistCombo.setCurrentIndex(maxDistIndex)

        maxIndex = self.tab2.maximumCombo.findText(str(self.maxDefault))
        self.tab2.maximumCombo.setCurrentIndex(maxIndex)

        minIndex = self.tab2.minimumCombo.findText(str(self.minDefault))

        # set to true as defaults for linear, cis and overlap off. Set trans off for now.
        self.tab2.minimumCombo.setCurrentIndex(minIndex)
        self.tab2.overlap.setChecked(True)
        self.tab2.concat.setChecked(True)
        self.tab2.cis.setChecked(True)
        self.tab2.linear.setChecked(True)
        self.tab2.trans.setEnabled(True)
        self.tab2.plusTwo.setChecked(True)
        self.tab1.byIonFlag.setChecked(True)

    # Second set of functions are those that run when the generate output is hit.

    def confirmationFunction(self):

        """
        Called on click of generate output button on tab2. Checks to ensure all input values are relevant and outputs
        message box summarising the inputs of the user. When yes is clicked on the message box, self.getOutputPath()
        is called, prompting the user to select a file location and name before the generation of spliced peptides
        begins.
        """

        # save the input variables as in the MyTableWidget class so that they can be accessed by all methods in the GUI.
        self.ppmVal, self.intensityThreshold, self.mined, self.maxed, self.maxDistance, self.overlapFlag, self.concatFlag, \
        self.transFlag, self.cisFlag, self.linearFlag, self.csvFlag, self.pepToProtFlag, self.protToPepFlag,  \
        self.modList, self.maxMod, self.outputFlag, self.chargeFlags, self.minSimBy, self.byIonAccuracy, \
        self.byIonFlag, self.useMgf = self.getInputParams()

        # if transFlag is selected, we check the size of the input to avoid the user unknowingly starting a huge computation.
        if self.transFlag:
            strng = ""
            maxAminos = 2000
            counter = 0
            for inputFile in self.fasta.inputFile:
                with open(inputFile, "rU") as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        # add to strng and counter
                        counter += 1
                        strng += str(record.seq)
                        # we do not want to measure the length of strng until we are sure that more than one protein has been uploaded.
                        # thus, we continue if counter == 1.
                        if counter == 1:
                            continue
                        # if running trans and the number of aminos in the fasta exceeds 2000, block input.
                        if len(strng) > maxAminos:
                            response = QMessageBox.question(self, 'Message', 'You have selected to compute trans splicing on a file containing over ' +
                            str(maxAminos) + ' amino acids. We do not recommend you persist with this input as it is likely to take a very long time to compute.' +
                            ' Do you still wish to continue with the input?')
                            if response == QMessageBox.Yes:
                                break
                            else:
                                return
            # block output if counter is less than two
            if counter < 2:
                QMessageBox.about(self, "Message", 'Trans output requires at least two proteins to have been uploaded! Please review the input accordingly.')
                return

        self.settingString = 'Minimum Peptide Length: ' + str(self.mined) + '\n' + \
                                     'Maximum Peptide Length: ' + str(self.maxed) + '\n' + \
                                     'Maximum Distance: ' + str(self.maxDistance) + '\n' + \
                                     'Modifications: ' + str(self.modList) + '\n' + \
                                     'Max Mods Per Pep: ' + str(self.maxMod) + '\n' + \
                                     'No Overlap: ' + str(self.overlapFlag) + '\n' + \
                                     'Linear Splicing: ' + str(self.linearFlag) + '\n' + \
                                     'Cis Splicing: ' + str(self.cisFlag) + '\n' + \
                                     'Trans Splicing: ' + str(self.transFlag) + '\n' + \
                                     'Print Intial Combinations: ' + str(self.csvFlag) + '\n' + \
                                     'Write Peptide to Protein Fasta: ' + str(self.pepToProtFlag) + '\n' + \
                                     'Write Protein to Peptide Fasta: ' + str(self.pepToProtFlag) + '\n' + \
                                     'Charge States: ' + str(self.chargeFlags) + '\n' + \
                                     'No MGF Comparison: ' + str(self.useMgf) + '\n' + \
                                     'PPM Value: ' + str(self.ppmVal) + '\n' + \
                                     'Intensity Threshold: ' + str(self.intensityThreshold) + '\n' + \
                                     'Apply b/y Ion Comparison: ' + str(self.byIonFlag) + '\n' + \
                                     'Min b/y Ion %: ' + str(self.minSimBy) + '\n' + \
                                     'b/y Ion Accuracy: ' + str(self.byIonAccuracy) + '\n'

        reply = QMessageBox.question(self, 'Message', 'Do you wish to confirm the following input?\n' + self.settingString,
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            # if the user confirms the input, we ned to run the getOutputPath function. It requires the user input the desired save folder
            # and name the output file. Once the user hits generate output on the file name dialog, the GUI initiates the splicing code.
            self.getOutputPath()

    def getInputParams(self):
        """
        Called by self.confirmationFunction(), this function simply returns all the values input into the GUI widgets
        when the user instructs the program to generate output.

        :return mined: the minimum length of an output peptide.
        :returnmaxed: the maximum length of an ouptut peptide.
        :return overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
        cis spliced peptides.
        :return concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.
        :return transFlag: if True, trans splicing should be completed as part of the output.
        :return cisFlag: if True, cis splicing should be completed as part of the output.
        :return linearFlag: if True, linear splicing should be completed as part of the output.
        :return csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :return pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :return protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath.
        :return modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :return maxMod: the max number of modifications allowable per peptide.
        :return maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
        that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
        :return outputPath: a dictionary holding the output file name of a trans, cis and linear splice output. The
        keys are the TRANS, LINEAR and CIS flags defined at the top of this script.
        :return chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :return mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.
        """
        byIonFlag = self.tab1.byIonFlag.isChecked()
        mgfFlag = self.mgfFlag.isChecked()

        if mgfFlag == False:
            ppmVal = float(self.tab1.ppmText.text())
            toleranceLevel = float(self.tab1.toleranceText.text())
            if byIonFlag == True:
                minByIon = int(self.tab1.minByIonText.text())
                byIonAccuracy = float(self.tab1.byIonAccText.text())
            else:
                minByIon = None
                byIonAccuracy = None
        else:
            minByIon = None
            byIonAccuracy = None
            ppmVal = None
            toleranceLevel = None

        mined = int(self.tab2.minimumCombo.currentText())
        maxed = int(self.tab2.maximumCombo.currentText())
        maxDistance = self.tab2.maxDistCombo.currentText()

        overlapFlag = self.tab2.overlap.isChecked()
        concatFlag = self.tab2.concat.isChecked()
        transFlag = self.tab2.trans.isChecked()
        cisFlag = self.tab2.cis.isChecked()
        linearFlag = self.tab2.linear.isChecked()

        csvFlag = self.tab2.csv.isChecked()
        pepToProtFlag = self.tab2.pepToProt.isChecked()
        protToPepFlag = self.tab2.protToPep.isChecked()

        outputFlag = cisFlag or linearFlag or transFlag

        plusOneFlag = self.tab2.plusOne.isChecked()
        plusTwoFlag = self.tab2.plusTwo.isChecked()
        plusThreeFlag = self.tab2.plusThree.isChecked()
        plusFourFlag = self.tab2.plusFour.isChecked()
        plusFiveFlag = self.tab2.plusFive.isChecked()

        chargeFlags = [plusOneFlag, plusTwoFlag, plusThreeFlag, plusFourFlag, plusFiveFlag]

        modList = [self.tab2.mod1Combo.currentText(), self.tab2.mod2Combo.currentText(),
                   self.tab2.mod3Combo.currentText()]
        if self.addModsFlag:
            modList += [self.tab2.mod4Combo.currentText(), self.tab2.mod5Combo.currentText(),
                        self.tab2.mod6Combo.currentText()]
        print(modList)
        maxMod = self.tab2.maxModCombo.currentText()

        return ppmVal, toleranceLevel, mined, maxed, maxDistance, overlapFlag, concatFlag, transFlag, cisFlag, \
               linearFlag, csvFlag, pepToProtFlag, protToPepFlag, modList, maxMod, outputFlag, chargeFlags, minByIon, \
               byIonAccuracy, byIonFlag, mgfFlag

    def getOutputPath(self):

        """
        Called after generate output is clicked and the user confirms their input. Opens a window to select a file location
        to save the output to, and if valid opens a window to input the file name.
        """
        # opens a window to select file location.
        self.outputPath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        # if no outout path is returned, simply return to the main GUI and the user can choose to recommence the file location
        # selection process if they desire.
        if self.outputPath == '':
            return
        # else if a valid path is selected, bring up a dialog to input the file name
        else:
            self.filePathDialog()

    def filePathDialog(self):
        """
        This function initialises and shows the file naming pop-up box. It is called by self.getOutputPath() if the
        user selects a valid output path.
        """
        self.outputNameBox = QGroupBox('Output Name')
        self.outputNameLayout = QFormLayout()
        self.outputNameLayout.addRow(QLabel("Add a name for the output file."))
        self.outputNameLayout.addRow(QLabel('Banned characters: \ / : * " < > |'))
        self.fileName = QLineEdit()
        self.fileName.textChanged[str].connect(self.nameChecker)
        self.button = QPushButton("Create Output")
        self.valid = QLabel("Valid")
        self.button.clicked.connect(self.returnPath)
        self.outputNameLayout.addRow(self.fileName, self.valid)
        self.outputNameLayout.addRow(self.button)
        self.outputNameBox.setLayout(self.outputNameLayout)
        self.outputNameBox.show()

    def nameChecker(self, input):
        """
        This function is called every time the file name lineEdit is updated. It takes the param input, which is the
        text in the lineEdit, and checks if it is a valid file name.

        :param input: the text in the line-edit which has to be confirmed as valid/invalid.
        """
        # assign bannedCharacters to variables.
        bannedCharacters = set('\/:*"<>|')
        # if the input has no intersection with the banned characters it is valid. If so, update the label validity label
        # and set ensure the generate output button is concurrently enabled/disabled.
        if len(set(input).intersection(bannedCharacters)) == 0:
            self.valid.setText("Valid")
            self.button.setEnabled(True)
        else:
            self.valid.setText("Invalid")
            self.button.setEnabled(False)

    def returnPath(self):
        """
        Called when create output button in the file name dialog is clicked. It takes self.outputPath and adds the
        name input by the user. It then creates the specific names for lin/cis/trans, adds the time and adds these
        output paths to a dictionary. It lastly calls the code in Mers.py to initiate the creation of spliced
        peptides.
        :return:
        """
        # create the base output file name which will be used to create the specific names for lin/cis/tras
        outputFile = self.outputPath + '/' + self.fileName.text()
        print(outputFile)
        # initialise the dictionary to store the splice type file names.
        outputFiles = {}
        now = datetime.now().strftime("%d%m%y_%H%M")

        # write settings to file
        settingPath = outputFile + '_' + 'Info' + now + '.txt'
        file = open(settingPath, 'w')
        file.write('SETTINGS' + '\n')
        file.write('FASTA Files: ' + str(self.fasta.inputFile) + '\n')
        # write to mgfPath to info file if one has been uploaded
        if self.mgfPath is not None:
            file.write('MGF File: ' + str(self.mgfPath) + '\n')
        # write all the settings to file
        file.write(self.settingString)

        # create the file name for each splice type and add to dictionary.
        if self.linearFlag:
            linPath = outputFile + '_' + LINEAR + now + ".fasta"
            outputFiles[LINEAR] = Path(linPath)
        if self.cisFlag:
            cisPath = outputFile + '_' + CIS + now + ".fasta"
            outputFiles[CIS] = Path(cisPath)
        if self.transFlag:
            transPath = outputFile + '_' + TRANS + now + ".fasta"
            outputFiles[TRANS] = Path(transPath)
            # print(outputPath[TRANS])

        # disable widgets here before the mgf upload commences
        self.disableWidgets()

        # if the mgfFlag is not checked, mgfGen imports the mgf then runs importedMGF() which generates the rest of the output.
        # note that all the input variables have been declared under the myTableWidget class and are thus callable from
        # this function.
        if self.mgfFlag.isChecked() == False:
            mgfGen = MGFImporter(self.uploadMgf, self.mgfPath, self.ppmVal, self.intensityThreshold, self.minSimBy,
                                 self.byIonAccuracy, self.byIonFlag, self.chargeFlags)
            mgfGen.signals.finished.connect(functools.partial(self.importedMGF, self.mined, self.maxed, self.overlapFlag,
                                                              self.concatFlag, self.transFlag, self.cisFlag, self.linearFlag, self.csvFlag,
                                                              self.pepToProtFlag, self.protToPepFlag, self.modList, self.maxMod,
                                                              self.maxDistance, outputFiles, self.chargeFlags))
            # add label informing the user that the mgf is uploading
            self.mgfLabel = QLabel("Uploading MGF. Please Wait!")
            self.tab2.layout.addWidget(self.mgfLabel, 19, 3, 1, 2)

            self.threadpool.start(mgfGen)
        # if mgfFlag is checked, no need to import the mgf, can skip straight to running importedMGF()
        else:
            self.importedMGF(self.mined, self.maxed, self.overlapFlag, self.concatFlag, self.transFlag, self.cisFlag, self.linearFlag, self.csvFlag, self.pepToProtFlag,
                             self.protToPepFlag, self.modList, self.maxMod, self.maxDistance, outputFiles, self.chargeFlags, True)

        # close the output name box.
        self.outputNameBox.close()

    def uploadMgf(self, input_path, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag, chargeFlags):
        """
        Called by self.returnPath() as the function to be run by MGFImporter if an MGF file has been uploaded by the
        user. This function extract the data needed for the splicing program from the MGF file and stores it in the
        object self.mgf.

        :param input_path: the path of the MGF file uploaded by the user.
        :param ppmVal: the value input by the user which details how close the monoisotpic mass a peptide must be to
        the precursor mass of an MGF spectra for them to be considered a match. Measured in parts per million.
        :param intensityThreshold: after the max intensity of each spectrum is stored and sorted, only those above
        the user input intensityThreshold are considered for comparison with the spliced peptides.
        :param minSimBy: the percentage of b/y ions that have a matching ion in the potential spectrum for the peptide
        and the spectrum to be considered a match.
        :param byIonAccuracy: how close the b/y ion m/z and ion m/z as per spectrum must be for the given b/y ion to
        be considered a match.
        :param byIonFlag: True if the user requires b/y ion comparison to occur.
        :param chargeFlags: the charge states from the MGF file that the user has selected to consider.
        :return:
        """
        mgfDfList, pepmassIonArrayList, mgfLen = readMGF(input_path, intensityThreshold, byIonFlag, chargeFlags)
        print(mgfDfList)

        maxMass, chargeMaxDict = self.maxMgfMass(mgfDfList, chargeFlags)

        self.mgf = MGF(mgfDfList, pepmassIonArrayList, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag,
                       maxMass, chargeMaxDict, mgfLen)

    def maxMgfMass(self, mgfDfList, chargeFlags):
        """
        Called from self.uploadMGF(), this returns the max mono-isotopic mass present in the MGF file and also the
        max m/z ratio present for each individual charge state.

        :param mgfDfList: a list of dictionaries, when each dictionary contains charges as keys and precursor masses
        from the MGF file with the given charge state. The list will only contain more than one of these dictionaries
        if the data was too large to pass into a multiprocessing.Process() as one dictionary.
        :param chargeFlags: the charge states that the user wishes to asses in the MGF file.

        :return maxMass: the max mono-isotopic mass (precursor mass must be converted from m/z to mono-isotopic)
        present in the MGF file given the selected charge states.
        :return chargeMaxDict: a dictionary containing charge states as keys and the corresponding max m/z ratio in
        the MGF file as values.
        """

        maxMass = 0
        chargeMaxDict = {}
        for mgfDf in mgfDfList:
            for z, masses in mgfDf.items():
                if chargeFlags[int(z) - 1]:
                    maxChargeMass = max(masses)

                    chargeMaxDict[z] = maxChargeMass
                    maxMassTemp = maxChargeMass * int(z) - int(z) * 1.00794
                    if maxMassTemp > maxMass:
                        maxMass = maxMassTemp

        return maxMass, chargeMaxDict

    def importedMGF(self, mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                    protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag=False):
        """
        Called by self.returnPath() once the output path has been created and all MGF data (if any) has been
        imported. This function then calls self.outputPreStep() which is the first step in initiating the code in
        Mer.py which creates the spliced peptides.

        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
        cis spliced peptides.
        :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.
        :param transFlag: if True, trans splicing should be completed as part of the output.
        :param cisFlag: if True, cis splicing should be completed as part of the output.
        :param linearFlag: if True, linear splicing should be completed as part of the output.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath.
        :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :param maxMod: the max number of modifications allowable per peptide.
        :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
        that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
        :param outputPath: a dictionary holding the output file name of a trans, cis and linear splice output. The
        keys are the TRANS, LINEAR and CIS flags defined at the top of this script.
        :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.

        :return:
        """

        print("MGF FILE UPLOADED")

        try:
            self.tab2.layout.removeWidget(self.mgfLabel)
            self.mgfLabel.deleteLater()
            self.mgfLabel = None
        except AttributeError:
            pass

        self.outputPreStep(mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                           protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag)

    def outputPreStep(self, mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                      protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag):

        """
        Called by self.importedMGF(), this function begins the spliced peptide generation by creating a new thread
        to keep the GUI responsive. It creates two threads, on via OutputGenerator and one via ProgressGenerator.
        The OutputGenerator thread runs the function self.output() while the ProgressGenerator thread controls the
        creating and updating of a progress bar.

        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
        cis spliced peptides.
        :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.
        :param transFlag: if True, trans splicing should be completed as part of the output.
        :param cisFlag: if True, cis splicing should be completed as part of the output.
        :param linearFlag: if True, linear splicing should be completed as part of the output.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath.
        :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :param maxMod: the max number of modifications allowable per peptide.
        :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
        that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
        :param outputPath: a dictionary holding the output file name of a trans, cis and linear splice output. The
        keys are the TRANS, LINEAR and CIS flags defined at the top of this script.
        :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.
        :return:
        """

        self.progressLabel = QLabel('Collating Combinations. Please Wait: ')
        self.tab2.layout.addWidget(self.progressLabel, 19, 3, 1, 2)
        self.progressBar = QProgressBar(self)
        self.tab2.layout.addWidget(self.progressBar, 20, 3, 1, 4)

        self.outputGen = OutputGenerator(self.output, mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag,
                                         csvFlag, pepToProtFlag, protToPepFlag,
                                         modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag)

        self.outputGen.signals.finished.connect(self.outputFinished)
        self.threadpool.start(self.outputGen)

        self.progressBarUpdate = ProgressGenerator()
        self.progressBarUpdate.signals.updateProgBar.connect(self.updateProgressBar)
        self.progressBarUpdate.signals.finished.connect(self.deleteTab2ProgressBar)
        #self.progressBarUpdate.signals.disableButtons.connect(self.disableWidgets)
        self.threadpool.start(self.progressBarUpdate)

    def output(self, mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag, protToPepFlag,
               modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag):

        """
        Sent as the worker function to an OutputGenerator thread, it runs the self.fasta.generateOutput() function
        contained within Mers.py. This ultimately begins the generation of spliced peptides.

        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
        cis spliced peptides.
        :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.
        :param transFlag: if True, trans splicing should be completed as part of the output.
        :param cisFlag: if True, cis splicing should be completed as part of the output.
        :param linearFlag: if True, linear splicing should be completed as part of the output.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath.
        :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :param maxMod: the max number of modifications allowable per peptide.
        :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
        that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
        :param outputPath: a dictionary holding the output file name of a trans, cis and linear splice output. The
        keys are the TRANS, LINEAR and CIS flags defined at the top of this script.
        :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.

        :return:
        """

        start = time.time()

        #self.tab2.stop.setEnabled(True)

        if maxDistance != 'None':
            maxDistance = int(maxDistance)

        self.fasta.generateOutput(mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                                  protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, self.mgf, mgfFlag)
        end = time.time()

        #self.tab2.stop.setEnabled(False)

        self.emptyProgQueues()

        print(end - start)

    def outputFinished(self):

        """
        Called once the OutputGenerator created in self.outputPreStep() has finished running self.output(). This
        function instructs the progress bar QRunnable to finish, enables all widgets on the window and creates a
        message box instructing the user that the output has finished.
        """
        print("IT'S DONE")
        # send flag to prog bar instructing it to break the update loop and finish running.
        self.progressBarUpdate.changeFlag()
        # enable all widgets and produce popup
        self.enableAllWidgets()
        QMessageBox.about(self, "Message", 'Output Complete')

    def emptyProgQueues(self):
        """
        Called after self.fasta.generateOutput() has finished, and thus the program has finished compliling the
        spliced peptides for the current user input. It empties any elements in the queue self.fasta.pepCompleted,
        and then clears the variables self.fasta.completedProcs, self.fasta.procGenCounter and self.fasta.totalProcs.
        If these variables aren't cleared, the progress bar and process generation will not run correctly if the user
        wished to run another splicing.
        :return:
        """
        while not self.fasta.pepCompleted.empty():
            clearQ = self.fasta.pepCompleted.get()
        self.fasta.completedProcs = 0
        self.fasta.procGenCounter = 0
        self.fasta.totalProcs = 0

    def updateProgressBar(self):
        """
        This function checks how many peptides have been completed. It does so by communicating with functions in
        Mers.py via the multiprocessing.Queue self.fasta.pepCompleted(). This function is called every 10ms to check
        if a 1 has been added to the queue, signifying that the output of an input protein has been completed. If
        a 1 has been added, the progress bar value is updated.
        :return:
        """
        if not self.fasta.pepCompleted.empty():
            self.fasta.completedProcs += self.fasta.pepCompleted.get()

        # plus 1 to avoid divide by 0 error
        value = self.fasta.completedProcs/(self.fasta.totalProcs+1)*100
        self.progressBar.setValue(value)
        #print(self.fasta.completedProcs)

    def deleteTab2ProgressBar(self):
        """
        This function is called when thread self.outputGen (a GenerateOutput instance) has finished after being
        initialised in self.outputPreStep(). The function removes the progress bar and the corresponding label from
        the window and deletes the reference to it.
        :return:
        """
        # Delete progress label and progress bar
        self.tab1.layout.removeWidget(self.progressLabel)
        self.progressLabel.deleteLater()
        self.progressLabel = None
        self.tab1.layout.removeWidget(self.progressBar)
        self.progressBar.deleteLater()
        self.progressBar = None

    # Functions called on click of upload fasta, upload mgf and next tab buttons

    def uploadMgfPreStep(self):
        """
        Called from the select MGF file button. Opens a window to select a file, and check if the file ends in MGF.
        It also creates a plot of the max intensities of each spectrum so that the user can make an informed decision
        regarding the intensity threshold input.
        :return:
        """

        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home')
        mgfTest = fname[0][-3:]

        if mgfTest == 'mgf':
            self.mgfPath = fname[0]
            if self.mgfPlotFlag.isChecked():
                self.progressLabel = QLabel('Creating Intensity Plot! Please Wait.')
                self.tab1.layout.addWidget(self.progressLabel, 10, 2, 1, 2)
                self.disableWidgets()
                self.mgfPlot = MGFPlotter(plotData, fname[0])
                self.mgfPlot.signals.plot.connect(self.onlyImportMGF)
                self.mgfPlot.signals.finished.connect(self.intensityPlotFin)
                self.threadpool.start(self.mgfPlot)
            else:
                QMessageBox.about(self, "Message", 'MGF file successfully selected!')
            self.enableControl()

        # Ensuring program does not crash if no file is selected
        elif fname[0] == '':
            print('')

        # Wrong extension selected! Try Again!
        else:
            QMessageBox.about(self, "Message", 'Please select a MGF file!')

    def onlyImportMGF(self, ms2Thresh, intensityPoints):
        """
        Called by uploadMgfPreStep() when the user has selected to produce the intensity plot. In this case, the MGF
        file is opened and a plot of the max intensities of each spectrum is produced and shown to the user.

        :param ms2Thresh: the intensity threshold values that are to be plotted on the x-axis
        :param intensityPoints: the percentage of spectra with an intensity above the intensity threshold of the same
        index. This is to be plotted on the y-axis.
        :return:
        """

        plot(ms2Thresh, intensityPoints)

    def intensityPlotFin(self):
        """
        Called when self.onlyImportMGF() finishes being run by the MGFPlotter thread. This function simply
        closes the progress bar by calling the changeFlag() function, and then enables all widgets as per
        self.enableControl().
        :return:
        """
        self.tab1.layout.removeWidget(self.progressLabel)
        self.progressLabel.deleteLater()
        self.progressLabel = None

        self.enableTab2Widgets()
        self.pushButton1.setEnabled(True)
        self.mgfButton.setEnabled(True)
        self.mgfPlotFlag.setEnabled(True)
        QMessageBox.about(self, "Message", 'MGF file successfully selected!')
        self.enableControl()

    def uploadFasta(self):
        """
        Called on click of the Upload Fasta File button. Opens a window to select a file, and check if the file ends in
        '.fasta'.

        :return:
        """
        sender = self.sender()

        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home/')

        fastaTest = fname[0][-5:]

        # Ensure opening fasta extension file by checking last five chars
        if fastaTest == 'fasta':
            if sender == self.pushButton1:
                self.fasta = Fasta(fname[0])
                #self.enableControl()
                self.controlMGFInput()
                QMessageBox.about(self, "Message", 'Fasta file imported.')
            else:
                self.fasta.inputFile.append(fname[0])
                QMessageBox.about(self, "Message", 'Additional Fasta file imported.')
            #print(self.fasta.inputFile)

        # Ensuring program does not crash if no file is selected
        elif fname[0] == '':
            pass

        # Wrong extension selected! Try Again!
        else:
            QMessageBox.about(self, "Message", 'Please select a Fasta file!')

    def nextTabFunc(self):
        """
        Called on click of the button self.nextTab, this function simply navigates to tab 2 of the window.
        :return:
        """
        self.tabs.setCurrentIndex(1)

    # ENABLE/DISABLE FUNCTIONS: control flow of program and which inputs are available to the user given the
    # inputs they have inserted.

    def disableWidgets(self):
        """
        Called from self.returnPath() before the peptide splicing computation begins. This function simply disables
        all input on the GUI until the output has finished.
        """
        self.tab2.minimumCombo.setEnabled(False)
        self.tab2.maximumCombo.setEnabled(False)
        self.tab2.maxDistCombo.setEnabled(False)
        self.tab2.overlap.setEnabled(False)
        self.tab2.concat.setEnabled(False)
        self.tab2.trans.setEnabled(False)
        self.tab2.cis.setEnabled(False)
        self.tab2.linear.setEnabled(False)
        self.tab2.plusOne.setEnabled(False)
        self.tab2.plusTwo.setEnabled(False)
        self.tab2.plusThree.setEnabled(False)
        self.tab2.plusFour.setEnabled(False)
        self.tab2.plusFive.setEnabled(False)
        self.tab2.mod1Combo.setEnabled(False)
        self.tab2.mod2Combo.setEnabled(False)
        self.tab2.mod3Combo.setEnabled(False)
        self.tab2.maxModCombo.setEnabled(False)
        self.pushButton1.setEnabled(False)
        self.mgfButton.setEnabled(False)
        self.mgfPlotFlag.setEnabled(False)
        self.mgfFlag.setEnabled(False)
        self.nextTab.setEnabled(False)
        self.tab1.ppmText.setEnabled(False)
        self.tab1.toleranceText.setEnabled(False)
        self.tab1.minByIonText.setEnabled(False)
        self.tab1.byIonAccText.setEnabled(False)
        self.tab1.byIonFlag.setEnabled(False)
        self.tab2.csv.setEnabled(False)
        self.tab2.pepToProt.setEnabled(False)
        self.tab2.protToPep.setEnabled(False)
        self.tab2.output.setEnabled(False)

    def enableAllWidgets(self):
        """
        Called from self.outputFinished() once the peptide splicing computation has been completed. The function
        enables all GUI widgets/inputs so the user can edit them and start another output if they desire.
        :return:
        """
        self.enableTab1Widgets()
        self.enableTab2Widgets()
        self.disableMaxDist()
        if self.mgfFlag.isChecked():
            self.controlMGFInput()

    def enableTab1Widgets(self):
        """
        Called from self.enableAllWidgets(), this function simply enables all the widgets on tab 1.
        :return:
        """
        self.pushButton1.setEnabled(True)
        self.mgfButton.setEnabled(True)
        self.mgfPlotFlag.setEnabled(True)
        self.mgfFlag.setEnabled(True)
        self.nextTab.setEnabled(True)
        self.tab1.ppmText.setEnabled(True)
        self.tab1.toleranceText.setEnabled(True)
        self.tab1.byIonFlag.setEnabled(True)
        if self.tab1.byIonFlag.isChecked() == True:
            self.tab1.minByIonText.setEnabled(True)
            self.tab1.byIonAccText.setEnabled(True)

    def enableTab2Widgets(self):
        """
        Called from self.enableAllWidgets(), this function simply enables all the widgets on tab 2.
        :return:
        """
        self.tab2.minimumCombo.setEnabled(True)
        self.tab2.maximumCombo.setEnabled(True)
        self.tab2.maxDistCombo.setEnabled(True)
        self.tab2.overlap.setEnabled(True)
        self.tab2.concat.setEnabled(True)
        self.tab2.trans.setEnabled(True)
        self.tab2.cis.setEnabled(True)
        self.tab2.linear.setEnabled(True)
        self.tab2.mod1Combo.setEnabled(True)
        self.tab2.mod2Combo.setEnabled(True)
        self.tab2.mod3Combo.setEnabled(True)
        self.tab2.maxModCombo.setEnabled(True)
        self.tab2.csv.setEnabled(True)
        self.tab2.pepToProt.setEnabled(True)
        self.tab2.protToPep.setEnabled(True)
        self.tab2.output.setEnabled(True)
        # only enable the charge checkboxes if mgfFlag is False
        if not self.mgfFlag.isChecked():
            self.tab2.plusOne.setEnabled(True)
            self.tab2.plusTwo.setEnabled(True)
            self.tab2.plusThree.setEnabled(True)
            self.tab2.plusFour.setEnabled(True)
            self.tab2.plusFive.setEnabled(True)
            self.tab2.chargeLabel.setEnabled(True)

    def firstTabValid(self):
        """
        Called by self.enableControl(), thus function checks that the user has input sufficient information to the
        first tab for them to be allowed to move to the second.
        :return:
        """
        if self.mgfFlag.isChecked() == True:
            return True
        if self.tab1.byIonFlag.isChecked():
            statusList = [self.tab1.byIonAccStatus.text(), self.tab1.ppmStatus.text(), \
                          self.tab1.toleranceStatus.text(), self.tab1.minByIonStatus.text()]
        else:
            statusList = [self.tab1.ppmStatus.text(), self.tab1.toleranceStatus.text()]
        for status in statusList:
            if status in ["Invalid", ""]:
                return False
        return True

    def secondTabValid(self):
        """
        Called by self.enableControl(), this function checks that all information in the second tab has a been
        correctly input by the user.
        :return:
        """
        transFlag = self.tab2.trans.isChecked()
        cisFlag = self.tab2.cis.isChecked()
        linearFlag = self.tab2.linear.isChecked()
        outputFlag = cisFlag or linearFlag or transFlag

        plusOneFlag = self.tab2.plusOne.isChecked()
        plusTwoFlag = self.tab2.plusTwo.isChecked()
        plusThreeFlag = self.tab2.plusThree.isChecked()
        plusFourFlag = self.tab2.plusFour.isChecked()
        plusFiveFlag = self.tab2.plusFive.isChecked()
        chargeFlag = plusOneFlag or plusTwoFlag or plusThreeFlag or plusFourFlag or plusFiveFlag

        if chargeFlag and outputFlag:
            return True
        return False

    def enableControl(self):
        """
        This function is called when an input widgets in the GUI is altered by the user. It checks what has been
        input and systematically allows the user access to the proceeding input. This ensures the user fills out all
        inputs and doesn't run the program without inputting all required data.

        :return:
        """
        # if a fasta file has been uploaded, enable tha self.addMultibleFasta button
        if self.fasta is not None:
            self.addMultipleFasta.setEnabled(True)
            # if an MGF file has been uploaded, or the user has selected not to select an MGF file, we enable access
            # to the inputs relevant to the MGF data.
            if self.mgfPath is not None or self.mgfFlag.isChecked() == True:
                self.tab1.toleranceText.setEnabled(True)
                self.tab1.toleranceLabel.setEnabled(True)
                self.tab1.ppmText.setEnabled(True)
                self.tab1.ppmLabel.setEnabled(True)
                self.tab1.byIonFlag.setEnabled(True)
                # if the b/y ion flag is checked, enable the b/y ion data. Otherwise, we disable them.
                if self.tab1.byIonFlag.isChecked():
                    self.tab1.byIonAccText.setEnabled(True)
                    self.tab1.byIonAccLabel.setEnabled(True)
                    self.tab1.minByIonText.setEnabled(True)
                    self.tab1.minByIonLabel.setEnabled(True)
                else:
                    self.tab1.byIonAccText.setEnabled(False)
                    self.tab1.byIonAccLabel.setEnabled(False)
                    self.tab1.minByIonText.setEnabled(False)
                    self.tab1.minByIonLabel.setEnabled(False)
                # if the first tab is valid, we enable the next tab. Otherwise, we ensure it is disabled.
                if self.firstTabValid():
                    self.tabs.setTabEnabled(1, True)
                    self.nextTab.setEnabled(True)
                    # if the second tab is valid, we enable the output button, otherwise we turn it off.
                    if self.secondTabValid():
                        self.tab2.output.setEnabled(True)
                    else:
                        self.tab2.output.setEnabled(False)
                else:
                    self.tabs.setTabEnabled(1, False)
                    self.nextTab.setEnabled(False)
            else:
                self.tabs.setTabEnabled(1, False)
                self.nextTab.setEnabled(False)
        else:
            self.addMultipleFasta.setEnabled(False)

    def controlMGFInput(self):
        """
        When self.mgfFlag is changed, this function is called before self.enableControl() to ensure all the MGF related inputs
        are either enabled or disabled prior to self.enableControl() being run.
        :return:
        """
        if self.mgfFlag.isChecked():
            # set values true/false before calling enableControl to avoid bug of it being called again later
            self.tab1.byIonFlag.setChecked(False)
            self.enableControl()
            self.mgfButton.setEnabled(False)
            self.mgfPlotFlag.setEnabled(False)
            self.tab1.ppmText.setEnabled(False)
            self.tab1.ppmLabel.setEnabled(False)
            self.tab1.toleranceText.setEnabled(False)
            self.tab1.toleranceLabel.setEnabled(False)
            self.tab1.minByIonText.setEnabled(False)
            self.tab1.minByIonLabel.setEnabled(False)
            self.tab1.byIonAccText.setEnabled(False)
            self.tab1.byIonAccLabel.setEnabled(False)
            self.tab1.byIonFlag.setEnabled(False)
            self.tab2.plusOne.setEnabled(False)
            self.tab2.plusTwo.setEnabled(False)
            self.tab2.plusThree.setEnabled(False)
            self.tab2.plusFour.setEnabled(False)
            self.tab2.plusFive.setEnabled(False)
            self.tab2.chargeLabel.setEnabled(False)
            #self.tab2.trans.setEnabled(True)
            # self.tab1.ppmLabel.setEnabled(False)
        if not self.mgfFlag.isChecked():
            # set values true/false before calling enableControl to avoid bug of it being called again later
            self.tab1.byIonFlag.setChecked(True)
            self.tab2.trans.setChecked(False)
            self.mgfButton.setEnabled(True)
            self.mgfPlotFlag.setEnabled(True)
            #self.tab2.trans.setEnabled(False)
            self.tab2.plusOne.setEnabled(True)
            self.tab2.plusTwo.setEnabled(True)
            self.tab2.plusThree.setEnabled(True)
            self.tab2.plusFour.setEnabled(True)
            self.tab2.plusFive.setEnabled(True)
            self.tab2.chargeLabel.setEnabled(True)
            self.tab1.ppmText.setEnabled(False)
            self.tab1.ppmLabel.setEnabled(False)
            self.tab1.toleranceText.setEnabled(False)
            self.tab1.toleranceLabel.setEnabled(False)
            self.enableControl()

    def disableMaxDist(self):
        """
        Called when the checkbox self.tab2.cis changes state. It ensures that the max distance combo box and the
        overlap checkbox are disabled when cis hasn't been selected, and enabled when it has.
        """
        if self.tab2.cis.isChecked():
            self.tab2.maxDistCombo.setEnabled(True)
            self.tab2.overlap.setEnabled(True)
            self.tab2.overlap.setEnabled(True)
        else:
            index = self.tab2.maxDistCombo.findText('None')
            self.tab2.maxDistCombo.setCurrentIndex(index)
            self.tab2.maxDistCombo.setEnabled(False)
            self.tab2.overlap.setChecked(True)
            self.tab2.overlap.setEnabled(False)

    def minMaxChanged(self, text):

        """
        Called when self.tab2.minimumCombo or self.tab2.maximumCombo values are changes. It alters the values available
        in the min, max and maxDistance combos to ensure a realistic input. That is, that min cannot be higher than max
        and max cannot be higher than maxDistance.
        """

        sender = self.tab2.sender()
        minChanged = sender == self.tab2.minimumCombo

        if minChanged:
            comboChange = self.tab2.maximumCombo
        else:
            comboChange = self.tab2.minimumCombo

        # current Value of combo to be changed
        value = int(comboChange.currentText())

        # Current Max Distance - convert 'None' to 0 so it can be used as a comparator later.
        maxDistValue = self.tab2.maxDistCombo.currentText()
        if maxDistValue == 'None':
            maxDistInt = 0
        else:
            maxDistInt = int(maxDistValue)

        # Clear combo box values for combo to be changed
        comboChange.clear()

        # Creates new values in max combo box which are greater than the min
        if minChanged:
            for i in range(int(text) - 1, 26):
                comboChange.addItem(str(i + 1))
            # Restores current value if it is greater than the min
            if value >= int(text):
                indexMax = comboChange.findText(str(value))
                comboChange.setCurrentIndex(indexMax)

        else:  # maxChanged
            # Creates new values in combo box which are less than than the max
            for i in range(2, int(text) + 1):
                comboChange.addItem(str(i))
            # Restores current value if it is less than the max
            if value <= int(text):
                indexMin = comboChange.findText(str(value))
                comboChange.setCurrentIndex(indexMin)

    def modSelected(self, text):
        """
        Called each time a modification is selected from the mod combo boxes. This function ensures that a mod selected
        in one combo box will not be available for selection in any of the others.
        :param text: the modification selected by the user, that is to be removed from the options in the other combo
        boxes.
        :return:
        """

        # if addmodsFlag as is checked, we have 6 mods not three. Update modCombos accordingly.
        if self.addModsFlag:
            modCombos = [self.tab2.mod1Combo, self.tab2.mod2Combo, self.tab2.mod3Combo, self.tab2.mod4Combo,
                         self.tab2.mod5Combo, self.tab2.mod6Combo]
        else:
            modCombos = [self.tab2.mod1Combo, self.tab2.mod2Combo, self.tab2.mod3Combo]

        modSender = []
        modChange = []

        # append True in element of modSender if corresponding element in modCombos is the sender.
        for combo in modCombos:
            modSender.append(combo == self.tab2.sender())

        # modChange contains all the combos which were not changed (not the sender).
        for i in range(0, len(modSender)):
            if not modSender[i]:
                modChange.append(modCombos[i])

        sender = self.tab2.sender()
        if sender.currentText() == 'Custom Modification':
            self.showCustomMod(sender)

        modValues = []
        for i in range(0, len(modChange)):
            modValues.append(modChange[i].currentText())
            modChange[i].clear()
            modChange[i].addItem('None')

        for modification in modTable.keys():
            # and modification != modValue1 and modification != modValue2:
            if modification != text and modification not in modValues:
                for i in range(0, len(modChange)):
                    modChange[i].addItem(modification)

        for i in range(0, len(modChange)):
            modChange[i].addItem("Custom Modification")

        for i in range(0, len(modChange)):
            if modValues[i] not in (text, 'None'):
                modChange[i].addItem(modValues[i])
                indexMod = modChange[i].findText(modValues[i])
                modChange[i].setCurrentIndex(indexMod)

    def showCustomMod(self, sender):
        """
        Called from self.modSelected() if the user has selected Custom Modification from one of the drop down boxes.
        This function creates a pop-up box which promts the user to input the details of a new modification that
        they wish to apply.
        :param sender: the name of the combo box which the user selected Custom Modification from.
        :return:
        """
        sender.setCurrentIndex(0)
        self.formGroupBox = QGroupBox('Custom Modification')
        self.formLayout = QFormLayout()
        self.formLayout.addRow(QLabel("Input modified amino acids without spaces: TGN \n" +
                                      "Input mass change as a decimal number. Place a minus (-) \n" +
                                        "directly before to signify mass loss."))
        self.custAminoInput = QLineEdit()
        self.custMassInput = QLineEdit()
        self.modName = QLineEdit()
        self.addModButton = QPushButton("Create Modification")
        # partial method allows variable to be passed to connected function on click
        self.addCustToModListSender = partial(self.addCustToModlist, sender)
        self.addModButton.clicked.connect(self.addCustToModListSender)
        self.formLayout.addRow(QLabel("Modified Amino Acids: "), self.custAminoInput)
        self.formLayout.addRow(QLabel("Mass Change: "), self.custMassInput)
        self.formLayout.addRow(QLabel("Modification Name: "), self.modName)
        self.formLayout.addRow(self.addModButton)
        self.formGroupBox.setLayout(self.formLayout)
        self.formGroupBox.show()

    def addCustToModlist(self, sender):
        """
        Called when self.addModButton is clicked from the custom modification pop-up window. This function checks that
        the inputs by the user are valid and then adds the modification to the global modTable variable.
        :param sender: the name of the combo box which the user selected Custom Modification from.
        :return:
        """
        # strip inputs so that leading or lagging whitespace does not void the tests.
        aminoAcids = self.custAminoInput.text().strip()
        massChange = self.custMassInput.text().strip()
        modName = self.modName.text().strip()
        modKey = "Custom " + modName + " (" + aminoAcids + ")"
        modValue = []

        # validity checks
        for char in aminoAcids:
            if char not in monoAminoMass.keys():
                QMessageBox.about(self, "Message",
                                  'Amino Acids are not valid. Please leave no spaces and use capitals!')
                return
            else:
                modValue.append(char)

        # if not float (accounts for +/- at start) return error message.
        try:
            float(massChange[1:])
        except ValueError:
            QMessageBox.about(self, "Message", 'Mass Change is not a valid decimal number!')
            return


        # update modTable and GUI
        modValue.append(float(massChange))
        modTable[modKey] = modValue

        sender.addItem(modKey)
        indexCustom = sender.findText(modKey)
        sender.setCurrentIndex(indexCustom)

        # delete pop up box
        self.formLayout.removeWidget(self.formGroupBox)
        self.formGroupBox.deleteLater()
        self.formGroupBox = None

        QMessageBox.about(self, "Message", 'Custom Modification successfully added!')

        return True

    def textBoxChanged(self, input):
        """
        Called when the input to any of the line edits on the first tab are changed. This function checks if the input
        to the altered line edit is valid and updates the validity labels next the text box to reflect this state.
        :param input: the text in the line edit that was changed.
        :return:
        """
        label, min, max = self.textBoxSender(self.sender())
        try:
            byIonAcc = float(input)
            if min <= float(input) and max >= float(input):
                label.setText("Valid")
            else:
                label.setText("Invalid")
        except ValueError:
            if input == "":
                label.setText("")
            else:
                label.setText("Invalid")

    def textBoxSender(self, sender):
        """
        Called from self.textBoxChanged(), this function returns the min and max values allowed in the input line edit,
        and also the validity label related to the input line edit.
        :param sender: the line edit which is to be checked for validity in self.textBoxChanged().
        :return:
        """
        if sender == self.tab1.byIonAccText:
            label = self.tab1.byIonAccStatus
            min = 0.01
            max = 0.2
        elif sender == self.tab1.ppmText:
            label = self.tab1.ppmStatus
            min = 0.1
            max = 1000
        elif sender == self.tab1.toleranceText:
            label = self.tab1.toleranceStatus
            min = 0
            max = 10000
        elif sender == self.tab1.minByIonText:
            label = self.tab1.minByIonStatus
            min = 0
            max = 100
        return label, min, max

    @pyqtSlot()
    def on_click(self):
        # print("\n")
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())


if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = QApplication(sys.argv)
    if len(sys.argv) == 1:
        ex = App()
        sys.exit(app.exec_())
    else:
        print('command line interface')





