# pyinstaller MersGUI --> this command from the relevant file location creates executable file
# Need to fix no parameter
import sys
import subprocess
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
from Mers import *
from MGFMain import *
import functools


class WorkerSignals(QObject):
    """
    Signals class that is used for the GUI when emitting custom signals
    """

    finished = pyqtSignal()
    updateProgBar = pyqtSignal(int)
    disableButtons = pyqtSignal()


class ProgressGenerator(QRunnable):
    """
    Progress Bar that shows up under the charge states once the peptide slicing begins
    """

    def __init__(self, *args, **kwargs):
        super(ProgressGenerator, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.flag = True

    def changeFlag(self):
        self.flag = False

    @pyqtSlot()
    def run(self):
        self.signals.disableButtons.emit()
        while self.flag:
            self.signals.updateProgBar.emit(0)
            time.sleep(0.5)
            self.signals.updateProgBar.emit(25)
            time.sleep(0.5)
            self.signals.updateProgBar.emit(50)
            time.sleep(0.5)
            self.signals.updateProgBar.emit(75)
            time.sleep(0.5)
            self.signals.updateProgBar.emit(100)
            time.sleep(0.5)
        self.signals.finished.emit()


class OutputGenerator(QRunnable):
    """

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
        self.fn(*self.args)
        self.signals.finished.emit()


class MGFImporter(QRunnable):
    """

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
        start = time.time()
        self.fn(*self.args)
        self.signals.finished.emit()
        end = time.time()
        print("Uploading mgf took: " + str(end - start))


class App(QMainWindow):
    """
    App serves as the parent class for the embedded MyTableWidget
    """

    def __init__(self):
        """
         Initialisation of main window class
        """

        super().__init__()
        self.title = 'Peptide Splicer'
        self.fastaTest = False
        self.outputPath = ""
        self.statusbar = self.statusBar()

        self.center()

        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)
        self.setWindowTitle(self.title)

        self.show()

    def center(self):
        """
        center function is called to centre the main window on the screen
        """

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def closeEvent(self, event):
        print('closed')
        sys.exit()


class MyTableWidget(QWidget):
    """
    MyTableWidget class is the child of the App class. It holds the tabs where most of the GUI functionality occurs
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

        # Init threading
        self.threadpool = QThreadPool()

        # Default values for the input parameters
        self.minDefault = '8'
        self.maxDefault = '12'
        self.maxDistDefault = '25'
        self.minByIonDefault = '50'
        self.byIonAccDefault = '0.1'

        # Initialisation of two tabs
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tabs.resize(300, 200)

        # Add tabs to table class (self)
        self.tabs.addTab(self.tab1, "Select File and Path")
        self.tabs.addTab(self.tab2, "Input Parameters")

        # Creation of tab layout and widgets within tab
        self.tab1.layout = QGridLayout(self)
        self.tab1.layout.setSpacing(10)
        self.createTab1ParameterWidgets()
        self.addTab1ParameterWidgets()

        self.tab1.setLayout(self.tab1.layout)

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



    def importedMGF(self):
        print("MGF FILE UPLOADED")
        QMessageBox.about(self, "Message", 'MGF file imported.')



    def uploadMgf(self, input_path, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag):
        print("UPLOADING MGF")
        mgfDf, pepmassIonArray = readMGF(input_path, intensityThreshold)

        self.mgf = MGF(mgfDf, pepmassIonArray, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag)


    def uploadMgfPreStep(self):
        """
        Called from the select Fasta file button. Opens a window to select a file, and check if the file ends in MGF
        """

        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home')
        mgfTest = fname[0][-3:]

        if mgfTest == 'mgf':
            self.mgfPath = fname[0]

            # mgfGen = MGFImporter(self.uploadMgf, fname[0])
            #
            # mgfGen.signals.finished.connect(self.importedMGF)
            # self.threadpool.start(mgfGen)


        # Ensuring program does not crash if no file is selected
        elif fname[0] == '':
            print('')

        # Wrong extension selected! Try Again!
        else:
            QMessageBox.about(self, "Message", 'Please select a MGF file!')

    def uploadFasta(self):

        """
        Called from the Upload Fasta File button. Opens a window to select a file, and check if the file ends in fasta
        """

        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home/')

        fastaTest = fname[0][-5:]

        # Ensure opening fasta extension file by checking last five chars
        if fastaTest == 'fasta':
            self.fasta = Fasta(addSequenceList(fname[0]))

            QMessageBox.about(self, "Message", 'Fasta file imported.'
                                               'There are ' + str(self.fasta.entries) + ' proteins in this file!')

        # Ensuring program does not crash if no file is selected
        elif fname[0] == '':
            print('')

        # Wrong extension selected! Try Again!
        else:
            QMessageBox.about(self, "Message", 'Please select a Fasta file!')

    def getOutputPath(self):

        """
        Called after generate output is clicked. Opens a window to select a file location to save the output to.
        Returns False if no path is selected, otherwise returns the selected path.
        """

        outputPath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        if outputPath == '':
            return False
        return outputPath

    def stopFunction(self):
        print('in stop function')
        for process in self.fasta.allProcessList:
            process.terminate()

    def confirmationFunction(self):

        """
        Called on click of generate output button on tab2. Checks to ensure all input values are relevant and outputs
        message box summarising the inputs of the user. When yes is clicked on the message box, the output function is
        called which begins generating results
        """

        if self.tab1.byIonAccStatus.text() in ["Invalid", ""] or self.tab1.ppmStatus.text() in ["Invalid",""]:
            QMessageBox.about(self, "Message", 'Please check that valid PPM and B/Y Ion Accuracy '
                                               'values have been selected')
        else:
            ppmVal, intensityThreshold, mined, maxed, maxDistance, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, \
            modList, outputFlag, chargeFlags, minSimBy, byIonAccuracy, byIonFlag = self.getInputParams()

        if self.fasta is None:

            QMessageBox.about(self, "Message", 'Please check that a valid Fasta file and output '
                                               'file location have been selected')

        elif not outputFlag:
            QMessageBox.about(self, "Message", 'Please select at least one output type; either trans, cis or linear')
        else:
            reply = QMessageBox.question(self, 'Message', 'Do you wish to confirm the following input?\n' +
                                         'Minimum Peptide Length: ' + str(mined) + '\n' +
                                         'Maximum Peptide Length: ' + str(maxed) + '\n' +
                                         'Maximum Distance: ' + str(maxDistance) + '\n' +
                                         'Mod List: ' + str(modList) + '\n' +
                                         'Overlap Flag: ' + str(overlapFlag) + '\n' +
                                         'Linear Flag: ' + str(linearFlag) + '\n' +
                                         'Cis Flag: ' + str(cisFlag) + '\n' +
                                         'Trans Flag: ' + str(transFlag) + '\n' +
                                         'CSV Flag: ' + str(csvFlag) + '\n' +
                                         'Charge States: ' + str(chargeFlags) + '\n' +
                                         'PPM Value: ' + str(ppmVal) + '\n' +
                                         'Intensity Threshold' + str(intensityThreshold) + '\n'
                                         'Min b/y Ion (%): ' + str(minSimBy) + '\n' +
                                         'b/y Ion Accuracy: ' + str(byIonAccuracy) + '\n' +
                                         'b/y Ion Flag: ' + str(byIonFlag),
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.Yes:

                print(self.mgfPath)
                # UPLOAD MGF FILE HERE
                mgfGen = MGFImporter(self.uploadMgf, self.mgfPath, ppmVal, intensityThreshold, minSimBy,
                                     byIonAccuracy, byIonFlag)

                if csvFlag:
                    outputPath = self.getOutputPath()
                    if outputPath is not False:
                        self.outputPreStep(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, modList,
                                           maxDistance, outputPath, chargeFlags)
                else:

                    outputPath = None
                    #mgfGen.signals.finished.connect(self.onlyImportMGF)
                    mgfGen.signals.finished.connect(functools.partial(self.importedMGF, mined, maxed, overlapFlag,
                                                                      transFlag, cisFlag, linearFlag, csvFlag, modList,
                                                                      maxDistance, outputPath, chargeFlags))
                    self.threadpool.start(mgfGen)

    def onlyImportMGF(self):
        print("MGF FILE Uploaded")

    def importedMGF(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, modList,
                    maxDistance, outputPath, chargeFlags):

        print("MGF FILE UPLOADED")


        self.outputPreStep(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, modList,
                           maxDistance, outputPath, chargeFlags)

    def finished(self):

        """
        Alerts when done!
        """
        print("IT'S DONE")

        self.progressBarUpdate.changeFlag()
        self.tab2.output.setEnabled(True)
        self.pushButton1.setEnabled(True)
        QMessageBox.about(self, "Message", 'Output Complete')

    def updateProgressBar(self, value):
        self.tab2.progressBar.setValue(value)

    def deleteProgressBar(self):
        # Delete progress label and progress bar
        self.tab2.layout.removeWidget(self.tab2.progressLabel)
        self.tab2.progressLabel.deleteLater()
        self.tab2.progressLabel = None
        self.tab2.layout.removeWidget(self.tab2.progressBar)
        self.tab2.progressBar.deleteLater()
        self.tab2.progressBar = None

    def disableButtons(self):
        """
        These buttons should not be being used when an output is being generated.
        """
        self.tab2.output.setEnabled(False)
        self.pushButton1.setEnabled(False)

    def outputPreStep(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, modList, maxDistance,
                      outputPath, chargeFlags):

        """
        Begins the output by creating a threadpool to keep gui responsive. Called by the confirmation function; also
        calls the actual output function to generate output
        """

        self.tab2.progressLabel = QLabel('Collating Combinations. Please Wait: ')
        self.tab2.layout.addWidget(self.tab2.progressLabel, 15, 3, 1, 2)
        self.tab2.progressBar = QProgressBar(self)
        self.tab2.layout.addWidget(self.tab2.progressBar, 16, 3, 1, 4)

        self.outputGen = OutputGenerator(self.output, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag,
                                         csvFlag,
                                         modList, maxDistance, outputPath, chargeFlags)

        self.outputGen.signals.finished.connect(self.finished)
        self.threadpool.start(self.outputGen)

        self.progressBarUpdate = ProgressGenerator()
        self.progressBarUpdate.signals.updateProgBar.connect(self.updateProgressBar)
        self.progressBarUpdate.signals.finished.connect(self.deleteProgressBar)
        self.progressBarUpdate.signals.disableButtons.connect(self.disableButtons)
        self.threadpool.start(self.progressBarUpdate)

    def disableMaxDist(self, state):
        """
        Called when trans is selected, it disables the use of the max distance function
        :param state: Whether trans flag is True or False
        :return:
        """
        if state == Qt.Checked:
            index = self.tab2.maxDistCombo.findText('None')
            self.tab2.maxDistCombo.setCurrentIndex(index)
            self.tab2.maxDistCombo.setEnabled(False)
        else:
            self.tab2.maxDistCombo.setEnabled(True)

    def output(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, modList,
               maxDistance, outputPath, chargeFlags):

        """
        called by output pre-step function, it runs the generateOutput function from Mers.py; This is shown to be
        running via the progress bar
        """

        start = time.time()

        if maxDistance != 'None':
            maxDistance = int(maxDistance)

        print("ABOUT TO DO OUTPUT")
        self.fasta.generateOutput(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, modList,
                                  maxDistance, outputPath, chargeFlags, self.mgf)
        end = time.time()

        print(end - start)

        # The following statements are used to open the output directory after output is created into file format
        # print(tryString)
        # replacedOutpath = outputPath.replace("/", '\\')
        # print(replacedOutpath)
        # openString = r'explorer "' + replacedOutpath + '"'
        # print(openString)
        # subprocess.Popen(openString)

    def minMaxChanged(self, text):

        """
        Data Validation Function.
        Called when minimumCombo value changes. It alters the values available in max and maxDistance combos to
        ensure a realistic input
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

            # wipe max distance values only if max is changed. Add None back to combo box.
            self.tab2.maxDistCombo.clear()
            self.tab2.maxDistCombo.addItem('None')
            # refill max distance combo box with items greater than or equal to max
            for i in range(int(text) - 1, 26):
                self.tab2.maxDistCombo.addItem(str(i + 1))
            if maxDistInt >= int(text):
                indexDist = self.tab2.maxDistCombo.findText(str(maxDistValue))
                self.tab2.maxDistCombo.setCurrentIndex(indexDist)

    def modSelected(self, text):

        """
        Ensures only one mod can be selected.
        """

        modCombos = [self.tab2.mod1Combo, self.tab2.mod2Combo, self.tab2.mod3Combo]
        modSender = []
        modChange = []

        for combo in modCombos:
            modSender.append(combo == self.tab2.sender())

        for i in range(0, len(modSender)):
            if not modSender[i]:
                modChange.append(modCombos[i])

        modValue1 = modChange[0].currentText()
        modValue2 = modChange[1].currentText()

        modChange[0].clear()
        modChange[1].clear()
        modChange[0].addItem('None')
        modChange[1].addItem('None')

        for modification in modTable.keys():
            # and modification != modValue1 and modification != modValue2:
            if modification not in (text, modValue1, modValue2):
                modChange[0].addItem(modification)
                modChange[1].addItem(modification)

        if modValue1 not in (text, 'None'):
            modChange[0].addItem(modValue1)
            indexMod1 = modChange[0].findText(modValue1)
            modChange[0].setCurrentIndex(indexMod1)

        if modValue2 not in (text, 'None'):
            modChange[1].addItem(modValue2)
            indexMod2 = modChange[1].findText(modValue2)
            modChange[1].setCurrentIndex(indexMod2)

    def getInputParams(self):

        ppmVal = int(self.tab1.ppmText.text())
        toleranceLevel = float(self.tab1.toleranceCombo.currentText())

        mined = int(self.tab2.minimumCombo.currentText())
        maxed = int(self.tab2.maximumCombo.currentText())
        maxDistance = self.tab2.maxDistCombo.currentText()

        overlapFlag = self.tab2.overlap.isChecked()
        transFlag = self.tab2.trans.isChecked()
        cisFlag = self.tab2.cis.isChecked()
        linearFlag = self.tab2.linear.isChecked()

        csvFlag = self.tab2.csv.isChecked()

        outputFlag = cisFlag or linearFlag or transFlag

        plusOneFlag = self.tab2.plusOne.isChecked()
        plusTwoFlag = self.tab2.plusTwo.isChecked()
        plusThreeFlag = self.tab2.plusThree.isChecked()
        plusFourFlag = self.tab2.plusFour.isChecked()
        plusFiveFlag = self.tab2.plusFive.isChecked()

        chargeFlags = [plusOneFlag, plusTwoFlag, plusThreeFlag, plusFourFlag, plusFiveFlag]

        # self.fasta = Fasta(addSequenceList('/Users/nicolaschapman/Documents/UROP/Code/MersProject/small.fasta'))
        # self.fasta = Fasta(addSequenceList('C:/Users/Arpit/Desktop/UROP/InputData/oneProtein.fasta'))
        # self.fasta = Fasta(addSequenceList('C:/Users/Administrator/Desktop/UROP/InputData/OneProtein.fasta'))
        # self.mgf = MGF(readMGF('C:/Users/Arpit/Desktop/UROP/InputData/mgf.mgf'))

        modList = [self.tab2.mod1Combo.currentText(), self.tab2.mod2Combo.currentText(),
                   self.tab2.mod3Combo.currentText()]


        minByIon = int(self.tab1.minByIonCombo.currentText())
        byIonAccuracy = float(self.tab1.byIonAccText.text())
        byIonFlag = self.tab1.byIonFlag.isChecked()

        return ppmVal, toleranceLevel, mined, maxed, maxDistance, overlapFlag, transFlag, cisFlag, \
               linearFlag, csvFlag, modList, outputFlag, chargeFlags, minByIon, byIonAccuracy, byIonFlag

    def addMinMaxAndDist(self):

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
        for i in range(int(self.maxDefault), 26):
            self.tab2.maxDistCombo.addItem(str(i))

    def addFlagChecks(self):

        # initialise overlap, trans, cis and linear check boxes
        self.tab2.overlap = QCheckBox('Overlap Off', self)
        self.tab2.trans = QCheckBox('Trans', self)
        self.tab2.trans.stateChanged.connect(self.disableMaxDist)  # connect trans check box to relevant function
        self.tab2.cis = QCheckBox('Cis', self)
        self.tab2.linear = QCheckBox('Linear', self)

    def addChargeStates(self):

        """
        Add charge state check boxes
        """
        self.tab2.chargeLabel = QLabel('Charge states (z): ')
        self.tab2.plusOne = QCheckBox('+1', self)
        self.tab2.plusTwo = QCheckBox('+2', self)
        self.tab2.plusThree = QCheckBox('+3', self)
        self.tab2.plusFour = QCheckBox('+4', self)
        self.tab2.plusFive = QCheckBox('+5', self)

    def addModifications(self):

        # Modifications combo boxes and labels
        self.tab2.mod1 = QLabel('Modification 1 : ')
        self.tab2.mod2 = QLabel('Modification 2 : ')
        self.tab2.mod3 = QLabel('Modification 3 : ')
        self.tab2.mod1Combo = QComboBox(self)
        self.tab2.mod1Combo.activated[str].connect(self.modSelected)
        self.tab2.mod2Combo = QComboBox(self)
        self.tab2.mod2Combo.activated[str].connect(self.modSelected)
        self.tab2.mod3Combo = QComboBox(self)
        self.tab2.mod3Combo.activated[str].connect(self.modSelected)

        # Adding values to modification combo boxes
        self.tab2.mod1Combo.addItem("None")
        self.tab2.mod2Combo.addItem("None")
        self.tab2.mod3Combo.addItem("None")
        for modification in modTable.keys():
            self.tab2.mod1Combo.addItem(modification)
            self.tab2.mod2Combo.addItem(modification)
            self.tab2.mod3Combo.addItem(modification)

    def ppmChanged(self, input):
        try:
            ppmVal = float(input)
            self.tab1.ppmStatus.setText("Valid")
        except ValueError:
            if input == "":
                self.tab1.ppmStatus.setText("")
            else:
                self.tab1.ppmStatus.setText("Invalid")


    def byIonAccChanged(self, input):
        try:
            byIonAcc = float(input)
            self.tab1.byIonAccStatus.setText("Valid")
        except ValueError:
            if input == "":
                self.tab1.byIonAccStatus.setText("")
            else:
                self.tab1.byIonAccStatus.setText("Invalid")

    def createTab1ParameterWidgets(self):
        self.pushButton1 = QPushButton("Select Fasta File")
        self.pushButton1.clicked.connect(self.uploadFasta)
        self.mgfButton = QPushButton("Select MGF File")
        self.mgfButton.clicked.connect(self.uploadMgfPreStep)

        self.tab1.ppmLabel = QLabel('PPM (decimal number): ')
        self.tab1.ppmText = QLineEdit(self)
        self.tab1.ppmText.textChanged[str].connect(self.ppmChanged)
        self.tab1.ppmStatus = QLabel("")

        self.tab1.toleranceLabel = QLabel('Intensity Threshold: ')
        self.tab1.toleranceCombo = QComboBox(self)

        self.tab1.minByIonLabel = QLabel('Minimum b/y Ion Matches(%): ')
        self.tab1.minByIonCombo = QComboBox(self)

        self.tab1.byIonAccLabel = QLabel('b/y Ion Accuracy (decimal number): ')
        self.tab1.byIonAccText = QLineEdit(self)
        self.tab1.byIonAccText.textChanged[str].connect(self.byIonAccChanged)
        self.tab1.byIonAccStatus = QLabel("")

        self.tab1.byIonFlag = QCheckBox('Apply b/y Ion Comparison: ')

        # for i in range(10, 110, 10):
        #     self.tab1.ppmCombo.addItem(str(i))

        intensities = [0, 10, 50, 100, 500, 1000, 5000, 10000]
        for intensity in intensities:
            self.tab1.toleranceCombo.addItem(str(intensity))

        for i in range(10, 100, 10):
            self.tab1.minByIonCombo.addItem(str(i))

        # ionAccuracies = [0.4, 0.2, 0.1, 0.05, 0.02, 0.01]
        # for accuracy in ionAccuracies:
        #     self.tab1.byIonAccCombo.addItem(str(accuracy))

    def addTab1ParameterWidgets(self):
        self.tab1.layout.setColumnStretch(0, 1)
        self.tab1.layout.setColumnStretch(5, 1)
        self.tab1.layout.setRowStretch(0, 1)
        self.tab1.layout.setRowStretch(5, 1)
        self.tab1.layout.addWidget(self.pushButton1, 1, 2, 1, 2)
        self.tab1.layout.addWidget(self.mgfButton, 2, 2, 1, 2)
        self.tab1.layout.addWidget(self.tab1.ppmLabel, 3, 2)
        self.tab1.layout.addWidget(self.tab1.ppmText, 3, 3)
        self.tab1.layout.addWidget(self.tab1.ppmStatus, 3, 4)
        self.tab1.layout.addWidget(self.tab1.toleranceLabel, 4, 2)
        self.tab1.layout.addWidget(self.tab1.toleranceCombo, 4, 3)

        self.tab1.layout.addWidget(self.tab1.minByIonLabel, 5, 2)
        self.tab1.layout.addWidget(self.tab1.minByIonCombo, 5, 3)
        self.tab1.layout.addWidget(self.tab1.byIonAccLabel, 6, 2)
        self.tab1.layout.addWidget(self.tab1.byIonAccText, 6, 3)
        self.tab1.layout.addWidget(self.tab1.byIonAccStatus, 6, 4)
        self.tab1.layout.addWidget(self.tab1.byIonFlag, 7, 2)

    def createTab2ParameterWidgets(self):

        self.addMinMaxAndDist()
        self.addModifications()
        self.addFlagChecks()
        self.addChargeStates()

        # AN EXTRA ADD "WRITE TO CSV FUNCTION CHECKBOX"
        self.tab2.csv = QCheckBox('Write To Csv')

        # create generate output push button
        self.tab2.output = QPushButton('Generate Output!', self)
        self.tab2.stop = QPushButton('Stop Process', self)
        self.tab2.output.clicked.connect(self.confirmationFunction)
        self.tab2.stop.clicked.connect(self.stopFunction)

        self.setDefaultParameters()

    def addTab2ParameterWidgets(self):

        """
        Add all widgets in the correct grid format where they should be on the screen
        """

        # All the labels added to grid layout of tab2
        self.tab2.layout.addWidget(self.tab2.minimum, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximum, 2, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistance, 3, 3)
        self.tab2.layout.addWidget(self.tab2.mod1, 4, 3)
        self.tab2.layout.addWidget(self.tab2.mod2, 5, 3)
        self.tab2.layout.addWidget(self.tab2.mod3, 6, 3)
        self.tab2.layout.addWidget(self.tab2.overlap, 7, 3)
        self.tab2.layout.addWidget(self.tab2.linear, 8, 3)
        self.tab2.layout.addWidget(self.tab2.cis, 9, 3)
        self.tab2.layout.addWidget(self.tab2.trans, 10, 3)
        self.tab2.layout.addWidget(self.tab2.csv, 11, 3)
        self.tab2.layout.addWidget(self.tab2.chargeLabel, 12, 3)

        # all dynamic elements added to the grid layout of tab 2
        self.tab2.layout.addWidget(self.tab2.minimumCombo, 1, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximumCombo, 2, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistCombo, 3, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod1Combo, 4, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod2Combo, 5, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod3Combo, 6, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.plusOne, 12, 4)
        self.tab2.layout.addWidget(self.tab2.plusTwo, 12, 5)
        self.tab2.layout.addWidget(self.tab2.plusThree, 12, 6)
        self.tab2.layout.addWidget(self.tab2.plusFour, 13, 4)
        self.tab2.layout.addWidget(self.tab2.plusFive, 13, 5)
        self.tab2.layout.addWidget(self.tab2.output, 14, 5, 1, 2)
        self.tab2.layout.addWidget(self.tab2.stop, 14, 3)

    def setDefaultParameters(self):

        """
        set default values
        """
        maxDistIndex = self.tab2.maxDistCombo.findText(str(self.maxDistDefault))
        self.tab2.maxDistCombo.setCurrentIndex(maxDistIndex)

        maxIndex = self.tab2.maximumCombo.findText(str(self.maxDefault))
        self.tab2.maximumCombo.setCurrentIndex(maxIndex)

        minIndex = self.tab2.minimumCombo.findText(str(self.minDefault))

        # set to true as defaults for linear, cis and overlap off. Set trans off for now.
        self.tab2.minimumCombo.setCurrentIndex(minIndex)
        self.tab2.overlap.setChecked(True)
        self.tab2.cis.setChecked(True)
        self.tab2.linear.setChecked(True)
        self.tab2.trans.setEnabled(False)
        self.tab2.plusTwo.setChecked(True)

        minByIonIndex = self.tab1.minByIonCombo.findText(self.minByIonDefault)
        self.tab1.minByIonCombo.setCurrentIndex(minByIonIndex)

        # byIonAccIndex = self.tab1.byIonAccCombo.findText(self.byIonAccDefault)
        # self.tab1.byIonAccCombo.setCurrentIndex(byIonAccIndex)

        self.tab1.byIonFlag.setChecked(True)

    @pyqtSlot()
    def on_click(self):
        print("\n")
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())


if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())




