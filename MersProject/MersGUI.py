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


class WorkerSignals(QObject):
    """
    Signals class that is used for the GUI when emitting custom signals
    """

    finished = pyqtSignal()
    updateProgBar = pyqtSignal()
    disableButtons = pyqtSignal()
    plot = pyqtSignal(list, list)


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
            self.signals.updateProgBar.emit()
            time.sleep(0.01)
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


class MGFPlotter(QRunnable):
    """

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
        start = time.time()
        ms2Thresh, intensityPoints = self.fn(*self.args)

        self.signals.plot.emit(ms2Thresh, intensityPoints)
        end = time.time()
        print("Getting plot data took: " + str(end - start))
        self.signals.finished.emit()


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
        self.outputPath = []
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

        # variables to count the total number of processes and those which have finished
        self.finishedPeptides = 0
        self.totalSize = 0

    def uploadMgf(self, input_path, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag, chargeFlags):
        mgfDf, pepmassIonArray = readMGF(input_path, intensityThreshold)

        maxMass, chargeMaxDict = self.maxMgfMass(mgfDf, chargeFlags)

        self.mgf = MGF(mgfDf, pepmassIonArray, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag,
                       maxMass, chargeMaxDict)

    def maxMgfMass(self, mgfDf, chargeFlags):
        maxMass = 0
        chargeMaxDict = {}
        for z, masses in mgfDf.items():
            if chargeFlags[int(z)-1]:
                maxChargeMass = max(masses)

                chargeMaxDict[z] = maxChargeMass
                maxMassTemp = maxChargeMass*int(z) - int(z)*1.00794
                if maxMassTemp > maxMass:
                    maxMass = maxMassTemp

        return maxMass, chargeMaxDict



    def onlyImportMGF(self, ms2Thresh, intensityPoints):

        plot(ms2Thresh, intensityPoints)

    def uploadMgfPreStep(self):
        """
        Called from the select Fasta file button. Opens a window to select a file, and check if the file ends in MGF
        """

        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home')
        mgfTest = fname[0][-3:]

        if mgfTest == 'mgf':
            self.mgfPath = fname[0]
            if self.mgfPlotFlag.isChecked():
                self.progressLabel = QLabel('Creating Intensity Plot. Please Wait: ')
                self.tab1.layout.addWidget(self.progressLabel, 8, 2, 1, 2)
                self.progressBar = QProgressBar(self)
                self.tab1.layout.addWidget(self.progressBar, 9, 2, 1, 4)

                self.progressBarUpdate = ProgressGenerator()
                self.progressBarUpdate.signals.updateProgBar.connect(self.updateProgressBar)
                self.progressBarUpdate.signals.finished.connect(self.deleteTab1ProgressBar)
                self.progressBarUpdate.signals.disableButtons.connect(self.disableWidgets)
                self.threadpool.start(self.progressBarUpdate)

                self.mgfPlot = MGFPlotter(plotData, fname[0])
                self.mgfPlot.signals.plot.connect(self.onlyImportMGF)
                self.mgfPlot.signals.finished.connect(self.intensityPlotFin)
                self.threadpool.start(self.mgfPlot)
            self.enableControl()
            QMessageBox.about(self, "Message", 'MGF file successfully selected!')

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
            self.fasta = Fasta(fname[0])
            #self.enableControl()
            self.controlMGFInput()
            QMessageBox.about(self, "Message", 'Fasta file imported.')

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

        outputFile = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        if outputFile == '':
            return False
        else:
            text, ok = QInputDialog.getText(self, 'Input Dialog',
                                            'Enter your file name:')

            if ok:
                outputPath = outputFile + '/' + text
            else:
                return False
        return outputPath

    def stopFunction(self):
        print('in stop function')
        for process in self.fasta.allProcessList:
            process.terminate()

    def nextTabFunc(self):
        self.tabs.setCurrentIndex(1)

    def firstTabValid(self):
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
        if self.fasta is not None:
            if self.mgfPath is not None or self.mgfFlag.isChecked() == True:
                self.tab1.toleranceText.setEnabled(True)
                self.tab1.toleranceLabel.setEnabled(True)
                self.tab1.ppmText.setEnabled(True)
                self.tab1.ppmLabel.setEnabled(True)
                self.tab1.byIonFlag.setEnabled(True)
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
                if self.firstTabValid():
                    self.tabs.setTabEnabled(1, True)
                    self.nextTab.setEnabled(True)
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

    def controlMGFInput(self):
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
            #self.tab2.trans.setEnabled(True)
            # self.tab1.ppmLabel.setEnabled(False)
        if not self.mgfFlag.isChecked():
            # set values true/false before calling enableControl to avoid bug of it being called again later
            self.tab2.trans.setChecked(False)
            self.enableControl()
            self.mgfButton.setEnabled(True)
            self.mgfPlotFlag.setEnabled(True)
            #self.tab2.trans.setEnabled(False)


    def confirmationFunction(self):

        """
        Called on click of generate output button on tab2. Checks to ensure all input values are relevant and outputs
        message box summarising the inputs of the user. When yes is clicked on the message box, the output function is
        called which begins generating results
        """

        ppmVal, intensityThreshold, mined, maxed, maxDistance, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, \
        pepToProtFlag, protToPepFlag,  modList, maxMod, outputFlag, chargeFlags, minSimBy, byIonAccuracy, \
        byIonFlag, mgfFlag = self.getInputParams()

        # if transFlag is selected, we check the size of the input to avoid the user unknowingly starting a huge computation.
        if transFlag:
            with open(self.fasta.inputFile, "rU") as handle:
                strng = ""
                maxAminos = 2000
                for record in SeqIO.parse(handle, 'fasta'):
                    strng += str(record.seq)
                    # if running trans and the number of aminos in the fasta exceeds 2000, block input.
                    if len(strng) > maxAminos:
                        response = QMessageBox.question(self, 'Message', 'You have selected to compute trans splicing on a file containing over ' +
                        str(maxAminos) + ' amino acids. We do not recommend you persist with this input as it is likely to take a very long time to compute.' +
                        ' Do you still wish to continue with the input?')
                        if response == QMessageBox.Yes:
                            break
                        else:
                            return

        reply = QMessageBox.question(self, 'Message', 'Do you wish to confirm the following input?\n' +
                                     'Minimum Peptide Length: ' + str(mined) + '\n' +
                                     'Maximum Peptide Length: ' + str(maxed) + '\n' +
                                     'Maximum Distance: ' + str(maxDistance) + '\n' +
                                     'Modifications: ' + str(modList) + '\n' +
                                     'Max Mods Per Pep: ' + str(maxMod) + '\n' +
                                     'No Overlap: ' + str(overlapFlag) + '\n' +
                                     'Linear Splicing: ' + str(linearFlag) + '\n' +
                                     'Cis Splicing: ' + str(cisFlag) + '\n' +
                                     'Trans Splicing: ' + str(transFlag) + '\n' +
                                     'Print Intial Combinations: ' + str(csvFlag) + '\n' +
                                     'Write Peptide to Protein Fasta: ' + str(pepToProtFlag) + '\n' +
                                     'Write Protein to Peptide Fasta: ' + str(pepToProtFlag) + '\n' +
                                     'Charge States: ' + str(chargeFlags) + '\n' +
                                     'No MGF Comparison: ' + str(mgfFlag) + '\n' +
                                     'PPM Value: ' + str(ppmVal) + '\n' +
                                     'Intensity Threshold: ' + str(intensityThreshold) + '\n' +
                                     'Apply b/y Ion Comparison: ' + str(byIonFlag) + '\n' +
                                     'Min b/y Ion %: ' + str(minSimBy) + '\n' +
                                     'b/y Ion Accuracy: ' + str(byIonAccuracy) + '\n',
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:

            outputPath = {}
            now = datetime.now().strftime("%d%m%y_%H%M")
            outputFile = self.getOutputPath()
            if outputFile is not False:
                if linearFlag:
                    linPath = outputFile + '-' + LINEAR + now + ".fasta"
                    outputPath[LINEAR] = Path(linPath)
                if cisFlag:
                    cisPath = outputFile + '-' + CIS + now + ".fasta"
                    outputPath[CIS] = Path(cisPath)
                if transFlag:
                    transPath = outputFile + '-' + TRANS + now + ".fasta"
                    outputPath[TRANS] = Path(transPath)
                    # print(outputPath[TRANS])

                if self.mgfFlag.isChecked() == False:
                    mgfGen = MGFImporter(self.uploadMgf, self.mgfPath, ppmVal, intensityThreshold, minSimBy,
                                         byIonAccuracy, byIonFlag, chargeFlags)
                    mgfGen.signals.finished.connect(functools.partial(self.importedMGF, mined, maxed, overlapFlag,
                                                                      transFlag, cisFlag, linearFlag, csvFlag,
                                                                      pepToProtFlag, protToPepFlag, modList, maxMod,
                                                                      maxDistance, outputPath, chargeFlags))
                    self.threadpool.start(mgfGen)
                else:
                    self.importedMGF(mined, maxed, overlapFlag,transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                                     protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, True)

    def importedMGF(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                    protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag=False):

        print("MGF FILE UPLOADED")

        self.outputPreStep(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                           protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag)
    def outputFinished(self):

        """
        Alerts when done!
        """
        print("IT'S DONE")

        self.progressBarUpdate.changeFlag()
        self.enableAllWidgets()
        QMessageBox.about(self, "Message", 'Output Complete')

    def intensityPlotFin(self):
        self.progressBarUpdate.changeFlag()
        self.enableTab2Widgets()
        self.pushButton1.setEnabled(True)
        self.mgfButton.setEnabled(True)
        self.mgfPlotFlag.setEnabled(True)
        self.enableControl()

    def updateProgressBar(self):
        if not self.fasta.pepTotal.empty():
            self.totalSize += self.fasta.pepTotal.get()
        if not self.fasta.pepCompleted.empty():
            self.finishedPeptides += self.fasta.pepCompleted.get()
        if self.totalSize is not 0 and self.finishedPeptides is not 0:
            value = self.finishedPeptides/self.totalSize*100
        else:
            value = 2
        self.progressBar.setValue(value)

    def deleteTab2ProgressBar(self):
        # Delete progress label and progress bar
        self.tab1.layout.removeWidget(self.progressLabel)
        self.progressLabel.deleteLater()
        self.progressLabel = None
        self.tab1.layout.removeWidget(self.progressBar)
        self.progressBar.deleteLater()
        self.progressBar = None

    def deleteTab1ProgressBar(self):
        # Delete progress label and progress bar
        self.tab1.layout.removeWidget(self.progressLabel)
        self.progressLabel.deleteLater()
        self.progressLabel = None
        self.tab1.layout.removeWidget(self.progressBar)
        self.progressBar.deleteLater()
        self.progressBar = None

    def disableWidgets(self):
        """
        These buttons should not be being used when an output is being generated.
        """
        self.tab2.minimumCombo.setEnabled(False)
        self.tab2.maximumCombo.setEnabled(False)
        self.tab2.maxDistCombo.setEnabled(False)
        self.tab2.overlap.setEnabled(False)
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
        self.nextTab.setEnabled(False)
        self.tab1.ppmText.setEnabled(False)
        self.tab1.toleranceText.setEnabled(False)
        self.tab1.minByIonText.setEnabled(False)
        self.tab1.byIonAccText.setEnabled(False)
        self.tab1.byIonFlag.setEnabled(False)
        self.tab2.csv.setEnabled(False)
        self.tab2.pepToProt.setEnabled(False)
        self.tab2.protToPep.setEnabled(False)

    def enableAllWidgets(self):

        self.enableTab1Widgets()
        self.enableTab2Widgets()
        self.disableMaxDist()
        if self.mgfFlag.isChecked():
            self.controlMGFInput()


    def enableTab1Widgets(self):
        self.pushButton1.setEnabled(True)
        self.mgfButton.setEnabled(True)
        self.mgfPlotFlag.setEnabled(True)
        self.nextTab.setEnabled(True)
        self.tab1.ppmText.setEnabled(True)
        self.tab1.toleranceText.setEnabled(True)
        self.tab1.byIonFlag.setEnabled(True)
        if self.tab1.byIonFlag.isChecked() == True:
            self.tab1.minByIonText.setEnabled(True)
            self.tab1.byIonAccText.setEnabled(True)

    def enableTab2Widgets(self):
        self.tab2.minimumCombo.setEnabled(True)
        self.tab2.maximumCombo.setEnabled(True)
        self.tab2.maxDistCombo.setEnabled(True)
        self.tab2.overlap.setEnabled(True)
        self.tab2.trans.setEnabled(True)
        self.tab2.cis.setEnabled(True)
        self.tab2.linear.setEnabled(True)
        self.tab2.plusOne.setEnabled(True)
        self.tab2.plusTwo.setEnabled(True)
        self.tab2.plusThree.setEnabled(True)
        self.tab2.plusFour.setEnabled(True)
        self.tab2.plusFive.setEnabled(True)
        self.tab2.mod1Combo.setEnabled(True)
        self.tab2.mod2Combo.setEnabled(True)
        self.tab2.mod3Combo.setEnabled(True)
        self.tab2.maxModCombo.setEnabled(True)
        self.tab2.csv.setEnabled(True)
        self.tab2.pepToProt.setEnabled(True)
        self.tab2.protToPep.setEnabled(False)
        self.tab2.output.setEnabled(True)

    def outputPreStep(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                      protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag):

        """
        Begins the output by creating a threadpool to keep gui responsive. Called by the confirmation function; also
        calls the actual output function to generate output
        """

        self.progressLabel = QLabel('Collating Combinations. Please Wait: ')
        self.tab2.layout.addWidget(self.progressLabel, 16, 3, 1, 2)
        self.progressBar = QProgressBar(self)
        self.tab2.layout.addWidget(self.progressBar, 17, 3, 1, 4)

        self.outputGen = OutputGenerator(self.output, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag,
                                         csvFlag, pepToProtFlag, protToPepFlag,
                                         modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag)

        self.outputGen.signals.finished.connect(self.outputFinished)
        self.threadpool.start(self.outputGen)

        self.progressBarUpdate = ProgressGenerator()
        self.progressBarUpdate.signals.updateProgBar.connect(self.updateProgressBar)
        self.progressBarUpdate.signals.finished.connect(self.deleteTab2ProgressBar)
        self.progressBarUpdate.signals.disableButtons.connect(self.disableWidgets)
        self.threadpool.start(self.progressBarUpdate)

    def disableMaxDist(self):
        """
        Called when trans is selected, it disables the use of the max distance function
        :param state: Whether trans flag is True or False
        :return:
        """
        if self.tab2.trans.isChecked():
            index = self.tab2.maxDistCombo.findText('None')
            self.tab2.maxDistCombo.setCurrentIndex(index)
            self.tab2.maxDistCombo.setEnabled(False)
            self.tab2.overlap.setChecked(True)
            self.tab2.overlap.setEnabled(False)
        else:
            self.tab2.maxDistCombo.setEnabled(True)
            self.tab2.overlap.setEnabled(True)

    def output(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag, protToPepFlag,
               modList, maxMod, maxDistance, outputPath, chargeFlags, mgfFlag):

        """
        called by output pre-step function, it runs the generateOutput function from Mers.py; This is shown to be
        running via the progress bar
        """

        start = time.time()

        if maxDistance != 'None':
            maxDistance = int(maxDistance)

        self.fasta.generateOutput(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                                  protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, self.mgf, mgfFlag)
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

            # # wipe max distance values only if max is changed. Add None back to combo box.
            # self.tab2.maxDistCombo.clear()
            # self.tab2.maxDistCombo.addItem('None')
            # # refill max distance combo box with items greater than or equal to max
            # for i in range(int(text) - 1, 26):
            #     self.tab2.maxDistCombo.addItem(str(i + 1))
            # if maxDistInt >= int(text):
            #     indexDist = self.tab2.maxDistCombo.findText(str(maxDistValue))
            #     self.tab2.maxDistCombo.setCurrentIndex(indexDist)

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

        sender = self.tab2.sender()
        if sender.currentText() == 'Custom Modification':
            self.showCustomMod(sender)

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
        modChange[0].addItem("Custom Modification")
        modChange[1].addItem("Custom Modification")

        if modValue1 not in (text, 'None'):
            modChange[0].addItem(modValue1)
            indexMod1 = modChange[0].findText(modValue1)
            modChange[0].setCurrentIndex(indexMod1)

        if modValue2 not in (text, 'None'):
            modChange[1].addItem(modValue2)
            indexMod2 = modChange[1].findText(modValue2)
            modChange[1].setCurrentIndex(indexMod2)

    def showCustomMod(self, sender):

        sender.setCurrentIndex(0)
        self.formGroupBox = QGroupBox('Custom Modification')
        self.formLayout = QFormLayout()
        self.formLayout.addRow(QLabel("Input modified amino acids without spaces: TGN \n" +
                                      "Input mass change as a decimal number. Place a minus (-) \n" +
                                        "directly before to signify mass loss."))
        self.custAminoInput = QLineEdit()
        self.custMassInput = QLineEdit()
        self.addModButton = QPushButton("Create Modification")
        # partial method allows variable to be passed to connected function on click
        self.addCustToModListSender = partial(self.addCustToModlist, sender)
        self.addModButton.clicked.connect(self.addCustToModListSender)
        self.formLayout.addRow(QLabel("Modified Amino Acids: "), self.custAminoInput)
        self.formLayout.addRow(QLabel("Mass Change: "), self.custMassInput)
        self.formLayout.addRow(self.addModButton)
        self.formGroupBox.setLayout(self.formLayout)
        self.formGroupBox.show()

    def addCustToModlist(self, sender):
        aminoAcids = self.custAminoInput.text()
        massChange = self.custMassInput.text()
        modKey = "Custom " + aminoAcids
        modValue = []

        # validity checks
        for char in aminoAcids:
            if char not in monoAminoMass.keys():
                QMessageBox.about(self, "Message",
                                  'Amino Acids are not valid. Please leave no spaces and use capitals!')
                return
            else:
                modValue.append(char)


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

    def getInputParams(self):

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

        maxMod = self.tab2.maxModCombo.currentText()


        return ppmVal, toleranceLevel, mined, maxed, maxDistance, overlapFlag, transFlag, cisFlag, \
               linearFlag, csvFlag, pepToProtFlag, protToPepFlag, modList, maxMod, outputFlag, chargeFlags, minByIon,\
               byIonAccuracy, byIonFlag, mgfFlag

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
        for i in range(2, 26):
            self.tab2.maxDistCombo.addItem(str(i))

    def addFlagChecks(self):

        # initialise overlap, trans, cis and linear check boxes
        self.tab2.overlap = QCheckBox('Cis Overlap Off', self)
        self.tab2.trans = QCheckBox('Trans', self)
        self.tab2.trans.stateChanged.connect(self.disableMaxDist)  # connect trans check box to relevant function
        self.tab2.trans.stateChanged.connect(self.enableControl)
        self.tab2.cis = QCheckBox('Cis', self)
        self.tab2.cis.stateChanged.connect(self.enableControl)
        self.tab2.linear = QCheckBox('Linear', self)
        self.tab2.linear.stateChanged.connect(self.enableControl)

    def addChargeStates(self):

        """
        Add charge state check boxes
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

    def textBoxSender(self, sender):
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

    def textBoxChanged(self, input):
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

    def createTab1ParameterWidgets(self):
        self.pushButton1 = QPushButton("Select Fasta File")
        self.pushButton1.clicked.connect(self.uploadFasta)
        self.mgfButton = QPushButton("Select MGF File")
        self.mgfButton.clicked.connect(self.uploadMgfPreStep)
        self.mgfPlotFlag = QCheckBox('Produce Intensity Plot')
        self.nextTab = QPushButton("Next Tab")
        self.nextTab.clicked.connect(self.nextTabFunc)
        # self.nextTab.setEnabled(False)
        # self.tabs.setTabEnabled(1, False)
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

        # for i in range(10, 110, 10):
        #     self.tab1.ppmCombo.addItem(str(i))

        # intensities = [0, 10, 50, 100, 500, 1000, 5000, 10000]
        # for intensity in intensities:
        #     self.tab1.toleranceCombo.addItem(str(intensity))

        # for i in range(10, 100, 10):
        #     self.tab1.minByIonCombo.addItem(str(i))

        # ionAccuracies = [0.4, 0.2, 0.1, 0.05, 0.02, 0.01]
        # for accuracy in ionAccuracies:
        #     self.tab1.byIonAccCombo.addItem(str(accuracy))

    def addTab1ParameterWidgets(self):
        self.tab1.layout.setColumnStretch(0, 1)
        self.tab1.layout.setColumnStretch(5, 1)
        self.tab1.layout.setRowStretch(0, 1)
        self.tab1.layout.setRowStretch(5, 1)
        self.tab1.layout.addWidget(self.pushButton1, 1, 2)
        self.tab1.layout.addWidget(self.mgfButton, 2, 2)
        self.tab1.layout.addWidget(self.mgfFlag, 1, 3)
        self.tab1.layout.addWidget(self.mgfPlotFlag, 2, 3)
        self.tab1.layout.addWidget(self.tab1.ppmLabel, 3, 2)
        self.tab1.layout.addWidget(self.tab1.ppmText, 3, 3)
        self.tab1.layout.addWidget(self.tab1.ppmStatus, 3, 4)
        self.tab1.layout.addWidget(self.tab1.toleranceLabel, 4, 2)
        self.tab1.layout.addWidget(self.tab1.toleranceText, 4, 3)
        self.tab1.layout.addWidget(self.tab1.toleranceStatus, 4, 4)

        self.tab1.layout.addWidget(self.tab1.minByIonLabel, 5, 2)
        self.tab1.layout.addWidget(self.tab1.minByIonText, 5, 3)
        self.tab1.layout.addWidget(self.tab1.minByIonStatus, 5, 4)
        self.tab1.layout.addWidget(self.tab1.byIonAccLabel, 6, 2)
        self.tab1.layout.addWidget(self.tab1.byIonAccText, 6, 3)
        self.tab1.layout.addWidget(self.tab1.byIonAccStatus, 6, 4)
        self.tab1.layout.addWidget(self.tab1.byIonFlag, 7, 2)
        self.tab1.layout.addWidget(self.nextTab, 7, 3)

    def createTab2ParameterWidgets(self):

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
        self.tab2.layout.addWidget(self.tab2.maxMod, 7, 3)
        self.tab2.layout.addWidget(self.tab2.linear, 8, 3)
        self.tab2.layout.addWidget(self.tab2.cis, 9, 3)
        self.tab2.layout.addWidget(self.tab2.trans, 10, 3)
        self.tab2.layout.addWidget(self.tab2.overlap, 11, 3)

        self.tab2.layout.addWidget(self.tab2.csv, 8, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.pepToProt, 9, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.protToPep, 10, 4, 1, 3)

        self.tab2.layout.addWidget(self.tab2.chargeLabel, 13, 3)

        # all dynamic elements added to the grid layout of tab 2
        self.tab2.layout.addWidget(self.tab2.minimumCombo, 1, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximumCombo, 2, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistCombo, 3, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod1Combo, 4, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod2Combo, 5, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod3Combo, 6, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maxModCombo, 7, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.plusOne, 13, 4)
        self.tab2.layout.addWidget(self.tab2.plusTwo, 13, 5)
        self.tab2.layout.addWidget(self.tab2.plusThree, 13, 6)
        self.tab2.layout.addWidget(self.tab2.plusFour, 14, 4)
        self.tab2.layout.addWidget(self.tab2.plusFive, 14, 5)
        self.tab2.layout.addWidget(self.tab2.output, 15, 5, 1, 2)
        self.tab2.layout.addWidget(self.tab2.stop, 15, 3)

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
        self.tab2.trans.setEnabled(True)
        self.tab2.plusTwo.setChecked(True)

        # minByIonIndex = self.tab1.minByIonCombo.findText(self.minByIonDefault)
        # self.tab1.minByIonCombo.setCurrentIndex(minByIonIndex)

        # byIonAccIndex = self.tab1.byIonAccCombo.findText(self.byIonAccDefault)
        # self.tab1.byIonAccCombo.setCurrentIndex(byIonAccIndex)

        self.tab1.byIonFlag.setChecked(True)

    @pyqtSlot()
    def on_click(self):
        # print("\n")
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())


if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())




