### Need to fix no parameter
import sys
import subprocess
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
                            QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
                            QProgressBar
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
from Mers import *

# pyinstaller MersGUI --> this command from the relevant file location creates executable file

class WorkerSignals(QObject):
    finished = pyqtSignal()
    updateProgBar = pyqtSignal(int)
    generateProgBar = pyqtSignal()

class ProgressGenerator(QRunnable):

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

    def __init__(self, fn, *args, **kwargs):
        super(OutputGenerator, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        #print(*self.args)
        #print(**self.kwargs)
        self.fn(*self.args)
        self.signals.finished.emit()

class App(QMainWindow):
    # App serves as the parent class for the embedded MyTableWidget

    # Initialisation of main window class
    def __init__(self):
        super().__init__()
        self.title = 'Peptide Splicer'
        self.left = 500
        self.fastaTest = False
        self.outputPath = ""
        self.statusbar = self.statusBar()

        self.center()

        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)

        self.show()

    # center function is called to centre the main window on the screen
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


# MyTableWidget class is the child of the App class. It holds the tabs where most of the GUI functionality occurs
class MyTableWidget(QWidget):

    # Initialisation of table child
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        self.fasta = None
        #Init threading
        self.threadpool = QThreadPool()
        self.minDefault = '8'
        self.maxDefault = '12'
        self.maxDistDefault = '25'

        # Initialisation of two tabs
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tabs.resize(300, 200)

        # Add tabs to table class (self)
        self.tabs.addTab(self.tab1, "Select File and Path")
        self.tabs.addTab(self.tab2, "Input Parameters")

        # Creation of tab layout and widgets within tab
        self.tab1.layout = QVBoxLayout(self)

        self.pushButton1 = QPushButton("Select File")
        self.pushButton1.clicked.connect(self.uploadFasta)

        self.tab1.layout.addWidget(self.pushButton1)

        self.tab1.setLayout(self.tab1.layout)

        # Create second tab layout and widgets within each tab
        self.tab2.layout = QGridLayout(self)
        self.tab2.layout.setSpacing(10)

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



        # set default values
        maxDistIndex = self.tab2.maxDistCombo.findText(str(self.maxDistDefault))
        self.tab2.maxDistCombo.setCurrentIndex(maxDistIndex)

        maxIndex = self.tab2.maximumCombo.findText(str(self.maxDefault))
        self.tab2.maximumCombo.setCurrentIndex(maxIndex)

        minIndex = self.tab2.minimumCombo.findText(str(self.minDefault))
        self.tab2.minimumCombo.setCurrentIndex(minIndex)


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
        for modification in modificationList:
            self.tab2.mod1Combo.addItem(modification)
            self.tab2.mod2Combo.addItem(modification)
            self.tab2.mod3Combo.addItem(modification)

        # initialise overlap, trans, cis and linear check boxes
        self.tab2.overlap = QCheckBox('Overlap Off', self)
        self.tab2.trans = QCheckBox('Trans', self)
        self.tab2.trans.stateChanged.connect(self.disableMaxDist) # connect trans check box to relevant function
        self.tab2.cis = QCheckBox('Cis', self)
        self.tab2.linear = QCheckBox('Linear', self)

        # set to true as defauls for linear, cis and overlap off. Set trans off for now.
        self.tab2.overlap.setChecked(True)
        self.tab2.cis.setChecked(True)
        self.tab2.linear.setChecked(True)
        self.tab2.trans.setEnabled(False)



        # create generate output push button
        self.tab2.output = QPushButton('Generate Output!', self)
        self.tab2.output.clicked.connect(self.confirmationFunction)

        # Add charge state check boxes
        self.tab2.chargeLabel = QLabel('Charge states (z): ')
        self.tab2.plusOne = QCheckBox('+1', self)
        self.tab2.plusTwo = QCheckBox('+2', self)
        self.tab2.plusThree = QCheckBox('+3', self)
        self.tab2.plusFour = QCheckBox('+4', self)
        self.tab2.plusFive = QCheckBox('+5', self)
        self.tab2.plusTwo.setChecked(True)

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
        self.tab2.layout.addWidget(self.tab2.chargeLabel, 11, 3)

        # all dynamic elements added to the grid layout of tab 2
        self.tab2.layout.addWidget(self.tab2.minimumCombo, 1, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximumCombo, 2, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistCombo, 3, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod1Combo, 4, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod2Combo, 5, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.mod3Combo, 6, 4, 1, 3)
        self.tab2.layout.addWidget(self.tab2.plusOne, 11, 4)
        self.tab2.layout.addWidget(self.tab2.plusTwo, 11, 5)
        self.tab2.layout.addWidget(self.tab2.plusThree, 11, 6)
        self.tab2.layout.addWidget(self.tab2.plusFour, 12, 4)
        self.tab2.layout.addWidget(self.tab2.plusFive, 12, 5)
        self.tab2.layout.addWidget(self.tab2.output, 13, 5, 1, 2)

        # set layout
        self.tab2.setLayout(self.tab2.layout)

        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    # called from the Upload Fasta File button. Opens a window to select a file, and check if the file ends in fasta
    def uploadFasta(self):
        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home')

        self.fastaTest = fname[0][-5:]
        if self.fastaTest == 'fasta':
            self.fasta = Fasta(addSequenceList(fname[0]))

            QMessageBox.about(self, "Message", 'Fasta file successfully imported!')
        elif fname[0] == '':
            print('')
        else:
            QMessageBox.about(self, "Message", 'Please select a Fasta file!')

    # called from the Select Output Path button. Opens a window to select a file location to save the output to.
    def getOutputPath(self):

        self.outputPath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        if self.outputPath == '':
            return False
        print(self.outputPath)
        return True
        # else:
            # convert to tool tip later
            # QMessageBox.about(self, "Message", 'Valid Path Selected')



    def confirmationFunction(self):

        """
        called on click of generate output button on tab2. Checks to ensure all input values are relevant and outputs
        message box summarising the inputs of the user. When yes is clicked on the message box, the output function is
        called which begins generating results
        """

        mined = int(self.tab2.minimumCombo.currentText())
        maxed = int(self.tab2.maximumCombo.currentText())
        overlapFlag = self.tab2.overlap.isChecked()
        transFlag = self.tab2.trans.isChecked()

        cisFlag = self.tab2.cis.isChecked()
        maxDistance = self.tab2.maxDistCombo.currentText()
        linearFlag = self.tab2.linear.isChecked()

        outputFlag = cisFlag or linearFlag or transFlag

        plusOneFlag = self.tab2.plusOne.isChecked()
        plusTwoFlag = self.tab2.plusTwo.isChecked()
        plusThreeFlag = self.tab2.plusThree.isChecked()
        plusFourFlag = self.tab2.plusFour.isChecked()
        plusFiveFlag = self.tab2.plusFive.isChecked()

        chargeFlags = [plusOneFlag, plusTwoFlag, plusThreeFlag, plusFourFlag, plusFiveFlag]

        # self.fasta = Fasta(addSequenceList('/Users/nicolaschapman/Documents/UROP/Code/MersProject/small.fasta'))
        # self.outputPath = '/Users/nicolaschapman/Desktop/Mers Output'
        # self.fasta = Fasta(addSequenceList('C:/Users/Arpit/Desktop/UROP/Example.fasta'))
        # self.outputPath = 'C:/Users/Arpit/Desktop/UROP'
        modList = [self.tab2.mod1Combo.currentText(), self.tab2.mod2Combo.currentText(),
                   self.tab2.mod3Combo.currentText()]

        if self.fasta is None:

            QMessageBox.about(self, "Message", 'Please check that a valid Fasta file and output '
                                               'file location have been selected')
        elif not outputFlag:
            QMessageBox.about(self, "Message", 'Please select at least one output type; either trans, cis or linear')
        else:
            reply = QMessageBox.question(self, 'Message', 'Do you wish to confirm the following input?\n' +
                                         'Minimum Length: ' + str(mined) + '\n' +
                                         'Maximum Length: ' + str(maxed) + '\n' +
                                         'Overlap Flag: ' + str(overlapFlag) + '\n' +
                                         'Trans Flag: ' + str(transFlag) + '\n' +
                                         'Linear Flag: ' + str(linearFlag) + '\n' +
                                         'Cis Flag: ' + str(cisFlag) + '\n' +
                                         'Mod List: ' + str(modList) + '\n' +
                                         'Maximum Distance: ' + str(maxDistance) + '\n',
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.Yes:
                print('reply is yes')

                if self.getOutputPath():
                    self.outputPreStep(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList, maxDistance,
                      self.outputPath, chargeFlags)
                    #self.output(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList,
                     #           maxDistance, self.outputPath, chargeFlags)

    def finished(self):
        print("ITS DONE")
        self.progressBarUpdate.changeFlag()
        #self.deleteProgressBar()
        QMessageBox.about(self, "Message", 'Output Complete')

    def updateProgressBar(self,int):
        self.tab2.progressBar.setValue(int)

    def deleteProgressBar(self):
        # Delete progress label and progress bar
        self.tab2.layout.removeWidget(self.tab2.progressLabel)
        self.tab2.progressLabel.deleteLater()
        self.tab2.progressLabel = None
        self.tab2.layout.removeWidget(self.tab2.progressBar)
        self.tab2.progressBar.deleteLater()
        self.tab2.progressBar = None

    def outputPreStep(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList, maxDistance,
                      outputPath, chargeFlags):

        self.tab2.progressLabel = QLabel('Collating Combinations. Please Wait: ')
        self.tab2.layout.addWidget(self.tab2.progressLabel, 14, 3, 1, 2)
        self.tab2.progressBar = QProgressBar(self)
        self.tab2.layout.addWidget(self.tab2.progressBar, 15, 3, 1, 4)

        outputGen = OutputGenerator(self.output, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList,
                               maxDistance, outputPath, chargeFlags)

        outputGen.signals.finished.connect(self.finished)
        self.threadpool.start(outputGen)

        self.progressBarUpdate = ProgressGenerator()
        self.progressBarUpdate.signals.updateProgBar.connect(self.updateProgressBar)
        self.progressBarUpdate.signals.finished.connect(self.deleteProgressBar)
        self.threadpool.start(self.progressBarUpdate)



    # called when trans is selected, it disables the use of the max distance function
    def disableMaxDist(self, state):

        if state == Qt.Checked:
            index = self.tab2.maxDistCombo.findText('None')
            self.tab2.maxDistCombo.setCurrentIndex(index)
            self.tab2.maxDistCombo.setEnabled(False)
        else:
            self.tab2.maxDistCombo.setEnabled(True)

    # called by confirmation function, it runs the generateOutput function from Mers.py while outputing small
    # bits of information to the user via a statusbar in the GUI
    def output(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList,
               maxDistance, outputPath, chargeFlags):
        start = time.time()

        #self.parent().statusbar.showMessage('Processing Data')

        if maxDistance != 'None':
            maxDistance = int(maxDistance)

        self.fasta.generateOutput(mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList,
                                  maxDistance, outputPath, chargeFlags)
        end = time.time()
        #self.parent().statusbar.hide()
        print(end - start)


        #print(tryString)
        #replacedOutpath = outputPath.replace("/", '\\')
        #print(replacedOutpath)
        #openString = r'explorer "' + replacedOutpath + '"'
        #print(openString)
        # subprocess.Popen(openString)


    # called when minimumCombo value changes. It alters the values available in max and maxDistance combos to
    # ensure a realistic input
    def minMaxChanged(self, text):

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
            for i in range(int(text)-1, 26):
                comboChange.addItem(str(i+1))
            # Restores current value if it is greater than the min
            if value >= int(text):
                indexMax = comboChange.findText(str(value))
                comboChange.setCurrentIndex(indexMax)

        else: #maxChanged
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
        modCombos = [self.tab2.mod1Combo, self.tab2.mod2Combo, self.tab2.mod3Combo]
        modSender = []
        modChange = []

        for combo in modCombos:
            modSender.append(combo == self.tab2.sender())

        for i in range (0,len(modSender)):
            if not modSender[i]:
                modChange.append(modCombos[i])

        modValue1 = modChange[0].currentText()
        modValue2 = modChange[1].currentText()

        modChange[0].clear()
        modChange[1].clear()
        modChange[0].addItem('None')
        modChange[1].addItem('None')

        for modification in modificationList:
            if modification not in (text, modValue1, modValue2): #and modification != modValue1 and modification != modValue2:
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



    @pyqtSlot()
    def on_click(self):
        print("\n")
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
