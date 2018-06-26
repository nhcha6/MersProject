import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout, \
                            QLineEdit, QTextEdit, qApp, QFileDialog, QAction, QGridLayout, QLabel, QComboBox, QCheckBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from Mers import *


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Peptide Splicer'
        self.left = 500
        self.top = 300
        self.width = 300
        self.height = 300
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)



        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)

        self.show()


class MyTableWidget(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        self.fasta = None

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()

        self.tabs.resize(300, 200)

        # Add tabs
        self.tabs.addTab(self.tab1, "Tab 1")
        self.tabs.addTab(self.tab2, "Tab 2")


        # Create first tab
        self.tab1.layout = QVBoxLayout(self)

        self.pushButton1 = QPushButton("Select File")
        self.pushButton1.clicked.connect(self.showDialog)

        self.pushButton2 = QPushButton("PRINT PATH")
        self.pushButton2.clicked.connect(self.createFasta)

        self.tab1.layout.addWidget(self.pushButton1)
        self.tab1.layout.addWidget(self.pushButton2)
        self.tab1.textEdit1 = QTextEdit()
        self.tab1.layout.addWidget(self.tab1.textEdit1)
        self.tab1.setLayout(self.tab1.layout)

        # Create second tab
        self.tab2.layout = QGridLayout(self)
        self.tab2.layout.setSpacing(10)

        # Minimum/maximum values
        self.tab2.minimum = QLabel('Minimum Value : ')
        self.tab2.minimumCombo = QComboBox(self)
        self.tab2.maximum = QLabel('Maximum Value : ')
        self.tab2.maximumCombo = QComboBox(self)

        self.tab2.maxDistance = QLabel('Maximum Distance : ')
        self.tab2.maxDistCombo = QComboBox(self)

        self.tab2.mod1 = QLabel('Modification 1 : ')
        self.tab2.mod2 = QLabel('Modification 2 : ')
        self.tab2.mod3 = QLabel('Modification 3 : ')

        self.tab2.overlap = QCheckBox('Overlap: ', self)
        self.tab2.cistrans = QCheckBox('Combine All ', self)

        self.tab2.output = QPushButton('Generate Output!', self)
        self.tab2.output.clicked.connect(self.printValues)


        for i in range(0, 26):
            self.tab2.minimumCombo.addItem(str(i))
            self.tab2.maximumCombo.addItem(str(i))
            self.tab2.maxDistCombo.addItem(str(i))

        # All the labels added
        self.tab2.layout.addWidget(self.tab2.minimum, 1, 3)
        self.tab2.layout.addWidget(self.tab2.maximum, 2, 3)
        self.tab2.layout.addWidget(self.tab2.maxDistance,3, 3)
        self.tab2.layout.addWidget(self.tab2.mod1, 4, 3)
        self.tab2.layout.addWidget(self.tab2.mod2, 5, 3)
        self.tab2.layout.addWidget(self.tab2.mod3, 6, 3)
        self.tab2.layout.addWidget(self.tab2.overlap, 7, 3)
        self.tab2.layout.addWidget(self.tab2.cistrans, 8, 3)

        self.tab2.layout.addWidget(self.tab2.minimumCombo, 1, 4)
        self.tab2.layout.addWidget(self.tab2.maximumCombo, 2, 4)
        self.tab2.layout.addWidget(self.tab2.maxDistCombo, 3, 4)
        self.tab2.layout.addWidget(self.tab2.output, 9, 4)

        self.tab2.setLayout(self.tab2.layout)



        #openFile = QAction(QIcon('open.png'), 'Open', self)
        #openFile.setShortcut('Ctrl+O')
        #openFile.setStatusTip('Open New File')
        #openFile.triggered.connect(self.showDialog)



        #self.setGeometry(300, 300, 350, 300)
        #self.setWindowTitle('File dialog')
        #self.show()



        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    def showDialog(self, parent):

        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home')
        self.fasta = Fasta(addSequenceList(fname[0]))

    def createFasta(self, parent):
        self.fasta.min = 3
        print(self.fasta.min)

    """
    def changeCombo(self):
        sender = self.sender()
        print(sender.text() + sender.currentText())"""

    #def genButton(self, parent, peptide, mined, maxed, overlapFlag, maxDistance):
     #   self.fasta.generateOutput(peptide, mined, maxed, overlapFlag, maxDistance)

    def printValues(self, parent):
        mined = int(self.tab2.minimumCombo.currentText())
        maxed = int(self.tab2.maximumCombo.currentText())
        maxDistance = int(self.tab2.maxDistCombo.currentText())
        overlapFlag = self.tab2.overlap.isChecked()
        combineAllFlag = self.tab2.cistrans.isChecked()
        peptide = "ABBCS"
        #fasta.generateOutput(peptide,mined,maxed,overlapFlag,maxDistance)





    @pyqtSlot()
    def on_click(self):
        print("\n")
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())