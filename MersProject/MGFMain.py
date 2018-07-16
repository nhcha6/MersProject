import numpy
import pandas as pd
from pyteomics import mgf
import math
from bisect import bisect_left


class MGF:

    def __init__(self, mgfDf):
        self.mgfDf = mgfDf
        self.mgfEntries = len(mgfDf)
        self.tempMgfDf = None

    def removeChargeStates(self, chargeFlags):
        """
        Remove all charges that are irrelevant, which is given by the chargeFlags params
        """
        firstTrue = True
        for i in range(0, len(chargeFlags)):

            if not chargeFlags[i]:
                # Comparing to i+1 because of charge state!

                # Resolve the bug where would have to reattach the mgf file everytime charge state is changed.
                if firstTrue:
                    self.tempMgfDf = self.mgfDf.drop(self.mgfDf[self.mgfDf.CHARGE_STATE == i+1].index)
                    firstTrue = False
                else:
                    self.tempMgfDf.drop(self.mgfDf[self.mgfDf.CHARGE_STATE == i + 1].index, inplace=True)

        pepmassList = self.tempMgfDf['PEPMASS'].tolist()
        pepmassList.sort()
        return pepmassList

def readMGF(input_path):
    """
    Creates a pandas dataframe based on mgf data
    """

    colNames = ['CHARGE_STATE', 'PEPMASS']
    mgfDf = pd.DataFrame(columns=colNames)
    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:

            mgfDf.loc[len(mgfDf)] = [spectrum['params']['charge'][0],
                                     spectrum['params']['pepmass'][0]]
    return mgfDf


def ppmToPercent(ppmVal):
    return ppmVal/10000


def ppmPosNeg(actualMass, ppmVal):
    ppmPercent = ppmToPercent(ppmVal)
    ppmPositive = actualMass + actualMass*ppmPercent
    ppmNegative = actualMass - actualMass*ppmPercent
    print("PPM Positive value is: " + str(ppmPositive))
    print("PPM Negative value is: " + str(ppmNegative))
    return ppmPositive, ppmNegative


def ppmCheck(pepmass, ppmPositive, ppmNegative, tolerance):
    negCheck = math.isclose(pepmass, ppmPositive, abs_tol = tolerance)
    posCheck = math.isclose(pepmass, ppmNegative, abs_tol = tolerance)
    return negCheck, posCheck


def addMass(listType, pepmass, actualMass, ppmVal):
    ppmPositive, ppmNegative = ppmPosNeg(actualMass, ppmVal)
    negCheck, posCheck = ppmCheck(ppmPositive, ppmNegative)
    if negCheck:
        listType.append(negCheck)
    if posCheck:
        listType.append(posCheck)


def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

actualMass = 495.25851750000004
pepmass = 495.7115
ppmVal = 90
tolerance = 0.000001
# ppmPositive, ppmNegative = ppmPosNeg(actualMass, ppmVal)
# ppmPosCheck, ppmNegCheck = ppmCheck(pepmass, ppmPositive, ppmNegative, tolerance)
#
# print('Actual mass: ' + str(actualMass))
#
# print('PPM PosCheck: ' + str(ppmPosCheck))
# print('PPM NegCheck: ' + str(ppmNegCheck))
#
# chargeFlags = [False, True, False, True, False]
# mgf = MGF(readMGF("MgfExample.mgf"))
# mgf.removeChargeStates(chargeFlags)
# print(mgf.mgfDf.head())



pepList = [559.28780, 559.2831, 560.2133, 600, 231.23]
pepList.sort()
actualMass = 100
print(takeClosest(pepList, actualMass))