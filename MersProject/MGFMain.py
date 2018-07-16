import numpy
import pandas as pd
from pyteomics import mgf
import math


class MGF:

    def __init__(self, mgfDf):
        self.mgfDf = mgfDf
        self.mgfEntries = len(mgfDf)

    def removeChargeStates(self, chargeFlags):
        """
        Remove all charges that are irrelevant, which is given by the chargeFlags params
        """
        for i in range(0, len(chargeFlags)):
            if not chargeFlags[i]:
                # Comparing to i+1 because of charge state!
                self.mgfDf.drop(self.mgfDf[self.mgfDf.CHARGE_STATE == i+1].index, inplace=True)

        print(self.mgfDf.head())


def readMGF(input_path):
    """
    Creates a pandas dataframe based on mgf data
    """

    colNames = ['TITLE', 'CHARGE_STATE', 'PEPMASS']
    mgfDf = pd.DataFrame(columns=colNames)
    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:

            mgfDf.loc[len(mgfDf)] = [spectrum['params']['title'], spectrum['params']['charge'][0],
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
