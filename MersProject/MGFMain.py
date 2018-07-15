import numpy
import pandas as pd
from pyteomics import mgf
import math


def readMGF(input_path):
    colNames = ['TITLE', 'PEPMASS']
    mgfDf = pd.DataFrame(columns=colNames)
    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:
            mgfDf.loc[len(mgfDf)] = [spectrum['params']['title'], spectrum['params']['pepmass'][0]]
    return mgfDf

def ppmToPercent(ppmVal):
    return ppmVal/10000

def ppmPosNeg(actualMass, ppmVal):
    ppmPercent = ppmToPercent(ppmVal)
    print(ppmPercent)
    ppmPositive = actualMass + actualMass*ppmPercent
    ppmNegative = actualMass - actualMass*ppmPercent
    return ppmPositive, ppmNegative

def ppmCheck(pepmass, ppmPositive, ppmNegative):
    negCheck = math.isclose(pepmass, ppmNegative, abs_tol = 0.00001)
    posCheck = math.isclose(pepmass, ppmNegative, abs_tol = 0.00001)
    return negCheck, posCheck

def addMass(listType, pepmass, actualMass, ppmVal):
    ppmPositive, ppmNegative = ppmPosNeg(actualMass, ppmVal)
    negCheck, posCheck = ppmCheck(ppmPositive, ppmNegative)
    if negCheck:
        listType.append(negCheck)
    if posCheck:
        listType.append(posCheck)

actualMass = 495.25851750000004
ppmVal = 10
ppmPositive, ppmNegative = ppmCheck(1, actualMass, ppmVal)

print('Actual mass: ' + str(actualMass))
print('PPM Positive:' + str(ppmPositive))
print('PPM Negative:' + str(ppmNegative))
#
# mgfDf = readMGF("MgfExample.mgf")
# mgfDf.set_index('TITLE', inplace=True, drop=True)
#
# print(mgfDf.head(5))