import numpy
import pandas as pd
from pyteomics import mgf
import math
from bisect import bisect_left

# 347805 entries - 6 minutes
class MGF:

    def __init__(self, mgfDf):
        self.mgfDf = mgfDf
        self.mgfEntries = len(mgfDf)
        self.ppmVal = None
        self.toleranceLevel = None
        self.tempMgfDf = None

    def initValues(self, ppmVal, toleranceLevel):
        self.ppmVal = ppmVal
        self.toleranceLevel = toleranceLevel

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

        self.tempMgfDf.sort_values(by=['CHARGE_STATE', 'PEPMASS'], inplace = True)
        groupedTempDf = self.tempMgfDf.groupby('CHARGE_STATE')
        self.tempMgfDf = groupedTempDf['PEPMASS'].unique()
        self.tempMgfDf = self.tempMgfDf.to_frame()
        self.tempMgfDf.reset_index(inplace=True)


def generateMGFList(mgfObj, massDict):

    if len(mgfObj.tempMgfDf.index) != 0:
        matchedPeptides = set()
        for key, value in massDict.items():

            for charge, chargeMass in value[2].items():

                chargeList = mgfObj.tempMgfDf.loc[mgfObj.tempMgfDf['CHARGE_STATE'] == charge, 'PEPMASS'].iloc[0]
                closest = takeClosest(chargeList, chargeMass)
                if pepMatch(chargeMass, chargeList[closest], mgfObj.ppmVal):
                    matchedPeptides.add(key)


                # if pepMatch(chargeMass, chargeList[closest], mgfObj.ppmVal):
                #     diffPpm = calcPpm(chargeMass, chargeList[closest])
                #
                #     print(key, charge, chargeMass, chargeList[closest], diffPpm)
        return matchedPeptides

def pepMatch(predictedMass, pepmass, ppmVal):

    currentPpm = calcPpm(predictedMass, pepmass)
    if int(round(currentPpm)) <= ppmVal:
        return True
    return False


def calcPpm(predictedMass, pepmass):
    a = (abs(predictedMass - pepmass) / predictedMass)*1000000
    return a


def readMGF(input_path):
    """
    Creates a pandas dataframe based on mgf data
    """
    uniqueSpec = set()
    colNames = ['CHARGE_STATE', 'PEPMASS']
    mgfDf = pd.DataFrame(columns=colNames)
    with mgf.read(input_path) as mgfReader:
        totalEntries = 0
        chargedEntries = 0
        for spectrum in mgfReader:

            if 'charge' in spectrum['params'].keys():
                chargePepmassTup = (spectrum['params']['charge'][0], spectrum['params']['pepmass'][0])
                # Add it to the dataframe if they are not already in the set
                if chargePepmassTup not in uniqueSpec:
                    mgfDf.loc[len(mgfDf)] = [spectrum['params']['charge'][0],
                                             spectrum['params']['pepmass'][0]]
                uniqueSpec.add(chargePepmassTup)
                chargedEntries+=1

            totalEntries+=1


    # print("There are " + str(totalEntries) + " total entries in this MGF file")
    # print("There are " + str(chargedEntries) + " entries with a charge in this MGF file")
    # print("There are " + str(len(uniqueSpec)) + " unique entries in this MGF file")
    # print("There are " + str(chargedEntries - len(uniqueSpec)) + " duplicate entries in this MGF file")

    return mgfDf

#readMGF('C:/Users/Arpit/Desktop/UROP/InputData/600MB.mgf')
#readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf')
#readMGF('C:/Users/Administrator/Desktop/UROP/InputData/918MB.mgf')
# readMGF('C:/Users/Administrator/Desktop/UROP/InputData/MgfExample.mgf')


def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber via index.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return -1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos-1

# actualMass = 495.25851750000004
# pepmass = 495.7115
# ppmVal = 90

#
# chargeFlags = [False, True, False, True, False]
# mgf = MGF(readMGF("MgfExample.mgf"))
# mgf.removeChargeStates(chargeFlags)
# print(mgf.mgfDf.head())

#
#
# pepList = [559.28780, 559.2831, 560.2133, 600, 231.23]
# pepList.sort()
# print(pepList)
# actualMass = 894
# print(takeClosest(pepList, actualMass))

