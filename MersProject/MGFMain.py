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

    if mgfObj.mgfDf:

        matchedPeptides = set()
        for key, value in massDict.items():

            # convert modified peptides to original form
            if not key.isalpha():
                alphaKey = modToPeptide(key)
            else:
                alphaKey = key

            for charge, chargeMass in value[2].items():
                if alphaKey not in matchedPeptides:

                    # chargeList = mgfObj.mgfDf[charge]
                    closest = takeClosest(mgfObj.mgfDf[charge], chargeMass)
                    if pepMatch(chargeMass, closest, mgfObj.ppmVal):

                        matchedPeptides.add(alphaKey)
                else:
                    break


                # if pepMatch(chargeMass, chargeList[closest], mgfObj.ppmVal):
                #     diffPpm = calcPpm(chargeMass, chargeList[closest])
                #
                #     print(key, charge, chargeMass, chargeList[closest], diffPpm)
        return matchedPeptides


def modToPeptide(moddedPeptide):
    peptide = ''.join(filter(lambda x: x.isalpha(), moddedPeptide))

    return peptide.upper()


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

    mgfDf = {}

    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:

            if 'charge' in spectrum['params'].keys():
                charge = spectrum['params']['charge'][0]
                pepmass = spectrum['params']['pepmass'][0]

                maxIntensity = max(spectrum['intensity array'])

                chargePepmassTup = (charge, pepmass)

                # Add it to the dataframe if they are not already in the set
                if chargePepmassTup not in uniqueSpec and maxIntensity > 5000:

                    if charge in mgfDf:
                        mgfDf[charge].append(pepmass)
                    else:
                        mgfDf[charge] = [pepmass]


                uniqueSpec.add(chargePepmassTup)

    sortMgfDf(mgfDf)
    return mgfDf

def sortMgfDf(mgfDf):
    for charge, masses in mgfDf.items():
        masses.sort()

#readMGF('C:/Users/Arpit/Desktop/UROP/InputData/600MB.mgf')
mgfObj = MGF(readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf'))
# print((1,mgfObj.mgfDf[1]))
# print(2,(mgfObj.mgfDf[2]))
# print((3,mgfObj.mgfDf[3]))

#readMGF('C:/Users/Administrator/Desktop/UROP/InputData/918MB.mgf')
# readMGF('C:/Users/Administrator/Desktop/UROP/InputData/MgfExample.mgf')


def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber via index.

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

