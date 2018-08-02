import numpy
import pandas as pd
from pyteomics import mgf
import math
from bisect import bisect_left


class MGF:

    """
    Class to represent MGF input data
    """

    def __init__(self, mgfDf, pepmassIonArray):
        # mgfDf looks like: {'charge': [list of masses]}
        self.mgfDf = mgfDf
        self.pepmassIonArray = pepmassIonArray
        self.mgfEntries = len(mgfDf)
        self.ppmVal = None
        self.toleranceLevel = None

    def initValues(self, ppmVal, toleranceLevel):

        """
        Add extra info required such as the ppmValue!
        """
        self.ppmVal = ppmVal
        self.toleranceLevel = toleranceLevel


def generateMGFList(mgfObj, massDict):
    """

    Generates the list of unique peptides that have masses that match within the specified
    """
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
    colNames = ['CHARGE_STATE', 'PEPMASS']
    mgfDf = {}
    pepmassIonArray = {}
    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:

            if 'charge' in spectrum['params'].keys():
                charge = spectrum['params']['charge'][0]
                pepmass = spectrum['params']['pepmass'][0]
                chargePepmassTup = (charge, pepmass)

                mzArray = spectrum['m/z array']

                # Add it to the dataframe if they are not already in the set
                if chargePepmassTup not in uniqueSpec:

                    if charge in mgfDf:
                        mgfDf[charge].append(pepmass)
                        pepmassIonArray[(charge,pepmass)] = mzArray
                    else:
                        mgfDf[charge] = [pepmass]
                        pepmassIonArray[(charge,pepmass)] = mzArray
                    # mgfDf.loc[len(mgfDf)] = [spectrum['params']['charge'][0],
                    #                          spectrum['params']['pepmass'][0]]



                uniqueSpec.add(chargePepmassTup)
                break

    sortDictValues(mgfDf)
    sortDictValues(pepmassIonArray)
    return mgfDf, pepmassIonArray


def sortDictValues(mgfDf):
    for charge, masses in mgfDf.items():
        masses.sort()
# print(readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf'))
# mgfDf, pepmassIonArray = readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf')
# mgfObj = MGF(readMGF(mgfDf, pepmassIonArray))
# print(2,(mgfObj.mgfDf[2][0][3]))
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

