import numpy
import pandas as pd
from pyteomics import mgf
import math
from bisect import bisect_left
from MonoAminoAndMods import *


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

# **redid
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

            # ion dict -> {'b/yion: mass'}

            #print(byIonArray)

            for charge, chargeMass in value[2].items():
                # Shift to outside for charge for loop
                if alphaKey not in matchedPeptides:

                    closest = takeClosest(mgfObj.mgfDf[charge], chargeMass)
                    if pepMatch(chargeMass, closest, mgfObj.ppmVal):

                        #mzArray = mgfObj.pepmassIonArray[(charge, closest)]
                        #simIons = findSimIons(mzArray, byIonArray, mgfObj.byIonAccuracy)

                        if True == False:
                            matchedPeptides.add(alphaKey)
                        else:
                            byIonDict = initIonMass(key)
                            byIonArray = sortBYDict(byIonDict)
                            mzArray = mgfObj.pepmassIonArray[(charge, closest)]
                            simIons = findSimIons(mzArray, byIonArray, 0.1)
                            if simIons >= 2:
                                matchedPeptides.add(alphaKey)

                        # count similar ions and add that to simComparisons. Note that mzArray have multiple lists
                        #matchedPeptides.add(alphaKey)
                else:
                    break
            # check it passes max simcomparisons and then add alphakey to matchedpeptides!

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

# ** redid
def readMGF(input_path):
    """
    Creates a pandas dataframe based on mgf data
    """
    uniqueSpec = set()

    mgfDf = {}

    pepmassIonArray = {}

    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:

            if 'charge' in spectrum['params'].keys():
                charge = spectrum['params']['charge'][0]
                pepmass = spectrum['params']['pepmass'][0]

                maxIntensity = max(spectrum['intensity array'])

                chargePepmassTup = (charge, pepmass)

                mzArray = spectrum['m/z array']

                # bring pepmassIonArray out here so that duplicates aren't ignored. If the the chargepepmasstuple
                # already exists, append the mzArray to the existing one as done below.

                # Add it to the dataframe if they are not already in the set
                if chargePepmassTup not in uniqueSpec and maxIntensity > 0:

                    if charge in mgfDf:
                        mgfDf[charge].append(pepmass)
                        #pepmassIonArray[(charge,pepmass)] = mzArray
                    else:
                        mgfDf[charge] = [pepmass]
                        #pepmassIonArray[(charge,pepmass)] = mzArray

                if chargePepmassTup in pepmassIonArray:
                    pepmassIonArray[chargePepmassTup].append(mzArray)
                    #print(pepmassIonArray[chargePepmassTup])
                else:
                    pepmassIonArray[chargePepmassTup] = [mzArray]

                    # mgfDf.loc[len(mgfDf)] = [spectrum['params']['charge'][0],
                    #                          spectrum['params']['pepmass'][0]]

                uniqueSpec.add(chargePepmassTup)

    sortMgfDFValues(mgfDf)
    sortPepmassIonArray(pepmassIonArray)
    return mgfDf, pepmassIonArray

# ** changed to two functions
def sortMgfDFValues(mgfDf):
    for charge, masses in mgfDf.items():
        masses.sort()

def sortPepmassIonArray(pepmassIonArray):
    for chargeMassTup, masses in pepmassIonArray.items():
        for list in masses:
            list.sort()

#readMGF('C:/Users/Arpit/Desktop/UROP/InputData/600MB.mgf')
#mgfObj = MGF(readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf'))
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

# ** added all function beyond here
def createBYIons(peptide):
    blist = []
    ylist = []
    for i in range (0,len(peptide)-1):
        b = peptide[0:i+1]
        y = peptide[i+1:]
        blist.append(b)
        ylist.append(y)
    return blist, ylist


def bMassCalc(peptide):
    mass = 1
    for char in peptide:
        mass += monoAminoMass[char]
    return mass

def yMassCalc(peptide):
    mass = H20_MASS + 1
    for char in peptide:
        mass += monoAminoMass[char]
    return mass

def ionMassDict(blist,ylist):
    dict = {}
    for i in range(0,len(blist)):
        pepB = blist[i]
        pepY = ylist[i]
        dict[pepB] = bMassCalc(pepB)
        dict[pepY] = yMassCalc(pepY)
    return dict

"""returns a dictionary holding the b and y ions and their correspondigg mass"""
def initIonMass(peptide):
    blist, ylist = createBYIons(peptide)
    dict = ionMassDict(blist, ylist)
    return dict

def sortBYDict(byIonDict):
    byIonArray = []
    for key, value in byIonDict.items():
        byIonArray.append(value)
    byIonArray.sort()
    return byIonArray

def findSimIons(mzArray, byIons, accuracy):
    simIonsArray = []
    for array in mzArray:
        simTemp = 0
        for mass in byIons:
            closest = takeClosest(array, mass)
            upperThresh = mass + accuracy
            lowerThresh = mass - accuracy
            if lowerThresh < closest < upperThresh:
                simTemp += 1
            simIonsArray.append(simTemp)
        simIons = max(simIonsArray)
    return simIons