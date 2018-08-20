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


    def __init__(self, mgfDf, pepmassIonArray, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag):

        # mgfDf looks like: {'charge': [list of masses]}
        self.mgfDf = mgfDf
        self.pepmassIonArray = pepmassIonArray
        self.mgfEntries = len(mgfDf)

        self.ppmVal = ppmVal
        self.intensityThreshold = intensityThreshold
        self.minSimBy = minSimBy

        self.byIonAccuracy = byIonAccuracy
        self.byIonFlag = byIonFlag


def generateMGFList(mgfObj, massDict, modList):
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

                        if mgfObj.byIonFlag == False:
                            matchedPeptides.add(alphaKey)
                        else:
                            byIonDict = initIonMass(key, modList)
                            byIonArray = sortBYDict(byIonDict)
                            mzArray = mgfObj.pepmassIonArray[(charge, closest)]
                            simIons = findSimIons(mzArray, byIonArray, mgfObj.byIonAccuracy)
                            if simIons >= mgfObj.minSimBy:
                                print(simIons)
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



def readMGF(input_path, intensityThreshold):

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

                if chargePepmassTup not in uniqueSpec and maxIntensity > intensityThreshold:


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

    return mgfDf, pepmassIonArray


def sortMgfDFValues(mgfDf):
    for charge, masses in mgfDf.items():
        masses.sort()



#readMGF('C:/Users/Arpit/Desktop/UROP/InputData/600MB.mgf
#  print(readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf'))
# mgfDf, pepmassIonArray = readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf')
# mgfObj = MGF(readMGF(mgfDf, pepmassIonArray))
# print(2,(mgfObj.mgfDf[2][0][3]))

# readMGF('C:/Users/Arpit/Desktop/UROP/InputData/600MB.mgf')


def sortMgfDf(mgfDf):
    for charge, masses in mgfDf.items():
        masses.sort()

#readMGF('C:/Users/Arpit/Desktop/UROP/InputData/600MB.mgf')

# mgfObj = MGF(readMGF('C:/Users/Arpit/Desktop/UROP/InputData/MgfExample.mgf'))
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

def createBYIons(peptide):
    blist = []
    ylist = []
    for i in range (0,len(peptide)-1):
        b = peptide[0:i+1]
        y = peptide[i+1:]
        blist.append(b)
        ylist.append(y)
    return blist, ylist

def createBYIonsMod(peptide):
    blist = []
    ylist = []

    if peptide[-1].isalpha():
        endNo = 1
    else:
        endNo = 2

    for i in range(0,len(peptide)-endNo):
        if peptide[i].islower():
            b = peptide[0:i+2]
            y = peptide[i+2:]
            blist.append(b)
            ylist.append(y)
        elif not peptide[i].isalpha():
            continue
        else:
            b = peptide[0:i + 1]
            y = peptide[i + 1:]
            blist.append(b)
            ylist.append(y)
    return blist, ylist

def bMassCalc(peptide, modlist = None):
    mass = 1
    for char in peptide:
        if char.isalpha():
            char = char.upper()
            mass += monoAminoMass[char]
        else:
            mod = modlist[int(char)-1]
            mass += modTable[mod][-1]
    return mass

def yMassCalc(peptide, modlist = None):
    mass = H20_MASS + 1
    for char in peptide:
        if char.isalpha():
            char = char.upper()
            mass += monoAminoMass[char]
        else:
            mod = modlist[int(char)-1]
            mass += modTable[mod][-1]
    return mass

def ionMassDict(blist,ylist):
    dict = {}
    for i in range(0,len(blist)):
        pepB = blist[i]
        pepY = ylist[i]
        dict[pepB] = bMassCalc(pepB)
        dict[pepY] = yMassCalc(pepY)
    return dict

def ionMassDictMod(blist,ylist,modlist):
    dict = {}
    for i in range(0,len(blist)):
        pepB = blist[i]
        pepY = ylist[i]
        dict[pepB] = bMassCalc(pepB, modlist)
        dict[pepY] = yMassCalc(pepY, modlist)
    return dict

"""returns a dictionary holding the b and y ions and their correspondigg mass"""
def initIonMass(peptide, modList):
    if peptide.isalpha():
        blist, ylist = createBYIons(peptide)
        dict = ionMassDict(blist, ylist)
    else:
        blist, ylist = createBYIonsMod(peptide)
        dict = ionMassDictMod(blist, ylist, modList)
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
            simIonsArray.append(simTemp*100/len(byIons))
        simIons = max(simIonsArray)

    return simIons


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
