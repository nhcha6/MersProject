import numpy
import pandas as pd
from Mers import *
from pyteomics import mgf
import math
from bisect import bisect_left
from MonoAminoAndMods import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import math

class MGF:

    """
    Class to represent MGF input data
    """


    def __init__(self, mgfDfList, pepmassIonArrayList, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag, maxMass,
                 chargeMaxDict, mgfLen):

        # mgfDf looks like: {'charge': [list of masses]}
        self.mgfDfList = mgfDfList
        self.pepmassIonArrayList = pepmassIonArrayList
        self.mgfEntries = mgfLen
        self.maxMass = maxMass
        self.chargeMaxDict = chargeMaxDict

        self.ppmVal = ppmVal
        self.intensityThreshold = intensityThreshold
        self.minSimBy = minSimBy

        self.byIonAccuracy = byIonAccuracy
        self.byIonFlag = byIonFlag


def generateMGFList(protId, mgfObj, massDict, modList):
    """

    Generates the list of unique peptides that have masses that match within the specified.

    """
    if mgfObj.mgfDfList:

        # convert the list of dictionaries into dictionaries.
        mgfDf = mgfObj.mgfDfList[0]
        pepmassIonArray = mgfObj.pepmassIonArrayList[0]
        for i in range(1, len(mgfObj.mgfDfList)):
            # declare the portions you are about to combine
            mgfDfIter = mgfObj.mgfDfList[i]
            pepmassIonIter = mgfObj.pepmassIonArrayList[i]
            # update mgfDf
            for key, value in mgfDfIter.items():
                if key in mgfDf.keys():
                    mgfDf[key] = mgfDf[key] + value
                else:
                    mgfDf[key] = value
            # update pepmassIonArray
            for key, value in pepmassIonIter.items():
                if key in pepmassIonArray.keys():
                    pepmassIonArray[key] = pepmassIonArray[key] + value
                else:
                    pepmassIonArray[key] = value

        # sort the final mgfDf
        sortMgfDFValues(mgfDf)

        modCountDict = Counter()
        matchedPeptides = {}
        for key, value in massDict.items():
            # convert modified peptides to original form
            if not key.isalpha():
                alphaKey = modToPeptide(key)
            else:
                alphaKey = key

            for charge, chargeMass in value[2].items():
                # Shift to outside for charge for loop
                if alphaKey not in matchedPeptides.keys():

                    # define required data in a temporary form
                    pepMasses = mgfDf[charge]
                    #closestIndex = takeClosest(pepMasses, chargeMass, True)
                    #pepMass = pepMasses[closestIndex]

                    steps = [1,-1]
                    matchAdded = False
                    closestMatched = False

                    # need to iterate up and down from the closest to ensure all b/y ion comparison is run
                    # on all relevant pepMasses
                    for step in steps:
                        # matchAdded is set to True if if the B/y ion check is passed in the iterations
                        # using step 1
                        if matchAdded == True:
                            break

                        # if the closest pepmass doesn't pass the pepMatch test, closestMatched will still
                        # be false and we should break the loop to avoid running the code again.
                        if step == -1 and closestMatched == False:
                            break

                        # Get the closest index, denoted by the bool flag as the last argument, think this
                        # should be outside the step loop, as done twice when traversing back.
                        index = takeClosest(pepMasses, chargeMass, True)

                        # Start traversing backwards
                        if step == -1:
                            index += step

                        # Get the pepmass at that index
                        pepMass = pepMasses[index]

                        # While there is a current match. pepMass is changed per iteration if previous
                        # pepmass didn't match.
                        while pepMatch(chargeMass, pepMass, mgfObj.ppmVal):
                            # Set closestMatched flag to True if the while loop is entered so that the second
                            # step can be broken if closest match fails
                            closestMatched = True
                            # Not super important but this should likely be at the very start and outside the step
                            # loop.
                            if mgfObj.byIonFlag == False:
                                # count modifications
                                if not key.isalpha():
                                    modCountDict += getModNumbers(key, modList)
                                # if it is trans, massDict[3] will exist and will hold the desired protId
                                try:
                                    # create the string, with peptides sorted so all permutations are matched as similar. There may be multiple
                                    # peptide locations in the list of tuples, hence the for loop. Tuples are listed in order, with consecutive
                                    # tuples relating to a pair of splice locations.
                                    string = ""
                                    # create seenOrigins list to store origins which have been added so that there is no double up of identical origins.
                                    seenOrigins = []
                                    # iterate through all the origins stored in the massDict for the given key
                                    for i in range(0, len(massDict[alphaKey][3]), 2):
                                        # sort the origin so that it can be compared to others accurately
                                        origProt = sorted(massDict[alphaKey][3][i:i + 2])
                                        # if origProt has already been seen (and in turn added to seenOrigins) continue iterating.
                                        if origProt in seenOrigins:
                                            continue
                                        # if origProt hasn't been seen, append it to seenOrigins and add it to the string which is to be added to allPeptidesDict.
                                        seenOrigins.append(origProt)
                                        string += origProt[0][0] + origProt[0][1] + '/' + origProt[1][0] + origProt[1][1] + ';'
                                    string = string[0:-1]
                                    matchedPeptides[alphaKey] = string
                                    matchAdded = True
                                except IndexError:
                                    matchedPeptides[alphaKey] = protId
                                    matchAdded = True
                                break
                            else:
                                # Check the similarity of the byIons as was being done previously
                                # if closestMatched == False: create byIonArray, ptherwise it has already been created
                                byIonArray = initIonMass(key, modList)
                                mzArray = pepmassIonArray[(charge, pepMass)]
                                # If they match in accordance with the input minimum requirement, add them to the list
                                if simIons(mzArray, byIonArray, mgfObj.byIonAccuracy, mgfObj.minSimBy):
                                    # count the modifications
                                    if not key.isalpha():
                                        modCountDict += getModNumbers(key, modList)
                                    # if it is trans, massDict[3] will exist and will hold the desired protId
                                    try:
                                        # create the string, with peptides sorted so all permutations are matched as similar. There may be multiple
                                        # peptide locations in the list of tuples, hence the for loop. Tuples are listed in order, with consecutive
                                        # tuples relating to a pair of splice locations.
                                        string = ""
                                        # create seenOrigins list to store origins which have been added so that there is no double up of identical origins.
                                        seenOrigins = []
                                        # iterate through all the origins stored in the massDict for the given key
                                        for i in range(0, len(massDict[key][3]), 2):
                                            # sort the origin so that it can be compared to others accurately
                                            origProt = sorted(massDict[key][3][i:i + 2])
                                            # if origProt has already been seen (and in turn added to seenOrigins) continue iterating.
                                            if origProt in seenOrigins:
                                                continue
                                            # if origProt hasn't been seen, append it to seenOrigins and add it to the string which is to be added to allPeptidesDict.
                                            seenOrigins.append(origProt)
                                            string += origProt[0][0] + origProt[0][1] + '/' + origProt[1][0] + \
                                                      origProt[1][1] + ';'
                                        string = string[0:-1]
                                        matchedPeptides[alphaKey] = string
                                        matchAdded = True

                                    except IndexError:
                                        matchedPeptides[alphaKey] = protId
                                        matchAdded = True
                                    break

                                # If they didn't match try the next one. Step will be 1 when traversing forward, -1
                                # when traversing backward thus will be able to go up and down.
                                else:
                                    index += step
                                    try:
                                        pepMass = pepMasses[index]
                                    except IndexError:
                                        print('index is outide range of pepMasses')
                                        break
                else:
                    break
        print(modCountDict)
        return matchedPeptides, modCountDict

def getModNumbers(peptide, modList):
    modNumberDict = {}
    for character in peptide:
        if character.islower():
            for mod in modList:
                if character.upper in modTable[mod]:
                    modType = mod + character.upper()
                    try:
                        modNumberDict[modType] += 1
                    except KeyError:
                        modNumberDict[modType] = 1
    return modNumberDict


def getModNumbers(peptide, modList):
    modCountDict = Counter()
    for i in range(0, len(peptide)):
        if peptide[i].isdigit():
            continue
        if peptide[i].islower():
            modNumber = int(peptide[i+1]) - 1
            mod = modList[modNumber]
            moddedAmino = peptide[i].upper()

            if moddedAmino in modTable[mod]:

                modType = mod.split(' ')[0] + ' '+ moddedAmino + ' modified'
                try:
                    modCountDict[modType] += 1
                except KeyError:
                    modCountDict[modType] = 1
    return modCountDict

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



def readMGF(input_path, intensityThreshold, byIonFlag):

    """
    Creates a pandas dataframe based on mgf data
    """
    uniqueSpec = set()

    mgfDf = {}
    mgfDfList = []

    pepmassIonArray = {}
    pepmassIonArrayList = []

    with mgf.read(input_path) as mgfReader:
        counter = 0
        for spectrum in mgfReader:
            counter += 1

            if counter == 200000:
                pepmassIonArrayList.append(pepmassIonArray)
                pepmassIonArray = {}
                sortMgfDFValues(mgfDf)
                mgfDfList.append(mgfDf)
                mgfDf = {}

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

                if byIonFlag:
                    if chargePepmassTup in pepmassIonArray:
                        pepmassIonArray[chargePepmassTup].append(mzArray)
                        #print(pepmassIonArray[chargePepmassTup])
                    else:
                        pepmassIonArray[chargePepmassTup] = [mzArray]


                        # mgfDf.loc[len(mgfDf)] = [spectrum['params']['charge'][0],
                        #                          spectrum['params']['pepmass'][0]]


                uniqueSpec.add(chargePepmassTup)

    pepmassIonArrayList.append(pepmassIonArray)
    sortMgfDFValues(mgfDf)
    mgfDfList.append(mgfDf)

    return mgfDfList, pepmassIonArrayList, counter

def readMgfInit(input_path):

    maxIntensityArray = []

    with mgf.read(input_path) as mgfReader:
        for spectrum in mgfReader:

            if 'charge' in spectrum['params'].keys():


                maxIntensity = max(spectrum['intensity array'])
                maxIntensityArray.append(maxIntensity)


    return sorted(maxIntensityArray)

def changeIntToPoints(maxIntensityArray):
    ms2Thresh = [0,1,5,10,20,50,60,70,80,90,100,200,300,400,500,600,700, 800,
                 900,1000,5000,10000,20000,30000]

    maxIntLen = len(maxIntensityArray)
    intensityPoints = []


    for point in ms2Thresh:
        closest = findLargeIndex(maxIntensityArray, point)
        numAbove = maxIntLen - closest

        numAbove = (numAbove*100)/maxIntLen
        intensityPoints.append(numAbove)
    return ms2Thresh, intensityPoints

def findLargeIndex(arr,x):

    closest = takeClosest(arr, x, True)

    if closest == len(arr)-1 or closest == -1:

        if arr[-1] < x:

            return len(arr)
        else:
            return len(arr)-1


    for i in range(closest, len(arr)-1):
        if arr[i] > x:
            return i
def plotData(input_path):
    maxIntensityArray = readMgfInit(input_path)
    ms2Thresh, intenistyPoints = changeIntToPoints(maxIntensityArray)
    return ms2Thresh, intenistyPoints

def plot(ms2Thresh, intensityPoints):
    #ms2Thresh, intensityPoints = plotData(input_path)
    plt.figure()
    #plt.xlim([0, ms2Thresh[-1]])
    plt.ylim([0, 110])
    plt.xscale('log')
    plt.xlabel("Max Intensity Threshold")
    plt.ylabel("Percentage (%)")
    plt.plot(ms2Thresh, intensityPoints)
    plt.show()


def takeClosest(myList, myNumber, indexBool = False):
    """
    Assumes myList is sorted. Returns closest value to myNumber via index.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        if indexBool:
            return 0

        return myList[0]
    if pos == len(myList):

        if indexBool:
            return -1

        return myList[-1]

    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        if indexBool:
            return pos
        return after

    else:
        if indexBool:
            return pos - 1
        return before


def sortMgfDFValues(mgfDf):
    for charge, masses in mgfDf.items():
        masses.sort()


#plot('C:/Users/Arpit/Desktop/UROP2/InputData/MgfExample.mgf')
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
    return dict.values()

def sortBYDict(byIonDict):
    byIonArray = []
    for key, value in byIonDict.items():
        byIonArray.append(value)
    byIonArray.sort()
    return byIonArray

def findSimIons(mzArray, byIons, accuracy):

    simIons = 0
    for array in mzArray:
        simTemp = 0
        for mass in byIons:
            closest = takeClosest(array, mass)
            upperThresh = mass + accuracy
            lowerThresh = mass - accuracy
            if lowerThresh < closest < upperThresh:
                simTemp += 1
        if simTemp > simIons:
            simIons = simTemp

    simIons = (simIons*100) / len (byIons)

    return simIons

def simIons(mzArray, byIons, accuracy, minSim):
    noSimReq = math.ceil(len(byIons)*minSim/100)
    for array in mzArray:
        simTemp = 0
        byIonsTested = 0
        for mass in byIons:
            byIonsTested += 1
            closest = takeClosest(array, mass)
            upperThresh = mass + accuracy
            lowerThresh = mass - accuracy
            if lowerThresh < closest < upperThresh:
                simTemp += 1
                if simTemp >= noSimReq:
                    return True
            elif (len(byIons) - byIonsTested) < noSimReq - simTemp:
                break
    return False

# mzArray = [[8, 10, 12, 16, 20], [2, 4, 6, 39, 35]]
# byIons = [8.01, 9.99, 12.01, 20, 22]
# accuracy = 0.05
# minSim = 80
#
# print(simIons(mzArray, byIons, accuracy, minSim))

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

