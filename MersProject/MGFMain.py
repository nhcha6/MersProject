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
import math

class MGF:

    """
    Class to represent MGF input data
    """


    def __init__(self, mgfDf, pepmassIonArray, ppmVal, intensityThreshold, minSimBy, byIonAccuracy, byIonFlag, maxMass,
                 chargeMaxDict):

        # mgfDf looks like: {'charge': [list of masses]}
        self.mgfDf = mgfDf
        self.pepmassIonArray = pepmassIonArray
        self.mgfEntries = len(mgfDf)
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
    if mgfObj.mgfDf:

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
                    pepMasses = mgfObj.mgfDf[charge]
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
                                # if it is trans, massDict[3] will exist and will hold the desired protId
                                try:
                                    string = ""
                                    for i in range(0, len(value[3]), 2):
                                        origProt = sorted(value[3][i:i + 2])
                                        string += origProt[0][0] + origProt[0][1] + '/' + origProt[1][0] + \
                                                  origProt[1][1] + ';'
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
                                mzArray = mgfObj.pepmassIonArray[(charge, pepMass)]
                                # If they match in accordance with the input minimum requirement, add them to the list
                                if simIons(mzArray, byIonArray, mgfObj.byIonAccuracy, mgfObj.minSimBy):
                                    # if it is trans, massDict[3] will exist and will hold the desired protId
                                    try:
                                        string = ""
                                        for i in range(0, len(value[3]), 2):
                                            origProt = sorted(value[3][i:i + 2])
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
                                    pepMass = pepMasses[index]
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
pepmassIonArray =
{(2, 1249.6011): [array([ 138.054,  143.152,  143.163,  169.14 ,  171.146,  171.162,
        171.169,  199.141,  199.148,  225.194,  229.115,  229.12 ,
        245.069,  245.079,  245.092,  253.19 ,  253.198,  268.167,
        272.156,  272.165,  284.127,  286.173,  292.128,  292.14 ,
        298.134,  298.141,  298.151,  302.143,  314.166,  314.178,
        316.142,  316.153,  326.137,  326.146,  328.187,  328.202,
        337.184,  337.207,  341.142,  344.145,  344.162,  344.171,
        345.175,  371.224,  371.239,  373.212,  373.228,  375.122,
        385.179,  395.154,  397.209,  398.208,  403.182,  409.178,
        409.194,  413.154,  413.168,  413.233,  415.212,  415.226,
        417.201,  425.19 ,  425.217,  429.203,  431.173,  431.187,
        443.211,  443.221,  443.231,  443.249,  453.115,  469.292,
        470.295,  471.232,  497.258,  500.268,  512.223,  512.236,
        514.283,  514.3  ,  516.266,  516.286,  524.279,  530.233,
        530.252,  542.265,  542.279,  542.301,  583.298,  583.311,
        583.328,  584.264,  593.303,  601.332,  603.241,  611.298,
        611.322,  629.298,  629.312,  629.325,  629.346,  629.357,
        635.258,  635.274,  635.308,  675.344,  680.299,  685.401,
        721.358,  728.345,  742.423,  770.378,  793.391,  815.405,
        815.468,  833.377,  833.391,  833.408,  833.42 ,  833.444,
        893.452,  930.48 , 1004.71 , 1006.535, 1007.999, 1066.461,
       1084.727, 1090.416, 1099.168, 1102.475, 1113.517, 1117.517,
       1119.471, 1120.604, 1123.576, 1129.49 , 1131.368, 1135.643,
       1137.445, 1138.578, 1142.465, 1142.67 , 1155.572, 1157.4  ,
       1159.417, 1161.594, 1168.552, 1170.65 , 1177.601, 1188.48 ,
       1193.071, 1205.416, 1213.525, 1220.573, 1224.489, 1231.663,
       1238.584, 1245.634, 1245.967, 1246.376, 1247.392, 1251.446,
       1259.629, 1259.749, 1263.522, 1264.6  , 1270.582, 1285.672,
       1286.613, 1295.641, 1296.57 , 1296.626, 1306.734, 1308.81 ,
       1324.649, 1327.618, 1330.529, 1332.171, 1333.706, 1336.969,
       1337.768, 1346.378, 1346.715, 1366.277, 1366.715, 1378.704,
       1379.763, 1380.774, 1383.649, 1384.651, 1390.71 , 1392.557,
       1393.758, 1398.507, 1398.591, 1400.787, 1406.231, 1419.774,
       1423.252, 1425.574, 1425.654, 1427.497, 1475.691, 1491.623,
       1496.672, 1501.681, 1524.801, 1526.856, 1536.695, 1544.596,
       1544.774, 1551.957, 1558.672, 1570.58 , 1579.812, 1591.761,
       1600.079, 1607.82 , 1611.844, 1624.931, 1650.793, 1658.615,
       1664.86 , 1665.965, 1690.764, 1710.217, 1714.882, 1728.356,
       1745.63 , 1745.978, 1746.007, 1753.791, 1762.842, 1774.755,
       1779.804, 1784.819, 2234.026, 2497.039])], (2, 1249.601): [array([  68.053,   70.069,   72.081,   72.096,   81.045,   86.098,
        126.053,  159.092,  159.115,  161.093,  167.085,  167.159,
        171.15 ,  185.162,  188.072,  199.142,  199.149,  199.159,
        217.087,  227.058,  227.068,  230.069,  241.05 ,  245.075,
        245.09 ,  251.245,  257.196,  261.155,  272.168,  275.097,
        292.127,  292.147,  292.162,  293.108,  300.194,  304.113,
        314.166,  314.177,  316.152,  316.169,  328.182,  328.192,
        331.16 ,  335.114,  335.129,  339.178,  344.14 ,  344.148,
        344.155,  344.165,  344.175,  345.172,  373.21 ,  379.2  ,
        395.148,  397.206,  398.189,  403.174,  403.179,  403.187,
        409.19 ,  409.202,  413.161,  413.168,  413.237,  415.216,
        415.224,  417.227,  423.196,  425.2  ,  431.171,  431.18 ,
        443.199,  443.207,  443.214,  443.223,  443.232,  443.241,
        443.247,  457.214,  457.232,  489.125,  496.269,  512.23 ,
        512.237,  514.277,  514.285,  514.291,  514.298,  516.268,
        530.242,  542.265,  542.275,  542.284,  542.293,  559.338,
        579.23 ,  579.276,  584.283,  584.309,  585.272,  586.321,
        601.293,  601.323,  603.257,  611.294,  611.305,  611.318,
        611.342,  615.295,  617.28 ,  619.282,  629.31 ,  629.317,
        629.331,  629.338,  629.358,  635.266,  635.276,  668.238,
        672.358,  683.271,  698.309,  743.348,  748.332,  765.372,
        785.424,  815.383,  815.402,  833.376,  833.39 ,  833.4  ,
        833.407,  833.422,  833.44 ,  922.419,  993.497, 1023.478,
       1027.484, 1050.651, 1186.578, 1231.659, 1257.586, 1310.667,
       1311.669, 1315.57 , 1319.549, 1319.564, 1330.522, 1337.757,
       1339.615, 1344.131, 1345.482, 1358.241, 1368.106, 1371.62 ,
       1372.675, 1374.699, 1380.569, 1380.794, 1392.809, 1397.198,
       1409.717, 1412.685, 1413.279, 1422.756, 1426.505, 1435.652,
       1438.571, 1439.588, 1443.825, 1444.463, 1446.668, 1454.679,
       1455.448, 1459.6  , 1459.757, 1464.806, 1469.659, 1469.734,
       1475.587, 1479.735, 1496.627, 1510.488, 1515.482, 1519.73 ,
       1524.808, 1524.837, 1524.871, 1524.918, 1525.807, 1525.84 ,
       1525.862, 1533.642, 1535.04 , 1546.763, 1555.792, 1559.69 ,
       1565.763, 1589.762, 1589.807, 1596.925, 1597.793, 1599.485,
       1612.811, 1617.681, 1617.789, 1624.878, 1628.747, 1636.579,
       1644.745, 1644.779, 1646.805, 1648.82 , 1649.898, 1655.788,
       1656.213, 1664.375, 1665.705, 1665.791, 1667.709, 1668.706,
       1669.841, 1671.68 , 1671.796, 1676.552, 1679.783, 1679.823,
       1680.234, 1722.647, 1725.81 , 1726.947, 1727.915, 1731.629,
       1732.791, 1732.956, 1736.699, 1741.782, 1744.892, 1744.945,
       1745.322, 1745.711, 1745.905, 1745.929, 1745.959, 1759.88 ,
       1768.691, 1769.777, 1780.773, 1789.266, 1795.687, 1797.732,
       2497.022, 2497.083, 2497.176, 2497.297, 2685.259, 2956.417])]}
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

