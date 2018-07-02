from Bio import SeqIO
import csv
from MonoAminoAndMods import *
import threading
import time

TRANS = "Trans"
LINEAR = "Linear"
CIS = "Cis"

class FileThread(threading.Thread):


    def __init__(self, spliceType, seqDict, mined, maxed, overlapFlag, modList, maxDistance, outputPath,
                 chargeFlags):
        threading.Thread.__init__(self)
        self.spliceType = spliceType
        self.seqDict = seqDict
        self.mined = mined
        self.maxed = maxed
        self.overlapFlag = overlapFlag
        self.modList = modList
        self.maxDistance = maxDistance
        self.outputPath = outputPath
        self.chargeFlags = chargeFlags
        self._is_running = True

    def run(self):
        print("Starting thread for: " + self.spliceType)

        while self._is_running:

            if self.spliceType == TRANS:

                finalPeptide = combinePeptides(self.seqDict)
                transOutput(finalPeptide, self.mined, self.maxed, self.overlapFlag, self.modList, self.maxDistance,
                            self.outputPath, self.chargeFlags)
                self.stop()

            elif self.spliceType == CIS:

                cisOutput(self.seqDict, self.mined, self.maxed, self.overlapFlag, self.modList, self.maxDistance,
                          self.outputPath, self.chargeFlags)
                self.stop()

            elif self.spliceType == LINEAR:
                linearOutput(self.seqDict, self.mined, self.maxed, self.modList, self.outputPath, self.chargeFlags)
                self.stop()
                
    def stop(self):
        self._is_running = False
        print('Closed Thread!')


class Fasta:

    def __init__(self, seqDict):
        self.seqDict = seqDict

    def generateOutput(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList,
                       maxDistance, outputPath, chargeFlags):

        """
           Function that literally combines everything to generate output
        """

        if transFlag:
            transThread = FileThread(TRANS, self.seqDict, mined, maxed, overlapFlag, modList, maxDistance, outputPath,
                                     chargeFlags)
            transThread.start()
            #finalPeptide = combinePeptides(self.seqDict)
            #transOutput(finalPeptide, mined, maxed, overlapFlag, modList, maxDistance, outputPath, chargeFlags)


            # combined = {'combined': combined}
            # with open('output.txt', 'wb') as file:
            # file.write(json.dumps(combined))  # use `json.loads` to do the reverse

        if cisFlag:
            cisThread = FileThread(CIS, self.seqDict, mined, maxed, overlapFlag, modList, maxDistance, outputPath,
                                   chargeFlags)
            cisThread.start()

        if linearFlag:

            linearThread = FileThread(LINEAR, self.seqDict, mined, maxed, overlapFlag, modList, maxDistance,
                                      outputPath, chargeFlags)
            linearThread.start()



def transOutput(finalPeptide, mined, maxed, overlapFlag, modList, maxDistance, outputPath, chargeFlags):
    massDict = genMassDict(finalPeptide, mined, maxed, overlapFlag, modList, maxDistance, chargeFlags)
    writeToCsv(massDict, 'w', 'Combined', outputPath, 'trans', chargeFlags)


def cisOutput(seqDict, mined, maxed, overlapFlag, modList, maxDistance, outputPath, chargeFlags):
    counter = 0
    for key, value in seqDict.items():
        massDict = genMassDict(value, mined, maxed, overlapFlag, modList, maxDistance, chargeFlags)

        if counter == 0:
            writeToCsv(massDict, 'w', key, outputPath, 'cis', chargeFlags)
            counter += 1

        else:
            writeToCsv(massDict, 'a', key, outputPath, 'cis', chargeFlags)

def linearOutput(seqDict, mined, maxed, modList, outputPath, chargeFlags):
    # linear dictionary function which converts splits and splits ref to the dictionary output desired
    counter = 0
    for key, value in seqDict.items():
        massDict = genMassLinear(value, mined, maxed, modList, chargeFlags)

        if counter == 0:
            writeToCsv(massDict, 'w', key, outputPath, 'linear', chargeFlags)
            counter += 1

        else:
            writeToCsv(massDict, 'a', key, outputPath, 'linear', chargeFlags)

def genMassDict(peptide, mined, maxed, overlapFlag, modList, maxDistance, chargeFlags):

    combined, combinedRef = outputCreate(peptide, mined, maxed, overlapFlag, maxDistance)
    massDict = combMass(combined, combinedRef)
    massDict = applyMods(massDict, modList)
    chargeIonMass(massDict, chargeFlags)
    return massDict


def genMassLinear(peptide, mined, maxed, modList, chargeFlags):
    linearFlag = True
    combined, combinedRef = splitDictPeptide(peptide, mined, maxed, linearFlag)
    combined, combinedRef = removeDupsQuick(combined, combinedRef)
    massDict = combMass(combined, combinedRef)
    massDict = applyMods(massDict, modList)
    chargeIonMass(massDict, chargeFlags)
    return massDict


def chargeIonMass(massDict, chargeFlags):

    """
    chargeFlags: [True, False, True, False, True]
    """

    for key, value in massDict.items():
        chargeAssoc = {}
        for z in range(0, len(chargeFlags)):

            if chargeFlags[z]:
                chargeMass = massCharge(value[0], z+1) # +1 for actual value
                chargeAssoc[z+1] = chargeMass
        value.append(chargeAssoc)


def massCharge(predictedMass, z):
    chargeMass = (predictedMass + (z * 1.00794))/z
    return chargeMass


def writeToCsv(massDict, writeFlag, header, outputPath, linkType, chargeFlags):

    finalPath = str(outputPath) + '/' + linkType + '.csv'

    chargeHeaders = getChargeIndex(chargeFlags)

    with open(finalPath, writeFlag, newline='') as csv_file:


        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow([header, ' ', ' '])
        headerRow = ['Peptide', 'Mass', 'Positions']

        for chargeIndex in chargeHeaders:

            headerRow.append('+' + str(chargeIndex+1))

        writer.writerow(headerRow)
        for key, value in massDict.items():
            infoRow = [key, value[0], value[1]]
            for chargeIndex in chargeHeaders:
                chargeMass = value[2][chargeIndex+1]
                infoRow.append(str(chargeMass))
            writer.writerow(infoRow)


def getChargeIndex(chargeFlags):
    chargeHeaders = [i for i, e in enumerate(chargeFlags) if e]
    return chargeHeaders


def outputCreate(peptide, mined, maxed, overlapFlag, maxDistance=None, linearFlag=False):

    # Splits eg: ['A', 'AB', 'AD', 'B', 'BD']
    # SplitRef eg: [[0], [0,1], [0,2], [1], [1,2]]
    # Produces splits and splitRef arrays which are passed through combined
    splits, splitRef = splitDictPeptide(peptide, mined, maxed, linearFlag)

    # combined eg: ['ABC', 'BCA', 'ACD', 'DCA']
    # combinedRef eg: [[0,1,2], [1,0,2], [0,2,3], [3,2,0]]
    # pass splits through combined overlap peptide and then delete all duplicates
    combined, combinedRef = combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance)
    combined, combinedRef = removeDupsQuick(combined, combinedRef)

    return combined, combinedRef


def applyMods(combineModlessDict, modList):

    """
    Calls the genericMod function and accesses the modification table to
    append modified combinations to the end of the combination dictionary
    """

    modNo = 0
    for mod in modList:
        # Keep track of which modification is taking place
        modNo += 1

        # Don't need to worry about it if no modification!
        if mod != 'None':
            # Get the list of modifications taking place
            aminoList = modTable[mod]
            # Go through each character in the modification one by one
            for i in range(0, len(aminoList) - 1):

                char = aminoList[i]
                massChange = aminoList[-1]
                # get the dictionary of mods and their mass
                modDict = genericMod(combineModlessDict, char, massChange, str(modNo))
                # Add it to the current list!
                combineModlessDict.update(modDict)
    return combineModlessDict


def genericMod(combineModlessDict, character, massChange, modNo):

    """
    From the modless dictionary of possible combinations, this function returns a
    dictionary containing all the modifications that arise from changing character. The
    key of the output is simply the modified peptide, and the value is the mass which
    results as set by massChange
    """
    # A, B, C  convert to ai, bi, ci where i is the modNo
    modDict = {}

    # Go through each combination and mod it if necessary
    for string in combineModlessDict.keys():

        currentMass = combineModlessDict[string][0]
        currentRef = combineModlessDict[string][1]

        # Only need to mod it if it exists (ie : A in ABC)
        if character in string:

            numOccur = string.count(character)
            # Generate all permutations with the mods
            for j in range(0, numOccur):
                temp = string
                for i in range(0, numOccur - j):
                    newMass = currentMass + (i + 1) * massChange
                    temp = nth_replace(temp, character, character.lower() + modNo, j + 1)

                    modDict[temp] = [newMass, currentRef]

    return modDict


def splitDictPeptide(peptide, mined, maxed, linearFlag):

    """
    Inputs: peptide string, max length of split peptide.
    Outputs: all possible splits that could be formed that are smaller in length than the maxed input
    """

    length = len(peptide)

    # splits will hold all possible splits that can occur
    splits = []
    # splitRef will hold a direct reference to the characters making up each split string: for starting peptide ABC,
    # the split AC = [0,2]
    splitRef = []

    # embedded for loops build all possible splits
    for i in range(0, length):
        character = peptide[i]
        toAdd = ""
        # add and append first character and add and append reference number which indexes this character

        toAdd += character
        ref = list([i])
        temp = list(ref)  # use list because otherwise shared memory overwrites

        # linear flag to ensure min is correct for cis and trans
        if linearFlag and minSize(toAdd, mined):
            splits.append(toAdd)
            splitRef.append(temp)

        elif not linearFlag:
            splits.append(toAdd)
            splitRef.append(temp)

        # iterates through every character after current and adds it to the most recent string if max size
        # requirement is satisfied
        for j in range(i + 1, length):
            toAdd += peptide[j]
            if linearFlag:
                ref.append(j)
                if maxSize(toAdd, maxed) and minSize(toAdd, mined):
                    splits.append(toAdd)
                    temp = list(ref)
                    splitRef.append(temp)

            else:
                if maxSize(toAdd, maxed-1):
                    splits.append(toAdd)
                    ref.append(j)
                    temp = list(ref)
                    splitRef.append(temp)

    return splits, splitRef


def combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance):

    """
    Input: splits: list of splits, splitRef: list of the character indexes for splits, mined/maxed: min and max
    size requirements, overlapFlag: boolean value true if overlapping combinations are undesired.
    Output: all combinations of possible splits which meets criteria
    """

    # initialise combinations array to hold the possible combinations from the input splits
    combModless = []
    combModlessRef = []
    # iterate through all of the splits and build up combinations which meet min/max/overlap criteria
    for i in range(0, len(splits)):

        # toAdd variables hold temporary combinations for insertion in final matrix if it meets criteria
        toAddForward = ""
        # addForwardRef = []
        toAddReverse = ""
        # addReverseRef = []
        for j in range(i + 1, len(splits)):
            # create forward combination of i and j
            toAddForward += splits[i]
            toAddForward += splits[j]
            addForwardRef = splitRef[i] + splitRef[j]
            toAddReverse += splits[j]
            toAddReverse += splits[i]
            addReverseRef = splitRef[j] + splitRef[i]

            # max, min and max distance checks combined into one function for clarity for clarity
            if combineCheck(toAddForward, mined, maxed, splitRef[i], splitRef[j], maxDistance):
                if overlapFlag:
                    if overlapComp(splitRef[i], splitRef[j]):
                        combModless.append(toAddForward)
                        combModless.append(toAddReverse)
                        combModlessRef.append(addForwardRef)
                        combModlessRef.append(addReverseRef)
                else:
                    combModless.append(toAddForward)
                    combModless.append(toAddReverse)
                    combModlessRef.append(addForwardRef)
                    combModlessRef.append(addReverseRef)

            toAddForward = ""
            toAddReverse = ""
            # addForwardRef = []
            # addReverseRef = []
    return combModless, combModlessRef


def maxDistCheck(ref1, ref2, maxDistance):
    if maxDistance == 'None':
        return True
    valid = ref2[-1] - ref1[0]

    if valid > maxDistance:
        return False
    return True


def maxSize(split, maxed):

    """
    ensures length of split is smaller than or equal to max
    """

    if len(split) > maxed:
        return False
    return True


def minSize(split, mined):

    """
    ensures length of split is greater than min
    """

    if len(split) < mined:
        return False
    return True


def combineCheck(split, mined, maxed, ref1, ref2, maxDistance):
    booleanCheck = maxSize(split, maxed) and minSize(split, mined) and maxDistCheck(ref1, ref2, maxDistance)
    return booleanCheck


def overlapComp(ref1, ref2):

    """
    checks if there is an intersection between two strings. Likely input it the splitRef data.
    Outputs True if no intersection
    overlapComp needs to delete paired reference number if being applied to the splits output
    """

    S1 = set(ref1)
    S2 = set(ref2)
    if len(S1.intersection(S2)) == 0:
        return True
    return False


def addSequenceList(input_file):

    """
    input_file is the file path to the fasta file.
    The function reads the fasta file into a dictionary with the format {proteinRef: protein}
    """

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    sequenceDictionary = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequenceDictionary[name] = sequence
    return sequenceDictionary


def combinePeptides(seqDict):

    """
    combines an array of strings into one string. Used for ultimately segments from multiple peptides
    """

    dictlist = []
    for key, value in seqDict.items():
        dictlist.append(value)

    finalPeptide = ''.join(dictlist)
    return finalPeptide


def removeDupsQuick(seq, seqRef):

    seen = set()
    seen_add = seen.add
    initial = []
    initialRef = []
    # initial = [x for x in seq if not (x in seen or seen_add(x))]
    for i in range(0, len(seq)):
        if not (seq[i] in seen or seen_add(seq[i])):
            initial.append(seq[i])
            initialRef.append(seqRef[i])

    return initial, initialRef


def combMass(combine, combineRef):
    massDict = {}
    for i in range(0, len(combine)):
        totalMass = 0
        for j in range(0, len(combine[i])):
            totalMass += monoAminoMass[combine[i][j]]
        totalMass += H20_MASS
        massRefPair = [totalMass, combineRef[i]]
        massDict[combine[i]] = massRefPair
    return massDict


def nth_replace(string, old, new, n=1, option='only nth'):

    # https://stackoverflow.com/questions/35091557/replace-nth-occurrence-of-substring-in-string
    """
    This function replaces occurrences of string 'old' with string 'new'.
    There are three types of replacement of string 'old':
    1) 'only nth' replaces only nth occurrence (default).
    2) 'all left' replaces nth occurrence and all occurrences to the left.
    3) 'all right' replaces nth occurrence and all occurrences to the right.
    """

    if option == 'only nth':
        left_join = old
        right_join = old
    elif option == 'all left':
        left_join = new
        right_join = old
    elif option == 'all right':
        left_join = old
        right_join = new
    else:
        print("Invalid option. Please choose from: 'only nth' (default), 'all left' or 'all right'")
        return None
    groups = string.split(old)
    nth_split = [left_join.join(groups[:n]), right_join.join(groups[n:])]
    return new.join(nth_split)
