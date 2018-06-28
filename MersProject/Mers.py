from Bio import SeqIO
import json
import csv
import time


class Fasta:

    def __init__(self, seqDict):
        self.seqDict = seqDict
        print(seqDict)

    min = 3
    '''
    Function that literally combines everything to generate output
    '''
    def generateOutput(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, modList, maxDistance, outputPath):
        if (transFlag):

            finalPeptide = combinePeptides(self.seqDict)
            massDict = genMassDict(finalPeptide, mined, maxed, overlapFlag, modList, maxDistance)
            writeToCsv(massDict, 'w', 'Combined',outputPath,'trans')

            # combined = {'combined': combined}
            #with open('output.txt', 'wb') as file:
             #   file.write(json.dumps(combined))  # use `json.loads` to do the reverse

        if (cisFlag):

            counter = 0
            for key, value in self.seqDict.items():
                massDict = genMassDict(value, mined, maxed, overlapFlag, modList, maxDistance)

                if (counter == 0):
                    writeToCsv(massDict, 'w', key, outputPath,'cis')
                    counter+=1

                else:
                    writeToCsv(massDict, 'a', key, outputPath,'cis')

        if (linearFlag):
            #linear dictionary function which converts splits and splits ref to the dictionary output desired
            counter = 0
            for key, value in self.seqDict.items():
                massDict = genMassLinear(value, mined, maxed, linearFlag, modList)

                if (counter == 0):
                    writeToCsv(massDict, 'w', key, outputPath, 'linear')
                    counter += 1

                else:
                    writeToCsv(massDict, 'a', key, outputPath, 'linear')

def genMassDict(peptide, mined, maxed, overlapFlag, modList, maxDistance):

    combined, combinedRef = outputCreate(peptide, mined, maxed, overlapFlag, maxDistance)
    massDict = combMass(combined, combinedRef)
    massDict = applyMods(massDict, modList)

    print(len(massDict))
    return massDict

def genMassLinear(peptide, mined, maxed, linearFlag, modList):

    combined, combinedRef = splitDictPeptide(peptide, mined, maxed, linearFlag)
    combined, combinedRef = removeDupsQuick(combined, combinedRef)
    massDict = combMass(combined, combinedRef)

    massDict = applyMods(massDict, modList)

    print(massDict)
    return massDict


def writeToCsv(massDict, writeFlag, header, outputPath, linkType):
    finalPath = str(outputPath) + '/' + linkType + '.csv'
    with open(finalPath, writeFlag, newline='') as csv_file:
        print(finalPath)
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow([header, ' ', ' '])
        writer.writerow(['Peptide', 'Mass', 'Positions'])
        for key, value in massDict.items():
            writer.writerow([key, value[0], value[1]])

# taking FASTA dictionary and passing through our splits and combine functions
#sequenceDictionary = addSequenceList("Example.fasta")

def outputCreate(peptide, mined, maxed, overlapFlag, maxDistance=None, linearFlag = False):

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

# MASS OF H20 is added for mass calculations
H20_MASS = 18.010565

# All possible modifications
modTable = {

    '4-hydroxynonenal (HNE)': ['C', 'H', 'K', 156.11504],
    'Acetylation (K)': ['K', 42.010567],
    'Beta-methylthiolation': ['C', 45.98772],
    'Carbamidomethylation': ['C', 57.021465],
    'Carboxylation (E)': ['E', 43.98983],
    'Carboxymethyl': ['C', 58.005478],
    'Citrullination': ['R', 0.984016],
    'Deamidation (NQ)': ['N', 'Q', 0.984016],
    'Dimethylation(KR)': ['K', 'R', 28.0313],
    'Dioxidation (M)': ['M', 31.989828],
    'FAD': ['C', 'H', 'Y', 783.1415],
    'Farnesylation': ['C', 204.1878],
    'Geranyl-geranyl': ['C', 272.2504 ],
    'Guanidination': ['K', 42.021797],
    'HexNAcylation (N)': ['N', 203.07938],
    'Hexose (NSY)': ['N', 'S', 'Y', 162.0528],
    'Lipoyl': ['K', 188.03296],
    'Methylation(KR)': ['K', 'R', 14.01565],
    'Methylation(others)': ['T', 'S', 'C', 'H', 'D', 'E', 14.01565],
    'Oxidation (HW)': ['H', 'W', 15.994915],
    'Oxidation (M)': ['M',15.994915],
    'Palmitoylation': ['C', 'S', 'T', 'K', 238.22966],
    'Phosphopantetheine': ['S', 340.0858],
    'Phosphorylation (HCDR)': ['H', 'C', 'D', 'R', 79.96633],
    'Phosphorylation (STY)': ['S', 'T', 'Y', 79.96633],
    'Propionamide': ['C', 71.03712],
    'Pyridoxal phosphate': ['K', 229.014],
    'S-pyridylethylation': ['C', 105.057846],
    'Sulfation': ['Y', 'S', 'T', 79.95682],
    'Sulphone': ['M', 31.989828],
    'Ubiquitin': ['T', 'S', 'C', 'K', 114.04293],
    'Ubiquitination': ['K', 383.2281],


}

# Mono-isotopic mass
monoAminoMass = {
    'A': 71.03711,
    'R': 156.10111,
    'N': 114.0493,
    'D': 115.02694,
    'C': 103.00919,
    'E': 129.04259,
    'Q': 128.05858,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'L': 113.08406,
    'K': 128.09496,
    'M': 131.04049,
    'F': 147.06841,
    'P': 97.05276,
    'S': 87.03203,
    'T': 101.04768,
    'W': 186.07931,
    'Y': 163.06333,
    'V': 99.06841,

}



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

"""
From the modless dictionary of possible combinations, this function returns a 
dictionary containing all the modifications that arise from changing character. The
key of the output is simply the modified peptide, and the value is the mass which
results as set by massChange
"""
def genericMod(combineModlessDict, character, massChange, modNo):
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

                    modDict[temp] = [newMass,currentRef]

    return modDict


"""
Inputs: peptide string, max length of split peptide. 
Outputs: all possible splits that could be formed that are smaller in length than the maxed input 
"""

def splitDictPeptide(peptide, mined, maxed, linearFlag):
    length = len(peptide)
    # splits will hold all possible splits that can occur
    splits = []
    # splitRef will hold a direct reference to the characters making up each split string: for starting peptide ABC,
    # the split AC = [0,2]
    splitRef = []

    # imbedded for loops build all possible splits
    for i in range(0, length):
        character = peptide[i]
        toAdd = ""
        # add and append first character and add and append reference number which indexes this character
        toAdd += character
        splits.append(toAdd)
        ref = []
        ref.append(i)
        temp = list(ref)  # use list because otherwise shared memory overwrites
        splitRef.append(temp)

        # iterates through every character after current and adds it to the most recent string if max size
        # requirement is satisfied
        for j in range(i + 1, length):
            toAdd += peptide[j]
            if linearFlag:
                if maxSize(toAdd, maxed) and minSize(toAdd, mined):
                    splits.append(toAdd)
                    ref.append(j)
                    temp = list(ref)

                    splitRef.append(temp)
                    temp = []
            else:
                if (maxSize(toAdd, maxed)):
                    splits.append(toAdd)
                    ref.append(j)
                    temp = list(ref)

                    splitRef.append(temp)
                    temp = []

        ref = []

    return splits, splitRef



"""Input: splits: list of splits, splitRef: list of the character indexes for splits, mined/maxed: min and max
   size requirements, overlapFlag: boolean value true if overlapping combinations are undesired.
   Output: all combinations of possible splits which meets criteria"""


def combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance):
    # initialise combinations array to hold the possible combinations from the input splits
    combModless = []
    combModlessRef = []
    # iterate through all of the splits and build up combinations which meet min/max/overlap criteria
    for i in range(0, len(splits)):

        # toAdd variables hold temporary combinations for insertion in final matrix if it meets criteria
        toAddForward = ""
        addForwardRef = []
        toAddReverse = ""
        addReverseRef = []
        for j in range(i + 1, len(splits)):
            # create forward combinaiton of i and j
            toAddForward += splits[i]
            toAddForward += splits[j]
            addForwardRef = splitRef[i] + splitRef[j]
            toAddReverse += splits[j]
            toAddReverse += splits[i]
            addReverseRef = splitRef[j] + splitRef[i]

            # max, min and max distance checks combined into one function for clarity for clarity
            if combineCheck(toAddForward,mined,maxed,splitRef[i],splitRef[j],maxDistance):
                if (overlapFlag == True):
                    if (overlapComp(splitRef[i], splitRef[j])):
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
            addForwardRef = []
            addReverseRef = []
    return combModless, combModlessRef


def maxDistCheck(ref1, ref2, maxDistance):
    if (maxDistance == 'None'):
        return True
    valid = ref2[-1] - ref1[0]

    if (valid > maxDistance):
        return False

    return True

# ensures length of split is smaller than or equal to max
def maxSize(split, maxed):
    if (len(split) >= maxed):
        return False
    return True


# ensures length of split is greater than min
def minSize(split, mined):
    if (len(split) < mined):
        return False
    return True

def combineCheck(split, mined, maxed, ref1, ref2, maxDistance):
    booleanCheck = maxSize(split,maxed) and minSize(split,mined) and maxDistCheck(ref1,ref2,maxDistance)
    return booleanCheck



# checks if there is an intersection between two strings. Likely input it the splitRef data.
# Outputs True if no intersection
# overlapComp needs to delete paired reference number if being applied to the splits output
def overlapComp(ref1, ref2):
    S1 = set(ref1)
    S2 = set(ref2)
    if (len(S1.intersection(S2)) == 0):
        return True
    return False





# opens FASTA file
def addSequenceList(input_file):
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    sequenceDictionary = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequenceDictionary[name] = sequence
    return sequenceDictionary


# combines an array of strings into one string. Used for ultimately segments from multiple peptides
def combinePeptides(dict):
    dictlist = []
    for key, value in dict.items():
        dictlist.append(value)

    finalPeptide = ''.join(dictlist)
    return finalPeptide


def removeDupsQuick(seq, seqRef):
    #seq = [''.join(sorted(s)) for s in seq]
    seen = set()
    seen_add = seen.add

    initial = []
    initialRef = []
    #initial = [x for x in seq if not (x in seen or seen_add(x))]
    for i in range(0, len(seq)):
        if not (seq[i] in seen or seen_add(seq[i])):
            initial.append(seq[i])
            initialRef.append(seqRef[i])
    print(initial)
    print(initialRef)


    return initial, initialRef




def combMass(combine, combineRef):
    massDict = {}
    for i in range(0, len(combine)):
        totalMass = 0
        for j in range(0, len(combine[i])):
            totalMass+=monoAminoMass[combine[i][j]]
        totalMass+=H20_MASS
        massRefPair = [totalMass, combineRef[i]]
        massDict[combine[i]] = massRefPair
    return massDict



# https://stackoverflow.com/questions/35091557/replace-nth-occurrence-of-substring-in-string
def nth_replace(string, old, new, n=1, option='only nth'):
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


# Analysis of removing duplicates at split level and combined level:
#     - Split level: 3.1 mil to 1.8 mil
#     - Combined level: 1.8 mil to 1.6 mil

