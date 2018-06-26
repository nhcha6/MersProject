from Bio import SeqIO




# taking FASTA dictionary and passing through our splits and combine functions
#sequenceDictionary = addSequenceList("Example.fasta")
'''
for key, value in sequenceDictionary.items():
    splits, splitRef = splitDictPeptide(value, maxed)

    splits = removeDupsQuick(splits)
    combine = combineOverlapPeptide(splits, splitRef, mined, maxed, overlap)
    combine = removeDupsQuick(combine)

    print(len(splits))
    print(len(combine))

    break;

'''

H20_MASS = 18.010565
maxed = 12
mined = 0
maxDistance = 25
overlap = True

# hello hello hello
# Monoisotopic mass
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

'''
Function that literally combines everything to generate output

'''
def generateOutput(peptide, mined, maxed, overlapFlag, maxDistance=None):

    # Produces splits and splitRef arrays which are passed through combined
    splits, splitRef = splitDictPeptide(peptide, maxed)
    #splits = removeDupsQuick(splits)
    print(splits)

    # pass splits through combined overlap peptide and then delete all duplicates
    combined = combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance)
    combined = removeDupsQuick(combined)

    return combined


"""
Inputs: peptide string, max length of split peptide. 
 Outputs: all possible splits that could be formed that are smaller in length than the maxed input
"""

"""Inputs: peptide string, max length of split peptide. 
   Outputs: all possible splits that could be formed that are smaller in length than the maxed input """



def splitDictPeptide(peptide, maxed):
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


def combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance=None):
    # initialise combinations array to hold the possible combinations from the input splits
    combine = []

    # iterate through all of the splits and build up combinations which meet min/max/overlap criteria
    for i in range(0, len(splits)):

        # toAdd variables hold temporary combinations for insertion in final matrix if it meets criteria
        toAddForward = ""
        toAddReverse = ""

        for j in range(i + 1, len(splits)):
            # create forward combinaiton of i and j
            toAddForward += splits[i]
            toAddForward += splits[j]
            toAddReverse += splits[j]
            toAddReverse += splits[i]

            # look to combine all checks together in a future for clarity
            if combineCheck(toAddForward,mined,maxed,splitRef[i],splitRef[j],maxDistance):

                if (overlapFlag == True):
                    if (overlapComp(splitRef[i], splitRef[j])):
                        combine.append(toAddForward)
                        combine.append(toAddReverse)
                else:
                    combine.append(toAddForward)
                    combine.append(toAddReverse)

            toAddForward = ""
            toAddReverse = ""
    return combine


def maxDistCheck(ref1, ref2, maxDistance):
    if (maxDistance == None):
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
        name, sequence = fasta.id, fasta.seq.tostring()
        sequenceDictionary[name] = sequence
    return sequenceDictionary


# combines an array of strings into one string. Used for ultimately segments from multiple peptides
def combinePeptides(peptideList):
    finalPeptide = ''.join(peptideList)
    return finalPeptide


def removeDupsQuick(seq):
    #seq = [''.join(sorted(s)) for s in seq]
    seen = set()
    seen_add = seen.add
    initial = [x for x in seq if not (x in seen or seen_add(x))]

    return initial


# generates most of the permutations possible when switching from A to a in all strings originally containing an A
# input is a list of all combinations
#need to look to create all possible combinations
def modTest(combine):
    # A, B, C  convert to a, b, c
    modComb = []
    for string in combine:
        if 'A' in string:

            numOccur = string.count('A')

            for i in range(0, numOccur):
                temp = string
                temp = temp.replace("A", "a", i + 1)

                modComb.append(temp)
    return modComb


def combMass(combine):
    massDict = {}
    for combination in combine:
        totalMass = 0
        for character in combination:
            totalMass += monoAminoMass[character]
        totalMass+=H20_MASS
        massDict[combination] = totalMass
    return massDict






# Analysis of removing duplicates at split level and combined level:
#     - Split level: 3.1 mil to 1.8 mil
#     - Combined level: 1.8 mil to 1.6 mil

print(generateOutput('ABC', mined, maxed, overlap))