# Adapted from https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
# removes duplicates given a list
def removeDups(seq):
    seen = set()
    seen_add = seen.add
    initial = [x for x in seq if not (x in seen or seen_add(x))]
    temp = []
    #Remove permutations
    for i in range(0, len(initial)):
        for j in range(i + 1, len(initial)):
            if (combRemove(initial[i], initial[j])):
                temp.append(initial[j])

    final = [x for x in initial if x not in temp]
    return final

"""
Returns true if two strings are a permutation of each other -> Called in removeDups
"""


def combRemove(ref1, ref2):
    return (len(ref1) == len(ref2) and sorted(ref1) == sorted(ref2))

"""" combine function which directly output a dictionary value with the mass incorporated. Further logic on the
dictionary was very slow so the function was disregarded"""


def combineMassDict(splits, splitRef, mined, maxed, overlapFlag, maxDistance):
    # initialise combinations array to hold the possible combinations from the input splits
    modlessComb = {}

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
            mass = combMass(toAddForward)
            # look to combine all checks together in a future for clarity
            if combineCheck(toAddForward, mined, maxed, splitRef[i], splitRef[j], maxDistance):

                if (overlapFlag == True):
                    if (overlapComp(splitRef[i], splitRef[j])):
                        modlessComb[toAddForward] = mass
                        modlessComb[toAddReverse] = mass

                else:
                    modlessComb[toAddForward] = mass
                    modlessComb[toAddReverse] = mass

            toAddForward = ""
            toAddReverse = ""
    return modlessComb

"""Mass calculation function which takes a string instead of a list of strings. Was used for combineMassDict before it was rejected"""

def combMass(peptide):
    mass = 0
    for character in peptide:
        mass += monoAminoMass[character]
    mass+=H20_MASS
    return mass

# concatenates references, converts each element to a string, and calls function which outputs false if
# the list is sequential
def linearCheckForCis(ref1, ref2):
    concRef = ref1 + ref2
    concRef = [str(i) for i in concRef]
    return checkSequential(concRef)

def checkSequential(l):
    a = [int(i, 16) for i in sorted(set(l))]
    return not (len(a) == (a[-1] - a[0] + 1))