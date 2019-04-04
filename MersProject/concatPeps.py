from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import os
import math

OUTPUT_PATH = 'example6-14_Linear_1_040419_0925_NoSubsets.fasta'
NO_RECORDS = 4000

def createPepList(OUTPUT_PATH):
    pepList = []
    with open(OUTPUT_PATH, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            pepList.append(str(record.seq))
        pepList.sort()
    return pepList

def overlapList(peptideList):
    i = 0
    while(True):
        try:
            createOverlap(peptideList, i)
            i+=1
        except IndexError:
            break
    return peptideList

def createOverlap(peptideList,i):
    # store peptide and check that a number hasn't been added to the front of it.
    peptide = peptideList[i]
    if not peptide.isalpha():
        print('woooo')
        return

    # loop through the different suffixes
    for j in range(1,len(peptide)):
        # extract suffix
        suffix = peptide[j:]
        # extract location of matching prefixIndexfrom the list, if there is one.
        prefixIndex = findSuff(suffix, peptideList, 0)
        if prefixIndex != False:
            # store the prefixPeptide, and replace its location in the list with None
            prefixPeptide = peptideList[prefixIndex]
            peptideList[prefixIndex] = prefixPeptide + '1'
            print(prefixPeptide)
            overlapPeptide = concatOverlapPep(peptide, j, prefixPeptide)
            # replace the current peptide with the new, overlapping peptide.
            peptideList[i] = overlapPeptide
            return

def findSuff(suffix, peptideList, zeroIndex):
    if len(peptideList) == 1:
        guessPeptide = peptideList[0]
        # check that the peptide we are up to is not None. If it is None, we want to return False as there are no
        # peptides left to iterate through
        if guessPeptide[0:len(suffix)] == suffix and guessPeptide.isalpha():
            return zeroIndex
        else:
            return False
    else:
        index = math.ceil(len(peptideList)/2)
        guessPeptide = peptideList[index]
        if guessPeptide[0:len(suffix)] == suffix and guessPeptide.isalpha():
            return zeroIndex+index
        guessPair = [suffix, guessPeptide]
        if sorted(guessPair)[0] == suffix:
            smallerPeptideList = peptideList[0:index]
            return findSuff(suffix, smallerPeptideList, zeroIndex)
        else:
            smallerPeptideList = peptideList[index:]
            return findSuff(suffix, smallerPeptideList, index+zeroIndex)


def concatOverlapPep(peptide, j, prefixPeptide):
    concatPep = peptide[0:j] + prefixPeptide
    return concatPep

def createSeqObj(seenPeptides):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []

    for sequence in seenPeptides:
        if not sequence.isalpha():
            continue
        finalId = "ipd|pep" + str(count) + ';'
        yield SeqRecord(Seq(sequence), id=finalId, description="")
        count += 1
    return seqRecords

def checkOutput(inputFile, outputFile):
    print("\nCHECKING OUTPUT")
    peptideList = []
    with open(outputFile, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            peptideList.append(str(record.seq))

    with open(inputFile, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            flag = True
            for peptide in peptideList:
                if str(record.seq) in peptide:
                    flag = False
                    break
            if flag:
                print(record.seq)

pepList = createPepList(OUTPUT_PATH)
print(len(pepList))
concatPepList = overlapList(pepList)
print(len(concatPepList))

# write to file
with open('concatOutput.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(concatPepList), output_handle, "fasta")

checkOutput(OUTPUT_PATH, 'concatOutput.fasta')

