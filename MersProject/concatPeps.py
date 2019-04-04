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
    # delete the peptide from the list
    peptide = peptideList[i]
    del peptideList[i]
    # loop through the different suffixes
    for j in range(0,len(peptide)):
        # extract suffix
        suffix = peptide[j:]
        #print(suffix)
        # extract matching prefixPeptide from the list, if there is one.
        prefixPeptide = findSuff(suffix, peptideList)
        #print(prefixPeptide)
        if prefixPeptide != False:
            peptideList.remove(prefixPeptide)
            overlapPeptide = concatOverlapPep(peptide, j, prefixPeptide)
            peptideList.insert(i,overlapPeptide)
            #print(overlapPeptide)
            return
    print(peptide)
    peptideList.insert(i, peptide)

def findSuff(suffix, peptideList):
    if len(peptideList) == 1:
        guessPeptide = peptideList[0]
        if guessPeptide[0:len(suffix)] == suffix:
            return guessPeptide
        else:
            return False
    else:
        index = math.ceil(len(peptideList)/2)
        guessPeptide = peptideList[index]
        if guessPeptide[0:len(suffix)] == suffix:
            return guessPeptide
        guessPair = [suffix, guessPeptide]
        if sorted(guessPair)[0] == suffix:
            smallerPeptideList = peptideList[0:index]
        else:
            smallerPeptideList = peptideList[index:]
        return findSuff(suffix, smallerPeptideList)

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
        finalId = "ipd|pep" + str(count) + ';'
        yield SeqRecord(Seq(sequence), id=finalId, description="")
        count += 1
    return seqRecords






pepList = createPepList(OUTPUT_PATH)
print(len(pepList))
concatPepList = overlapList(pepList)
print(len(concatPepList))

# write to file
with open('concatOutput.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(concatPepList), output_handle, "fasta")
