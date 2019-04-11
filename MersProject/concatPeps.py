from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import os
import math

OUTPUT_PATH1 = 'example6-14_Linear_1_040419_0925_NoSubsets.fasta'
OUTPUT_PATH = 'concatOutputStop.fasta'
OUTPUT_PATH2 = 'C:/Users/Administrator/Desktop/Remove Subseqs/a2Maxmods3-Linear050219_2324_NoSubsets.fasta'
NO_RECORDS = 4000

class ConcatList:

    def __init__(self, pepList):
        self.peptideList = pepList
        self.pepListLen = len(pepList)
        self.pepListLenOld = 0

    def overlapList(self):
        i = 0
        while(True):
            try:
                if i % 10000 == 0:
                    print(i)
                self.createOverlap(i)
                i+=1
            except IndexError:
                break
        return

    def createOverlap(self,i):
        # store peptide and check that a number hasn't been added to the front of it.
        peptide = self.peptideList[i]
        if peptide[-1]=='1':
            #print('woooo')
            return

        # as concated peptides get longer we do not want to iterate through all suffixes. Max suffix length is set
        # to be 20
        startIter = len(peptide) - 20
        if startIter < 1:
            startIter = 1

        # loop through the different suffixes
        for j in range(startIter,len(peptide)):
            # extract suffix
            suffix = peptide[j:]
            # extract location of matching prefixIndexfrom the list, if there is one.
            prefixIndex = self.findSuff(suffix, (0,self.pepListLen-1))

            if prefixIndex != None:
                # store the prefixPeptide, and replace its location in the list with None
                prefixPeptide = self.peptideList[prefixIndex]
                self.peptideList[prefixIndex] = prefixPeptide + '1'
                overlapPeptide = concatOverlapPep(peptide, j, prefixPeptide)
                # replace the current peptide with the new, overlapping peptide.
                self.peptideList[i] = overlapPeptide
                return

    # need to rewrite this function so that we do not create a subset of the peptideList each time. This is the only
    # place i can find that is causing inefficiencies.
    # find suff is essentially a binary search algorithm.
    def findSuff(self, suffix, range):
        if range[0] == range[1]:
            index = range[0]
            guessPeptide = self.peptideList[index]
            # check that the peptide we are up to is not None. If it is None, we want to return False as there are no
            # peptides left to iterate through
            if guessPeptide[0:len(suffix)] == suffix and guessPeptide[-1]!=1:
                return index
            else:
                return None
        else:
            index = math.ceil((range[1] - range[0]) / 2) + range[0]
            guessPeptide = self.peptideList[index]
            if guessPeptide[0:len(suffix)] == suffix:
                if guessPeptide[-1]!='1':
                    return index
                else:
                    # iterate forward through self.peptideList and try find a match without the 1 on the end.
                    j = 1
                    while True:
                        guessPeptide = self.peptideList[index+j]
                        if guessPeptide[0:len(suffix)]==suffix:
                            if guessPeptide[-1]!='1':
                                return index+j
                            else:
                                j+=1
                        else:
                            break
                    # iterate backwards if the forward iteration failed.
                    j = -1
                    while True:
                        guessPeptide = self.peptideList[index + j]
                        if guessPeptide[0:len(suffix)] == suffix:
                            if guessPeptide[-1]!='1':
                                return index + j
                            else:
                                j -= 1
                        else:
                            break
                    # if neither loop yields a match, we will not find this split and return False.
                    return None

            # if the suffix didn't match, simply continue.
            guessPair = [suffix, guessPeptide]
            if sorted(guessPair)[0] == suffix:
                newRange = (range[0], index - 1)
            else:
                newRange = (index, range[1])
            return self.findSuff(suffix, newRange)

    def updatePepList(self):
        newList = []
        for peptide in self.peptideList:
            if peptide[-1]!='1':
                newList.append(peptide)
        self.peptideList = newList
        self.pepListLenOld = self.pepListLen
        self.pepListLen = len(self.peptideList)

    def concatRemaining(self):
        noSeqPerConcat = math.ceil(self.pepListLen/NO_RECORDS)
        newSeq = ""
        finalSeq = []
        for i in range(0, self.pepListLen):
            if i % noSeqPerConcat == 0:
                finalSeq.append(newSeq)
                newSeq = ""
            newSeq += self.peptideList[i]
        self.peptideList = finalSeq

    def createOutput(self):
        self.overlapList()
        print("concat loop finished")
        while(self.pepListLen > NO_RECORDS*2):
            self.updatePepList()
            if self.pepListLenOld == self.pepListLen:
                self.concatRemaining()
                break
            print("new length: " + str(self.pepListLen))
            self.overlapList()
            print("concat loop finished")


def createPepList(OUTPUT_PATH):
    pepList = []
    with open(OUTPUT_PATH, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            pepList.append(str(record.seq))
        pepList.sort()
    return pepList

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

def concatPepsFromFile():
    pepList = createPepList(OUTPUT_PATH)
    concatListObject = ConcatList(pepList)
    concatListObject.createOutput()
    # write to file
    with open('concatOutput.fasta', "w") as output_handle:
        SeqIO.write(createSeqObj(concatListObject.peptideList), output_handle, "fasta")
    #checkOutput(OUTPUT_PATH, 'concatOutput.fasta')

def concatPepsFromSet(pepSet, outputPath):
    pepList = list(pepSet)
    pepList.sort()
    concatListObject = ConcatList(pepList)
    concatListObject.createOutput()
    # write to file
    with open(outputPath, "w") as output_handle:
        SeqIO.write(createSeqObj(concatListObject.peptideList), output_handle, "fasta")

# old find suff function
def findSuffOld(suffix, peptideList, zeroIndex):
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
            return findSuffOld(suffix, smallerPeptideList, zeroIndex)
        else:
            smallerPeptideList = peptideList[index:]
            return findSuffOld(suffix, smallerPeptideList, index+zeroIndex)

concatPepsFromFile()
#checkOutput(OUTPUT_PATH, "concatOutput.fasta")

