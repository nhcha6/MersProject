from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import os
import math

# define the input file containing peptides that you wish to concatenate. If this script is being run via the splicing
# program, this variable will not be needed.
OUTPUT_PATH = 'example6-14_Linear_1_040419_0925_NoSubsets.fasta' #'C:/Users/Administrator/Desktop/Remove Subseqs/a2Maxmods3-Linear050219_2324_NoSubsets.fasta'

# the max number of sequences you want the input peptides to be concatenated to. The actual number of output peptides
# will likely be slightly less.
NO_RECORDS = 4000

# define peptides you want to check are passing through the code correctly. They are directed to a separate function
# which enables targeted error checking.
# **can be deleted from final code.
CHECK_PEPS = []#['QHAAAAAAAAAA', 'RKEEAAAAAAAAAAA', 'KGGAAAAAAAAAA']

class ConcatList:
    """
        Class that handles the concatenation of peptide sequences that have shared prefixes and suffixes. It allows
        a single list of peptides to shared across a series of functions.

        Methods: self.overlapList(), self.createOverlap(), self.findSuff(), self.updatePepList(),
                 self.concatRemaining(), self.createOutput().

        Class Variables:
            self.peptideList = a list of peptides which are to be concatenated. This list is reduced in a series of
                               concatenations until it is under the NO_RECORDS.
            self.pepListLen = the length of the current self.peptideList.
            self.pepListLenOld = the length of the previous self.peptideList. When this equals self.pepList,
                                 concatenation using prefix/suffix is no longer having an impact, and the code
                                 turns to direct concatenation.

    """

    def __init__(self, pepList):
        """
        A list of peptides is required to create a ConcatList object.

        :param pepList: a list of peptides.
        """
        self.peptideList = pepList
        self.pepListLen = len(pepList)
        self.pepListLenOld = 0

    def createOutput(self):
        """
        Called by concatPepsFromFile() or concatPepsFromSet(). This function controls a series of concatenation cycles,
        and updates the peptide list after each cycle. It also checks if self.pepListLen == self.pepListLenOld and
        concatenates the remaining peptides without looking for peptides with similar prefixes/suffixes.

        :return:
        """

        # call self.overlapList to initiate the first round of concatenation
        self.overlapList()
        print("concat loop finished")

        # while self.pepListLen is more than twice the prescribed NO_RECORDS, we want to run another round of
        # concatenation.
        while(self.pepListLen > NO_RECORDS*2):
            # update the peptide list (this consists of removing peptides which have been flagged for deletion by
            # inserting a 1 on the end.
            self.updatePepList()
            # check that the prefix/suffix concatenation technique is still concatenating peptides.
            if self.pepListLenOld == self.pepListLen:
                # if it has stopped working, concatenate the remaining peptides arbitrarily and finish the program.
                self.concatRemaining()
                break
            # otherwise, print the new length of self.peptideList and run another round of concatenation.
            print("new length: " + str(self.pepListLen))
            self.overlapList()
            print("concat loop finished")

    def overlapList(self):
        """
        Called by self.createOutput(), this function controls a single round of concatenation by finding peptides
        with similar prefix and suffix sequences. It indexes i and passed it into self.createOverlap(), allowing
        self.createOverlap() to move through each peptide in self.peptideList and concatenate it with another peptide.
        :return:
        """
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
        """
        Called by overlapList(), this function takes an index which references a peptide in self.peptideList. This
        peptide is extracted, and each potential suffix of it is iterated through from longest to shortest.
        When a suffix is found as a prefix in another peptide the two peptides are concatenated together. This
        concatenated peptide replace the original peptide at self.peptideList[i], and the second peptide is combined
        with gets a 1 added to the end of it to flag it for later deletion.

        :param i: the index location of the peptide to be concatenated.
        :return:
        """
        # store peptide and check that a number hasn't been added to the front of it.
        peptide = self.peptideList[i]

        for missing in CHECK_PEPS:
            if missing in peptide:
                self.createOverlap2(i)

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
            if guessPeptide[0:len(suffix)] == suffix and guessPeptide[-1]!='1':
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
        print(self.pepListLen)
        print(noSeqPerConcat)
        newSeq = ""
        finalSeq = []
        for i in range(0, self.pepListLen):
            if (i+1) % noSeqPerConcat == 0:
                print(i)
                finalSeq.append(newSeq)
                newSeq = ""
            newSeq += self.peptideList[i]
        finalSeq.append(newSeq)
        self.peptideList = finalSeq

    # **delete these functions once output is confirmed to be a ok.
    def createOverlap2(self,i):
        # store peptide and check that a number hasn't been added to the front of it.
        peptide = self.peptideList[i]

        print("\nPeptide: " + peptide + '\n')

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
            prefixIndex = self.findSuff2(suffix, (0,self.pepListLen-1))
            print("\npref index: " + str(prefixIndex) + '\n')

            if prefixIndex != None:
                # store the prefixPeptide, and replace its location in the list with None
                prefixPeptide = self.peptideList[prefixIndex]
                print('prefixPeptide: ' + prefixPeptide)
                self.peptideList[prefixIndex] = prefixPeptide + '1'
                overlapPeptide = concatOverlapPep(peptide, j, prefixPeptide)
                print("conatPeptide " + overlapPeptide)
                # replace the current peptide with the new, overlapping peptide.
                self.peptideList[i] = overlapPeptide
                return

    # need to rewrite this function so that we do not create a subset of the peptideList each time. This is the only
    # place i can find that is causing inefficiencies.
    # find suff is essentially a binary search algorithm.
    def findSuff2(self, suffix, range):
        print("SUffix: " + suffix)
        print("range: " + str(range))
        if range[0] == range[1]:
            index = range[0]
            guessPeptide = self.peptideList[index]
            print('guess peptide: ' + guessPeptide)
            # check that the peptide we are up to is not None. If it is None, we want to return False as there are no
            # peptides left to iterate through
            if guessPeptide[0:len(suffix)] == suffix and guessPeptide[-1]!='1':
                return index
            else:
                return None
        else:
            index = math.ceil((range[1] - range[0]) / 2) + range[0]
            guessPeptide = self.peptideList[index]
            print("guess peptide: " + guessPeptide)
            if guessPeptide[0:len(suffix)] == suffix:
                if guessPeptide[-1]!='1':
                    return index
                else:
                    print("matched to a 1, checking either side.")
                    # iterate forward through self.peptideList and try find a match without the 1 on the end.
                    j = 1
                    while True:
                        guessPeptide = self.peptideList[index+j]
                        print("guess peptide: " + guessPeptide)
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
                        print("guess peptide: " + guessPeptide)
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
            print("guess pair: " + str(sorted(guessPair)))
            if sorted(guessPair)[0] == suffix:
                newRange = (range[0], index - 1)
            else:
                newRange = (index, range[1])
            return self.findSuff2(suffix, newRange)

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

# ** old find suff function and checkOutput function to be deleted
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

#concatPepsFromFile()
#checkOutput(OUTPUT_PATH, "concatOutput.fasta")

