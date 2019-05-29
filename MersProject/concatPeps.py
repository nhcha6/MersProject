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

        # call self.overlapList to initiate the first round of concatenation if the length of the peptide list is less
        # than N)_RECORDS.
        if self.pepListLen > NO_RECORDS:
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

        # if the peptide has a 1 on the end, it has already been concatenated this iteration and should be skipped.
        if peptide[-1]=='1':
            #print('woooo')
            return

        # as concated peptides get longer we do not want to iterate through all suffixes. Max suffix length is set
        # to be 20
        startIter = len(peptide) - 20
        if startIter < 1:
            startIter = 1

        # loop through the different suffixes in order from longest to shortest
        for j in range(startIter, len(peptide)):
            # extract suffix from peptide
            suffix = peptide[j:]
            # extract index of peptide with matching prefix from the list. Return None if the suffix cannot be found.
            prefixIndex = self.findSuff(suffix, (0,self.pepListLen-1))

            # if a prefix was found, reconfigure self.peptideList and return from the function.
            if prefixIndex != None:
                # store the prefixPeptide
                prefixPeptide = self.peptideList[prefixIndex]
                # add a 1 to the end of the prefix peptide so that it is skipped this iteration, and deleted from
                # self.peptideList before the next concat iteration.
                self.peptideList[prefixIndex] = prefixPeptide + '1'
                # create the concatenated peptide.
                overlapPeptide = concatOverlapPep(peptide, j, prefixPeptide)
                # replace the current peptide with the new, overlapping peptide.
                self.peptideList[i] = overlapPeptide
                return

    def findSuff(self, suffix, range):
        """
        This recursive function is initially called by self.createOvelap(). It is essentially based on the binary
        search algorith, taking a suffix and searching through a sorted list for a matching prefix. It does so by
        guessing the middle element of the input range, and then halving the range and calling the function again.
        When a match is found the index is returned, and if the final prefix does not match it returns None.

        :param suffix: the sequence that is being searched for as a prefix in one of the peptides.
        :param range: the range of self.peptideList that is to be assessed. This is halved each time the function is
        recursively called as per binary search.

        :return: None if the suffix cannot be found, the index of the matched prefix within self.peptideList if it
        is found.
        """
        # check if the range only encompasses one peptide.
        if range[0] == range[1]:
            # store the index being guessed and the guessed peptide.
            index = range[0]
            guessPeptide = self.peptideList[index]
            # if the guessPeptide has a prefix matching suffix, and guessPeptide hasn't already been combined (detailed
            # by the addition of a 1 to the end) we return the index.
            if guessPeptide[0:len(suffix)] == suffix and guessPeptide[-1]!='1':
                return index
            # if either of the above conditions fail, return None.
            else:
                return None
        # if there are at least two peptides in the range, continue as normal.
        else:
            # calculate the guess index and guessPeptide
            index = math.ceil((range[1] - range[0]) / 2) + range[0]
            guessPeptide = self.peptideList[index]

            # if the guessPeptide contains the suffix, continue to next check.
            if guessPeptide[0:len(suffix)] == suffix:
                # if the guess peptide doesn't have the 1 flag, return the index as we have a match.
                if guessPeptide[-1]!='1':
                    return index
                # however if there is the 1 flag, check the peptides either side to see if a match that hasn't been
                # flagged exists.
                else:
                    # iterate forward through self.peptideList and try find a match without the 1 on the end.
                    j = 1
                    while True:
                        # guess the peptide after the one that matched.
                        guessPeptide = self.peptideList[index+j]
                        # if it matches return index+j
                        if guessPeptide[0:len(suffix)]==suffix:
                            if guessPeptide[-1]!='1':
                                return index+j
                            # if it fails based on another 1 being there, index j and try the next peptide.
                            else:
                                j+=1
                        # if the prefix doesn't match the suffix, break the while loop.
                        else:
                            break

                    # iterate backwards in the same manner if the forward iteration failed.
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

                    # if neither loop yields a match, we will not find this split so we return False.
                    return None

            # if the suffix didn't match, we reduce the size of the range and call self.findSuff again.
            guessPair = [suffix, guessPeptide]
            # check if the suffix should sit to the left or the right of the guessPeptide.
            if sorted(guessPair)[0] == suffix:
                # if it sits to the left, take the lower half of the current range.
                newRange = (range[0], index - 1)
            else:
                # if it sits to the right, take the upper half.
                newRange = (index, range[1])
            # call self.findSuff on the new range and same suffix and return the ultimate return value.
            return self.findSuff(suffix, newRange)

    def updatePepList(self):
        """
        Called by self.createOutput() after self.overlapList() runs an iteration of concatenation. It iterates through
        self.peptideList and deletes all peptides which have been flagged with a 1 on the end of them. It then updates
        the self.pepListLen and self.pepListLenOld variables.
        :return:
        """
        newList = []
        for peptide in self.peptideList:
            if peptide[-1]!='1':
                newList.append(peptide)
        self.peptideList = newList
        self.pepListLenOld = self.pepListLen
        self.pepListLen = len(self.peptideList)

    def concatRemaining(self):
        """
        Called by self.createOverlap() when self.pepListLen and self.pepListLenOld are equal after an iteration of
        concatenation, and thus combining based on suffix/prefix matching is no longer having any affect. This
        function reduces the number of sequences in self.peptideList to under the NO_RECORDS by arbitrarily
        concatenating them.
        :return:
        """
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

def createPepList(OUTPUT_PATH):
    """
    Called by concatPepsFromFile() when this script is being used independently of the splicing program. This function
    simply stores the peptides in the input file path in a list, and sorts them alphabetically.

    :param OUTPUT_PATH: the a input fasta file of peptides.
    :return pepList: the peptides in the input file sorted alphabetically in a list.
    """
    pepList = []
    with open(OUTPUT_PATH, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            pepList.append(str(record.seq))
        pepList.sort()
    return pepList

def concatOverlapPep(peptide, j, prefixPeptide):
    """
    Called by self.createOverlap(), this function takes two peptides which have an identical prefix/suffix sequence
    and combines them around this like sequence. Eg: ABCDE + DEFGH = ABCDEFGH

    :param peptide: the peptide with matching suffix sequence
    :param j: the length of the matching suffix
    :param prefixPeptide: the peptide with matching prefix seqeunce.

    :return concatPep: the peptide resulting from concatenation around the matching prefix/suffix.
    """
    concatPep = peptide[0:j] + prefixPeptide
    return concatPep

def createSeqObj(seenPeptides):
    """
    Called by concatPepsFromFile() and concatPepsFromSet, this function takes a list of peptides and configures them
    for writing to fasta. Given the list of matchedPeptides, converts all of them into SeqRecord objects and passes
    back a generator.

    :param seenPeptides: the peptides to be written to file.
    :return:
    """
    count = 1
    seqRecords = []

    for sequence in seenPeptides:
        if not sequence.isalpha():
            continue
        finalId = "ipd|pep" + str(count) + ';'
        yield SeqRecord(Seq(sequence), id=finalId, description="")
        count += 1

def concatPepsFromFile():
    """
    Called when this script is to be run independently of the splicing program. It initiates the concatenation of all
    the peptides in the file path declared globally as OUTPUT_PATH to reduce the number of sequences to under the number
    set by NO_RECORDS.
    :return:
    """
    pepList = createPepList(OUTPUT_PATH)
    concatListObject = ConcatList(pepList)
    concatListObject.createOutput()
    # write to file
    with open('concatOutput.fasta', "w") as output_handle:
        SeqIO.write(createSeqObj(concatListObject.peptideList), output_handle, "fasta")
    #checkOutput(OUTPUT_PATH, 'concatOutput.fasta')

def concatPepsFromSet(pepSet, outputPath):
    """
    Called from removeSubsets.py when the concatenated output has been selected in the splicing program. It takes a set
    of peptides and an outputPath and initiates the concatenation of all the peptides to reduce the number of sequences
    to under the number set by NO_RECORDS.

    :param pepSet: a set of peptides to be concatenated.
    :param outputPath: the output path that the concatenated peptides are to be written to.
    :return:
    """
    pepList = list(pepSet)
    pepList.sort()
    concatListObject = ConcatList(pepList)
    concatListObject.createOutput()
    # write to file
    with open(outputPath, "w") as output_handle:
        SeqIO.write(createSeqObj(concatListObject.peptideList), output_handle, "fasta")

#concatPepsFromFile()
#checkOutput(OUTPUT_PATH, "concatOutput.fasta")


