from concatPeps import *

# input file name and sortedFile name. This is used if running this script independantly of the rest of the splicing
# program.
PEPTIDE_FASTA = 'example6-14_Linear_1_040419_0925'

def removeSubsetSeq(outputPath):
    """
    This function manages the flow of the script, calling the functions required to remove from the input file
    all peptides which exist as subsets in longer peptides.

    :param outputPath: the input fasta file from which subset peptides are to be deleted.
    :return:
    """
    # edit the outputPath to create the desired inputPath, sortedPath, noSubseqPath and concatPath
    inputPath = outputPath + ".fasta"
    sortedPath = outputPath + "_Sorted.fasta"
    noSubseqPath = outputPath + "_NoSubsets.fasta"
    concatPath = outputPath + "_concat.fasta"

    # multiple timing variables used to time separate parts of the code
    time1 = time.time()

    # create list of all peptides recorded in the file
    seenPeptides, origNo = seenPepList(inputPath)

    time2 = time.time()

    print('Number of original entries: ' + str(origNo))
    print('Time to read all entries to list: ' + str(time2-time1))

    # sort list of peptides and write to new file
    sortList(sortedPath, seenPeptides)

    time3 = time.time()
    print('All sequences sorted and written to new dictionary in: ' + str(time3-time2))

    # convert seenPeptides to a set
    seenPeptides = set(seenPeptides)

    time4 = time.time()
    print('Second read of input file took: ' + str(time4-time3))

    # open sorted fasta to iterate through the sequences in order from largest to smallest. We do so because
    # a large peptide will not be a subset of a smaller one, and thus we can delete peptides sooner and reduce
    # the runtime of the algorithm by starting with the largest.
    with open(sortedPath, "rU") as handle:
        seenPeptides = pepRemoveNoOrigin(handle, seenPeptides)

    # remove sorted.fasta from where it is saved
    os.remove(sortedPath)

    print('Time to delete subset sequences: ' + str(time.time()-time4))
    print('No. of sequences reduced from ' + str(origNo) + ' to ' + str(len(seenPeptides)))

    # send list of peptides without subsets to be concatenated by the script concatPeps.py. If we want to also print
    # this list of peps, simply uncomment the SeqIO.write line.
    with open(noSubseqPath, "w") as output_handle:
        #SeqIO.write(createSeqObj(seenPeptides), output_handle, "fasta")
        concatPepsFromSet(seenPeptides, concatPath)

def createSeqObj(seenPeptides):
    """
    Called by removeSubsetSeq(). Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a function which
    writes them to file.

    :param seenPeptides: the list of peptides that is to be written to file
    :return:
    """
    count = 1
    seqRecords = []

    for sequence in seenPeptides:
        finalId = "ipd|pep" + str(count) + ';'
        yield SeqRecord(Seq(sequence), id=finalId, description="")
        count += 1

def seenPepList(filePath):
    """
    Called by Called by removeSubsetSeq(). From the file path input, this function adds all peptides to a list and
    returns the list along with its size.

    :param filePath: the path of the input peptide file.

    :return seenPeptides: a list of all the peptide sequences in the input file.
    :return origNo: the number of peptides in the list.
    """
    seenPeptides = []
    # iterate through the input peptide file and add all sequences to a list
    with open(filePath, "rU") as handle:
        origNo = 0
        for record in SeqIO.parse(handle, 'fasta'):
            origNo += 1
            seenPeptides.append(str(record.seq))
    return seenPeptides, origNo

def sortList(savePath, seenPeptides):
    """
    Called by removeSubsetSeq(), this function takes an output file location and a list of peptides, and write the
    peptides to file in order of longest to shortest.

    :param savePath: the location that the sorted peptides are to be written to.
    :param seenPeptides: the unsorted list of peptides to be written to sorted file.
    :return:
    """
    seenPeptides.sort(key=len)
    seenPeptides.reverse()
    # write sorted list to fasta file so it need not be stored in memory
    with open(savePath, "w") as output_handle:
        SeqIO.write(createSeqObj(seenPeptides), output_handle, "fasta")

def pepRemoveNoOrigin(handle, seenPeptides):
    """
    Called by removeSubsetSeq(), this function iterates through each peptide in the sorted peptideList (from longest
    to shortest) and deletes the peptides from seenPeptides which are subsets of the current peptide.

    :param handle: the filePath of the fasta file containing the sorted peptides.
    :param seenPeptides: the set of all peptides from the initial input file.

    :return seenPeptides: the edited set of peptides which has all subsets removed.
    """
    # iterate through each record in the sorted.fasta
    counter = 1
    for record in SeqIO.parse(handle, 'fasta'):
        # print the counter ever 10000 records to track progress of code without slowing it by printing too much.
        if counter % 10000 == 0:
            print(counter)
        counter += 1
        pep = str(record.seq)
        # if the current peptide has already been deleted from seenPeptides because it is a subset of a longer
        # peptide, we simply move on as any subsets of this peptide will have been deleted also.
        if pep not in seenPeptides:
            continue
        # for peptides which are not subsets of any others, we create all the possible subset sequences (same as
        # linear splits) in order to check if they match any other sequence in seenPeptides
        for i in range(len(pep)):
            for j in range(i + 1, len(pep) + 1):
                # if a given split/subsequences is in seenPeptides it is deleted from seenPeptides. Additional check
                # is to ensure that we do not delete the original peptide which has been split up from seenPeptides
                if pep[i:j] in seenPeptides and pep[i:j] is not pep:
                    seenPeptides.remove(pep[i:j])
    return seenPeptides

# call removeSubsetSeq to run the script in isolation on the file set by outputPath. If this script is being run
# in conjunction with the splicing program, it should be commented out.
#removeSubsetSeq(PEPTIDE_FASTA)