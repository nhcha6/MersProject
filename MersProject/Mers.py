from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from TransPlaceholder import *
import csv
from MonoAminoAndMods import *
import multiprocessing
from multiprocessing import Queue
from collections import Counter
import time
import sys
# import h5py
import json
import logging
from MGFMain import *
import atexit
import os
import psutil
import tempfile
from queue import Queue
import io
import traceback
from pathlib import Path
import math
from removeSubsets import *

# define our the spliceType flags
TRANS = "Trans"
LINEAR = "Linear"
CIS = "Cis"

# MEMORY_THRESHOLD is the percentage of RAM being used at which all computation is paused, the current output data is
# written to file, after which the output recommences.
MEMORY_THRESHOLD = 85
# NUM_PROC_TOTAL is the total number of processes generated per trans/cis/linear splicing computation.
NUM_PROC_TOTAL = 1000
# The generation of processes is slowed to avoid an overload of memory due to spawned but uncompleted processes.
# MAX_PROC_ALIVE defines the maximum number of processes that have been started (stored in Fasta.procGenCounter) and
# not yet completed (as counted by Fasta.completedProcs).
MAX_PROC_ALIVE = 40

# define the flags which are passed between the the writer() function and genMassDict/transProcess to communicate
# details of the computation process.
MEMFLAG = 'mem'
STOPFLAG = 'stop'
PROC_FINISHED = 'Process Finished'

# logging module to be used as an alternative to print statements.
logging.basicConfig(level=logging.DEBUG, format='%(message)s')
# logging.disable(logging.INFO)

# initialise mgfData as none.
mgfData = None

# massDict syntax: PEPTIDE: [monoisotopic mass, [referenceLocation], {charge: m/z}]
# for example: {'DQQ': [389.15466499999997, ['121', '117-118'], {2: 195.58527249999997}],
# 'QDQ': [389.15466499999997, ['118-119', '117'], {2: 195.58527249999997}]......


class Fasta:

    """
    Class that represents the input from a fasta file.
    Methods: generateOutput, transOutput, cisAndLinearOutput
    Class Variables:
        self.inputFile = a list of the all fasta files uploaded containing input proteins.
        self.procGenCounter = a counter which stores the total number of processes generated.
        self.pepCompleted = a multiprocessing queue which is added to every time a the output from a process is pulled
                            from the toWriteQueue in the writer() function.
        self.completedProcs = a counter which is indexed everytime a process is finished, which the writer() function
                              and self.pepCompleted keep track of.
        self.totalProcs = the number of processes generated in total. This is estimated based on NUM_PROC_TOTAL, and
                          adjusted if the number of input proteins per process is calculated to be 1.
    """

    def __init__(self, inputFile):
        """
        A Fasta object simply requires the input path to a fasta file for it to be created.

        :param inputFile: the input path to a fasta file.
        """

        self.inputFile = [inputFile]
        self.procGenCounter = 0
        self.pepCompleted = multiprocessing.Queue()
        self.completedProcs = 0
        self.totalProcs = 0

    def generateOutput(self, mined, maxed, overlapFlag, concatFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                       protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfObj, mgfFlag):

        """
        This function is called from MersGUI.py to initiate the output process once the user input has been recieved
        and confirmed suitable.

        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
        cis spliced peptides.
        :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.
        :param transFlag: if True, trans splicing should be completed as part of the output.
        :param cisFlag: if True, cis splicing should be completed as part of the output.
        :param linearFlag: if True, linear splicing should be completed as part of the output.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath.
        :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :param maxMod: the max number of modifications allowable per peptide.
        :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
        that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
        :param outputPath: a dictionary holding the output file name of a trans, cis and linear splice output. The
        keys are the TRANS, LINEAR and CIS flags defined at the top of this script.
        :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :param mgfObj: the object that contains the mgf data uploaded in the format required. This is processed before
        generate ouptut is reached.
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.
        :return:
        """
        # estimatethe total number of processes that will be generated (NUM_PROC_TOTAL for every splice type running)
        if transFlag:
            self.totalProcs += NUM_PROC_TOTAL
        if cisFlag:
            self.totalProcs += NUM_PROC_TOTAL
        if linearFlag:
            self.totalProcs += NUM_PROC_TOTAL

        # sequentially check if each splice type has been selected by the user, and run the relevant function which
        # will generate processes for this splice type.
        if transFlag:
            self.transOutput(self.inputFile, mined, maxed, modList, maxMod,outputPath[TRANS], chargeFlags, mgfObj,
                        mgfFlag, csvFlag, pepToProtFlag,protToPepFlag, concatFlag)

        if cisFlag:
            self.cisAndLinearOutput(self.inputFile, CIS, mined, maxed,overlapFlag, csvFlag, pepToProtFlag,
                                    protToPepFlag, modList, maxMod, maxDistance, outputPath[CIS], chargeFlags, mgfObj,
                                    mgfFlag, concatFlag)

        if linearFlag:
            self.cisAndLinearOutput(self.inputFile, LINEAR, mined, maxed, overlapFlag, csvFlag, pepToProtFlag,
                                    protToPepFlag, modList, maxMod, maxDistance, outputPath[LINEAR], chargeFlags,
                                    mgfObj, mgfFlag, concatFlag)


    def transOutput(self, inputFile, mined, maxed, modList, maxMod, outputPath, chargeFlags, mgfObj, mgfFlag, csvFlag,
                    pepToProtFlag, protToPepFlag, concatFlag):
        """
        This function is called Fasta.generateOutput() if transFlag == True. It controls the creation of all processes
        required to conduct trans splicing.

        :param inputFile: the fasta files containing proteins input by the user.
        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :param maxMod: the max number of modifications allowable per peptide.
        :param outputPath: the file name of the output fasta file.
        :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :param mgfObj: the object that contains the mgf data uploaded in the format required. This is processed before
        generate ouptut is reached.
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath
        :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.
        :return:
        """

        # initialise the final path
        finalPath = None

        # Open the csv file if the csv file is selected
        if csvFlag:
            finalPath = getFinalPath(outputPath)
            open(finalPath, 'w')

        # initialise the dictionary to store the protein sequences as values and the protein names as keys
        seqDict = {}

        # check number of input files and pass flag into addSequenceList so that the seqId can updated to include file
        # name if there are more than 1 files imported.
        if len(inputFile) > 1:
            multiFileFlag = True
        else:
            multiFileFlag = False

        # loop through files and add to seqDict
        for file in inputFile:
            seqDict.update(addSequenceList(file, multiFileFlag))

        # combine all the proteins into one long sequence, stored in finalPeptide (** rename to finalProtein??).
        # We also need protIndexList to store the range of each initial protein in this longer protein. Stored as a
        # list of start and end references: [[0,150],[151,200]....]
        # protList stores the initial proteins in the same order as the references in protIndexList. That is, the
        # references at index i of of protIndexList corresponds to the location of protein at reference i in protList.
        finalProtein, protIndexList, protList = combinePeptides(seqDict)

        # create all the splits and their reference with respect to finalPeptide
        splits, splitRef = splitTransPeptide(finalProtein, mined, maxed, protIndexList)

        # define the length of splits for later calculations
        splitLen = len(splits)

        # define the number of workers to be input into the multiprocessing.Pool() based on the number of CPUs in the
        # computer.
        num_workers = multiprocessing.cpu_count()

        # Used to lock write access to file when writing from separate processes
        lockVar = multiprocessing.Lock()

        # initialise the queue which will feed the outputs generated in each process to the writer process.
        toWriteQueue = multiprocessing.Queue()
        # initialise the queue which will feed all cis and linear peptides dynamically created to the writer queue.
        # we do not want any cis or linear spliced peptides to appear in the final trans splice output.
        linCisQueue = multiprocessing.Queue()

        # declare the pool.
        pool = multiprocessing.Pool(processes=num_workers, initializer=processLockTrans, initargs=(lockVar, toWriteQueue,
                                                                                                   splits,splitRef, mgfObj,
                                                                                                   modTable, linCisQueue))

        # declare and start the writer process.
        writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, linCisQueue, pepToProtFlag,
                                                                     protToPepFlag, self.pepCompleted, concatFlag, True))
        writerProcess.start()

        # Create a process for groups of splits, pairing element 0 and 1 with -1 and -2 and so on. The indexes of the
        # splits to be computed in this fashion is stored in splitsIndex.
        # the size of each process (2 splits from the front and 2 from the end in the above example) is set by
        # dividing the number of splits by NUM_PROC_TOTAL*2.
        splitsIndex = []
        procSize = math.ceil(splitLen / (NUM_PROC_TOTAL*2))

        # if procSize = 1, NUM_PROC_TOTAL may be much larger than the number of procs actually generated. Thus, to
        # ensure the progress bar still works readjust self.totalProcs to be splitLen/2 instead.
        if procSize == 1:
            self.totalProcs = self.totalProcs - NUM_PROC_TOTAL + math.ceil(splitLen/2)

        # these nested create a list of indexes which denote a single process. This list is stored in splitsIndex.
        # splitsIndex is altered each iteration.
        # continuing with the example where procSize = 2: splitsIndex = [0,1,-1,-2] and then [2,3,-3-4] and so on.
        for i in range(0, math.ceil(splitLen / 2), procSize):
            if i + procSize > math.floor(splitLen / 2):
                for j in range(i, splitLen - i):
                    splitsIndex.append(j)
            else:
                for j in range(i, i + procSize):
                    splitsIndex.append(j)
                    splitsIndex.append(splitLen - 1 - j)

            # wait until the number processes started but not completed is less than MAX_PROC_ALIVE to begin a new
            # process.
            self.procGenCounter += 1
            while True:
                if self.procGenCounter - self.completedProcs < MAX_PROC_ALIVE:
                    break

            # once the above while loop is broken, we can start a new process.
            pool.apply_async(transProcess, args=(splitsIndex, mined, maxed, modList, maxMod,
                                                 finalPath, chargeFlags, mgfFlag, csvFlag, protIndexList, protList))
            # reset splitsIndex to [] so that the next iteration can begin.
            splitsIndex = []

            # after a process has been completed, check that the memory being used has not exceded MEMORY_THRESHOLD.
            # if the memory limit has been reached, close the pool and put MEMFLAG to toWriteQueue to tell the
            # writer() to output all the data that is stored in memory. We then wait for the writing to complete,
            # restart the pool and continue creating the remaining processes.
            if memory_usage_psutil() > MEMORY_THRESHOLD:
                print('Memory usage exceded. Waiting for processes to finish.')
                pool.close()
                pool.join()
                toWriteQueue.put(MEMFLAG)
                comp = self.completedProcs
                # wait for writer to communicate back that it is done by adding one to self.completedProcs
                while self.completedProcs == comp:
                    continue
                self.completedProcs += -1
                print('Restarting Pool')
                pool = multiprocessing.Pool(processes=num_workers, initializer=processLockTrans,
                                            initargs=(lockVar, toWriteQueue, splits, splitRef, mgfObj, modTable,
                                                      linCisQueue))

        # wait for all the processes to finish before continuing.
        pool.close()
        pool.join()

        # send the flag to the toWriteQueue to tell the writer() function that the process is complete.
        toWriteQueue.put(STOPFLAG)
        # wait for the writer process to finish computing.
        writerProcess.join()
        # print that the trans computation is completed.
        logging.info("All Trans joined")

    def cisAndLinearOutput(self, inputFile, spliceType, mined, maxed, overlapFlag, csvFlag, pepToProtFlag, protToPepFlag,
                           modList, maxMod, maxDistance, outputPath, chargeFlags, mgfObj, mgfFlag, concatFlag):

        """
        This function is called Fasta.generateOutput() if cisFlag == True and if linFlag == True. It controls
        the creation of all processes required to conduct trans splicing.

        :param inputFile: the fasta files containing proteins input by the user.
        :param spliceType: holds either cis or linear depending on what splice type is being run.
        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
        cis spliced peptides.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath
        :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
        which can be found in the file MonoAminoAndMods.py.
        :param maxMod: the max number of modifications allowable per peptide.
        :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
        that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
        :param outputPath: the file name of the output fasta file.
        :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
        mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
        :param mgfObj: the object that contains the mgf data uploaded in the format required. This is processed before
        generate ouptut is reached.
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.
        :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
        into approximately 6000 sequences.

        :return:
        """

        # declare final path.
        finalPath = None

        # Open the csv file if the csv file is selected
        if csvFlag:
            finalPath = getFinalPath(outputPath)
            open(finalPath, 'w')

        # declare the num_workers, which is to be used when declaring the pool, based on the number of cores available
        # in the computer.
        num_workers = multiprocessing.cpu_count()

        # Used to lock write access to files
        lockVar = multiprocessing.Lock()

        # declare toWriteQueue to feed the output of individual processes to the writer queue for processing and
        # printing
        toWriteQueue = multiprocessing.Queue()
        # declare linSetQueue to feed linear spliced peptides created during the cis process to the writer queue so
        # that they can be deleted from the final output.
        linSetQueue = multiprocessing.Queue()

        # start the pool and the writerProcess.
        pool = multiprocessing.Pool(processes=num_workers, initializer=processLockInit,
                                    initargs=(lockVar, toWriteQueue, mgfObj, modTable, linSetQueue))
        writerProcess = multiprocessing.Process(target=writer,
                                                args=(toWriteQueue, outputPath, linSetQueue, pepToProtFlag,
                                                      protToPepFlag, self.pepCompleted, concatFlag))
        writerProcess.start()

        # calculate total size of input fasta by iterating through each record.
        totalProt = 0
        for file in inputFile:
            with open(file, "rU") as handle:
                for entry in SeqIO.parse(handle, 'fasta'):
                    totalProt += 1

        # calculate the number of proteins that need too be in each process to approximately create the number of
        # processes set by NUM_PROC_TOTAL.
        # ** change to protPerProc
        pepPerProc = math.ceil(totalProt / NUM_PROC_TOTAL)
        # if pepPerProc == 1, it is likely there are less proteins than NUM_PROC_TOTAL. Therefore, we adjust
        # the self.totalProcs counter (used in the progress bar).
        if pepPerProc == 1:
            self.totalProcs = self.totalProcs - NUM_PROC_TOTAL + totalProt

        # initialise the counter of the number of proteins iterated through, and the protDict which stores the proteins
        # to be put into each process.
        counter = 0
        protDict = {}
        # iterate through each record in each file.
        for file in inputFile:
            with open(file, "rU") as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    # add to the protein counter and extract the seq and seqId from the record.
                    counter += 1
                    seq = str(record.seq)
                    seqId = record.name

                    seqId = seqId.split('|')[1]
                    seqId = seqId.split(';')[0]
                    # add the filename to the seqId if more than one file has been included
                    if len(inputFile) > 1:
                        fileName = file.split('/')[-1]
                        seqId = fileName.split('.')[0] + '_' + seqId

                    # add seqId and seq to protDict
                    protDict[seqId] = seq

                    # if counter % pepPerProc == 0, we are ready to create a process which will compute the cis and/or
                    # linear peptides which can be generated by the proteins in protDict.
                    if counter % pepPerProc == 0:
                        # add 1 to self.procGenCounter which will be used to control the creation of processes.
                        self.procGenCounter += 1
                        # the code is stalled in the while loop until the number of processes that have been created
                        # but not had their output retrieved from the toWriteQueue is less tha MAX_PROC_ALIVE.
                        while True:
                            if self.procGenCounter - self.completedProcs < MAX_PROC_ALIVE:
                                break
                        # Once the while loop is passed, start the processes for this set of processes with the
                        # target function being genMassDict
                        pool.apply_async(genMassDict, args=(spliceType, protDict, mined, maxed, overlapFlag,
                                                            csvFlag, modList, maxMod, maxDistance, finalPath,
                                                            chargeFlags, mgfFlag))
                        # reset protDict to be empty, ready for the next set of proteins for the next process.
                        protDict = {}

                        # check if the percentage of RAM used is approaching the upper limit set by MEMORY_THRESHOLD
                        # if the threshold has been exceded, close the pool and wait for all processes to finish via
                        # join(). Then put the MEMFLAG to the toWriteQueue, which instructs the writer() function
                        # to write all the computed peptides to file.
                        if memory_usage_psutil() > MEMORY_THRESHOLD:
                            print('Memory usage exceded. Waiting for processes to finish.')
                            pool.close()
                            pool.join()
                            toWriteQueue.put(MEMFLAG)
                            comp = self.completedProcs
                            # wait for writer to communicate back that it is done by adding 1 to self.completedProcs.
                            # when we see self.completedProcs has increased, we are ready to restart the pool (after
                            # we have decremented self.completedProcs).
                            while self.completedProcs == comp:
                                continue
                            self.completedProcs += -1
                            print('Restarting Pool')
                            pool = multiprocessing.Pool(processes=num_workers, initializer=processLockInit,
                                                        initargs=(lockVar, toWriteQueue, mgfObj, modTable,
                                                                  linSetQueue))

        # once all records are iterated through, a process is created for the remaining proteins in protDict.
        if protDict:
            self.procGenCounter += 1
            pool.apply_async(genMassDict, args=(spliceType, protDict, mined, maxed, overlapFlag,
                                                csvFlag, modList, maxMod, maxDistance, finalPath,
                                                chargeFlags, mgfFlag))

        # wait for all processes to finish before continuing.
        pool.close()
        pool.join()

        # send STOPFLAG to the writer queue to communicate that all processes have been completed.
        toWriteQueue.put(STOPFLAG)
        # wait for the writerQueue to finish via the .join() command.
        writerProcess.join()
        logging.info("All " + spliceType + " !joined")

def combinePeptides(seqDict):

    """
    This function is called by Fasta.transOutput(). It takes a dictionary with protein names as keys and proteins sequences as values and concatenates the sequences
    to form one sequence. It returns a list of the initial protein names and a list of references denoting where
    the initial proteins appear in the long, concatenated protein.

    :param seqDict: a dictionary with protein names as keys and proteins sequences as values.
    :return finalPeptide: a single sequence of all the individual protein sequences concatenated together.
    :return protIndexList: a list of index pairs denoting the start and end position of each individual protein.
    List is structures as follows: [[0,150], [151,205] ... ]
    :return protList: a list of the protein names of the individual proteins, where the index of each name corresponds
    to the index of its location data in protIndexList.
    """

    # declare variables to be created within for loop.
    dictlist = []
    protIndexList = []
    protList = []
    # initial ind is 0, and is updated by the length of the protein added each iteration.
    ind = 0
    for key, value in seqDict.items():
        dictlist.append(value)
        protIndexList.append([ind,ind + len(value) - 1])
        protList.append(key)
        ind += len(value)

    # combine all protein sequences in dictList to form a long, concatenated sequence.
    finalPeptide = ''.join(dictlist)
    return finalPeptide, protIndexList, protList

def splitTransPeptide(proteinSeq, mined, maxed, protIndexList):

    """
    Called from Fasta.transOutput(), this function creates all the possible splits (peptide cleavages) from the input protein sequence. This differs from
    splitDictPeptide due to functionality required specifically by trans. The input protein to this function
    is all the inividual input proteins concatenated together, creating potential splits across the border between
    two proteins which cannot be formed given the input data. The protIndexList data is thus included to ensure that
    such cleavages are ignored from the fonal splits output.

    :param proteinSeq: the input proteins sequence, which will be all individual proteins concatenated together.
    :param mined: the min length of a potential output peptide.
    :param maxed: the max length of a potential output peptide.
    :param protIndexList: a list of location pairs denoting the start and end point of original, individual proteins.
    :return splits: all the potential cleavages which could be recombined to create a trans splice peptide.
    :return splitRef: where each split originated in the input protein.
    """

    # Makes it easier to integrate with earlier iteration where linearFlag was being passed as an external flag
    # instead of spliceType
    linearFlag = False
    length = len(proteinSeq)

    # splits will hold all possible splits that can occur
    splits = []
    # splitRef will hold a direct reference to the characters making up each split string: for starting peptide ABC,
    # the split AC = [0,2]
    splitRef = []

    # embedded for loops build all possible splits
    for i in range(0, length):

        character = proteinSeq[i]
        toAdd = ""

        # figure out which protein the splits starts in, and the max index the splits can reach before it becomes
        # a part of a second peptide
        initProt, protInd = findInitProt(i, protIndexList)

        # add and append first character and add and append reference number which indexes this character
        toAdd += character
        ref = list([i+1])
        temp = list(ref)  # use list because otherwise shared memory overwrites

        # check first character is a valid amino acid.
        if toAdd in monoAminoMass.keys():
            splits.append(toAdd)
            splitRef.append(temp)
        else:
            continue

        # iterates through every character after current and adds it to the most recent string if max size
        # requirement is satisfied
        for j in range(i + 1, length):
            if j > initProt[1]:
                break
            toAdd += proteinSeq[j]
            if linearFlag:
                ref.append(j+1)
                if maxSize(toAdd, maxed):
                    if minSize(toAdd, mined):

                        # ensure all characters in toAdd are valid amino acids.
                        if not aminoCheck(toAdd):
                            break
                        splits.append(toAdd)
                        temp = list(ref)
                        splitRef.append(temp)
                else:
                    break

            else:
                if maxSize(toAdd, maxed-1):
                    # ensure all characters in toAdd are valid amino acids.
                    if not aminoCheck(toAdd):
                        break
                    splits.append(toAdd)
                    ref.append(j+1)
                    temp = list(ref)
                    splitRef.append(temp)
                else:
                    break

    return splits, splitRef

def transProcess(splitsIndex, mined, maxed, modList, maxMod, finalPath,
                 chargeFlags, mgfFlag, csvFlag, protIndexList, protList):
    """
    Called as the worker process to the mulitprocessing.Pool() created in Fasta.transOutput(). This function controls
    the flow of computing the trans output of a given subset of splits. It calls functions to create
    the trans splice peptides, calculate the relevant mass and charge data, apply modifications and conduct
    required mgf comparison and put the output the toWriteQueue so that the writer queue can process it.

    :param splitsIndex: stores the indexes of the splits which this process will compute the output for.
    :param mined: the minimum length of an output peptide.
    :param maxed: the maximum length of an ouptut peptide.
    :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
    which can be found in the file MonoAminoAndMods.py.
    :param maxMod: the max number of modifications allowable per peptide.
    :param finalPath: the specific output path of the trans output.
    :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
    mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
    :param mgfFlag: if True, the user has selected not to compare to any MGF data.
    :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
    or precursor mass comparison is conducted.
    :param protIndexList: a list of index pairs denoting the start and end position of each individual protein.
    List is structures as follows: [[0,150], [151,205] ... ]
    :param protList: a list of the protein names of the individual proteins, where the index of each name corresponds
    to the index of its location data in protIndexList.
    :return:
    """

    # try and except are to ensure that an error message is reported if something doesn't work within an individual
    # process.
    try:
        # Look to produce only trans spliced peptides - not linear or cis. Do so by not allowing combination of peptides
        # which originate from the same protein as opposed to solving for Cis and Linear and not including that
        # in the output
        combined, combinedRef, linCisSet = combineTransPeptide(splits, splitRef, mined, maxed, splitsIndex, protIndexList)

        # Put linCisSet to linCisQueue:
        transProcess.linCisQueue.put(linCisSet)

        # update combineRef to include information on where the peptide originated from
        origProtTups = findOrigProt(combinedRef, protIndexList, protList)

        # Convert it into a dictionary that has a mass
        massDict = combMass(combined, combinedRef, origProtTups)

        # Apply mods to the dictionary values and update the dictionary
        massDict = applyMods(massDict, modList, maxMod)

        # Add the charge information along with their masses
        massDict = chargeIonMass(massDict, chargeFlags)

        # Get the positions in range form, instead of individuals (0,1,2) -> (0-2)
        massDict = editRefMassDict(massDict)

        if mgfFlag:
            allPeptides = massDict.keys()
            allPeptidesDict = {}
            for peptide in allPeptides:
                # create the string, with peptides sorted so all permutations are matched as similar. There may be multiple
                # peptide locations in the list of tuples, hence the for loop. Tuples are listed in order, with consecutive
                # tuples relating to a pair of splice locations.
                string = ""
                # create seenOrigins list to store origins which have been added so that there is no double up of identical origins.
                seenOrigins = []
                # iterate through all the origins stored in the massDict for the given key
                for i in range(0, len(massDict[peptide][3]), 2):
                    # sort the origin so that it can be compared to others accurately
                    origProt = sorted(massDict[peptide][3][i:i + 2])
                    # if origProt has already been seen (and in turn added to seenOrigins) continue iterating.
                    if origProt in seenOrigins:
                        continue
                    # if origProt hasn't been seen, append it to seenOrigins and add it to the string which is to be added to allPeptidesDict.
                    seenOrigins.append(origProt)
                    string += origProt[0][0] + origProt[0][1] + '/' + origProt[1][0] + origProt[1][1] + ';'
                string = string[0:-1]
                allPeptidesDict[peptide] = string
            transProcess.toWriteQueue.put((allPeptidesDict,False))

        # If there is an mgf file AND there is a charge selected
        elif mgfData is not None and True in chargeFlags:
            matchedPeptides, modCountDict = generateMGFList(TRANS, mgfData, massDict, modList)
            transProcess.toWriteQueue.put((matchedPeptides, modCountDict))

        # If csv is selected, write to csv file
        if csvFlag:
            logging.info("Writing locked :(")
            lock.acquire()
            writeToCsv(massDict, TRANS, finalPath, chargeFlags)
            lock.release()
            logging.info("Writing released!")


    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in worker process: ' + str(splitsIndex) + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e

def combineTransPeptide(splits, splitRef, mined, maxed, splitsIndex, protIndexList):

    """
    Called from transProcess(). Takes splits index from the multiprocessing pool and computes all the possible peptides which can be created
    from these splits. The entire splits and splitRef are global variables within the pool and can thus be accessed
    by all processes. This function differs from the lin/cis version as it must check that a created peptide is not
    a linear or cis peptide, and that it must only iterate through the splits denoted by splits index, as opposed to
    all of them.

    ** do all the splits need to be passed into the function? Surely this is consuming a lot of unneccesary memory
    when they are redeclared.

    :param splits: all cleaved peptides which can be recombined to form a trans spliced peptide.
    :param splitRef: the locations of these cleavages within the concatenated input protein.
    :param mined: the minimum length of an output peptide.
    :param maxed: the maximum length of an ouptut peptide.
    :param splitsIndex: the indexes of the splits to be combined with all other splits in this given process.
    :param protIndexList: a list of index pairs denoting the start and end position of each individual protein.
    List is structures as follows: [[0,150], [151,205] ... ]
    :return combModless: all the trans peptides (combined cleavage/splits) created in this process.
    :return combModlessRef: the location data of where this trans peptide originated in the concatenated protein
    sequence
    :return linCisSet: a set of all the linear and cis spliced peptides that were dynamically produced during
    the combination process. These will be deleted from the final output later.
    """
    # initialise linCisVariable holder.
    linCisSet = set()
    # initialise combinations array to hold the possible combinations from the input splits
    combModless = []
    combModlessRef = []

    # iterate through all of the splits and build up combinations which meet min/max/overlap criteria
    for i in splitsIndex:
        # toAdd variables hold temporary combinations for insertion in final matrix if it meets criteria
        toAddForward = ""

        toAddReverse = ""

        for j in range(i, len(splits)):
            # create forward combination of i and j
            toAddForward += splits[i]
            toAddForward += splits[j]
            addForwardRef = splitRef[i] + splitRef[j]
            toAddReverse += splits[j]
            toAddReverse += splits[i]
            addReverseRef = splitRef[j] + splitRef[i]

            # max, min and max distance checks combined into one function for clarity for clarity
            if combineCheck(toAddForward, mined, maxed, splitRef[i], splitRef[j], 'None'):
                # overlap is forced in trans hence why the flag has not been included.
                # linCisPepCheck assesses if the splits used to create the peptide were from the same protein.
                # If so, its a cis or lin protein and not trans, and we add it to linCisSet not trans
                if linCisPepCheck(addForwardRef, protIndexList):
                    linCisSet.add(toAddForward)
                    linCisSet.add(toAddReverse)
                # if the splits are from different proteins, it is a trans peptide and can be added.
                else:
                    combModless.append(toAddForward)
                    combModlessRef.append(addForwardRef)
                    combModless.append(toAddReverse)
                    combModlessRef.append(addReverseRef)

            toAddForward = ""
            toAddReverse = ""

    return combModless, combModlessRef, linCisSet

def findOrigProt(combinedRef, protIndexList, protList):
    """
    Called from transProcess(). Iterates through each combineRef and returns the origin protein names and the index within each protein for
    every protein.

    :param combinedRef: the location data of where this trans peptide originated in the concatenated protein
    sequence. Of the form: [[1,2,3,4,5,101,102,103], [6,7,8,9,10,110,111,112,113]....]
    :param protIndexList: a list of indexes, which store the start and end value within the concatenated protein
    of each input protein.
    :param protList: a list of the names of the original input proteins. Ordered to match protIndexList.
    :return proteinTups: a list of pairs of tuples that stores the name of the original protein a peptide was from,
    and if the sub-sequence is more than 6 amino acids in length, it stores the location also. This list is also
    ordered to match the indexes of protList andprotIndexList.
    Format: proteinTups = [(prot1, "(1-6)"), (prot2,""), (prot3, "(1-6)"), (prot4, "(10-16)")....]
    """
    proteinTups = []
    for i in range(0, len(combinedRef)):
        # declare the protRefs for this iteration
        protRef1 = ""
        protRef2 = ""
        # declare the ref for this iteration
        ref = combinedRef[i]
        # call findInitProt to return the index of the protein that the first split is from (protIter1), and the relevant reference
        # to the location of this protein within the overall concatenated protein (protIndex1).
        protIndex1, protIter1 = findInitProt(ref[0] - 1, protIndexList)
        # return the origin protein via the index protIter1
        prot1 = protList[protIter1]

        # special check if peptide is overlap spliced. If so, let protRef1 equal the subset of peptides that this
        # splicing is from, and append the tuple pair [(prot1, protRef1),('Overlap',"")] to proteinTups.
        if len(set(ref)) != len(ref):
            protRef1 += ('(' + str(min(ref) - protIndex1[0]))
            protRef1 += ('-' + str(max(ref) - protIndex1[0]) + ')')
            proteinTups.append([(prot1, protRef1),('Overlap',"")])
            continue

        # iterate through all the values in ref to split it into two independent references relating to the two
        # cleavages which were combined to make the peptide.
        for j in range(1,len(ref)):
            # check if the current ref is outside the range of the origin protein the first split was found to
            # belong to.
            if ref[j] - 1 > protIndex1[1] or ref[j] - 1 < protIndex1[0]:
                # check to see if the first split is at least 6 amino acids in length.
                # if so append the location of the split within the peptide to prot1
                if j > 5:
                    protRef1 += ('(' + str(ref[0] - protIndex1[0]))
                    protRef1 += ('-' + str(ref[j-1] - protIndex1[0]) + ')')

                # if we have entered the loop, we are at a ref which belongs to the second split.
                # We know need to find the index of the protein that this split is from (protIter2), and the relevant reference
                # to the location of this protein within the overall concatenated protein (protIndex2).
                protIndex2, protIter2 = findInitProt(ref[j] - 1, protIndexList)
                prot2 = protList[protIter2]
                # same as above, check if second split is at least 6 amino acids long, and if so update protRef2
                # with the location data.
                if len(ref) - j > 5:
                    protRef2 += ('(' + str(ref[j] - protIndex2[0]))
                    protRef2 += ('-' + str(ref[-1] - protIndex2[0]) + ')')

                # append the pair of protein tups to proteinTups and break the inner loop.
                proteinTups.append([(prot1,protRef1),(prot2,protRef2)])
                break
    return proteinTups

def findInitProt(index, protIndexList):
    """
    Called from a number of functions, but most notably from findOrigProt(). Takes an index to an amino acid
    within the concatenated input protein, and locates which range within protIndexList this reference
    belongs to. It then return this range of values (relating to the range the origin protein occupies
    in the concatenated sequence) and the index of this range within protIndexList which can be used to return the
    name of the original protein given by this range.

    :param index: the input location that we are trying to locate within protIndexList
    :param protIndexList: a list of index pairs denoting the start and end position of each individual protein.
    List is structures as follows: [[0,150], [151,205] ... ]
    :return protIndexList[protIter]: the specific range within the concatenated protein sequence that the
    input reference is contained within.
    :return protIter: the index of the located range within protIndexList. This index is used to extract the name
    of the origin protein which corresponds to the range.
    """
    # return the last ref of the last range to acquire the length.
    length = protIndexList[-1][-1]
    # Find the average length. Plus 1 needed for when Protein length is perfectly divisible by protein index length
    aveLen = math.ceil(length/len(protIndexList)) + 1
    # the starting iter that the input index will be at is predicted using the average length and assigned to
    # protIter.
    protIter = math.floor(index/aveLen)
    # protIter may be predicted to be the outside the range of indexes (max in can be outside by is 1).
    # If so reduce it by 1 to bring it within the range.
    if protIter == len(protIndexList):
        protIter -= 1
    # check if the index is within the range at the initially predicted protIter. If it is, return the index and range.
    # if it is not, either increase or decrease protIter depending on if the input index was higher or lower than
    # the current prediction.
    while True:
        lower = protIndexList[protIter][0]
        upper = protIndexList[protIter][1]
        if lower <= index:
            if upper >= index:
                return protIndexList[protIter], protIter
            else:
                protIter += 1
        else:
            protIter -= 1

def genMassDict(spliceType, protDict, mined, maxed, overlapFlag, csvFlag, modList, maxMod,
                maxDistance, finalPath, chargeFlags, mgfFlag):

    """
    Called as the worker function to the multiprocessing.Pool() in Fasta.cisAndLinearOutput(). This is the worker
    function which controls the computation and output of a given process for cis and linear splicing. It calls
    functions to create the cis or linear peptides, calculates relevant masses and m/z ratios, applies modifications
    and compares this data to an MGF file is desire.

    :param spliceType: holds either the CIS or TRANS flag so that the appropriate splicing can be computed.
    :param protDict:
    :param mined: the minimum length of an output peptide.
    :param maxed: the maximum length of an ouptut peptide.
    :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
    permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
    cis spliced peptides.
    :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
    or precursor mass comparison is conducted.
    :param modList: a list of the modifications input by the user. The modifications match the keys in modTable,
    which can be found in the file MonoAminoAndMods.py.
    :param maxMod: the max number of modifications allowable per peptide.
    :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
    that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.
    :param finalPath: the specific output path/file name for the cis or linear output file
    :param chargeFlags: a list of bools which defines which charge flags are to be assessed when comparing to the
    mgf file. The first bool corresponds to +1, the second to +2 and so on up to +5.
    :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
    is to be ouptut to Fasta.

    :return:
    """

    # try and except included so that miscellaneous bugs in the code can be picked up within processes generated by
    # the pool.
    try:
        # iterate through each entry in the protDict.
        for protId, protSeq in protDict.items():

            # Get the initial peptides and their positions, and the set of linear peptides produced for this protein
            # if conducting cis splicing.
            combined, combinedRef, linSet = outputCreate(spliceType, protSeq, mined, maxed, overlapFlag, maxDistance)

            # add this set of linear proteins to the linSetQueue
            genMassDict.linSetQueue.put(linSet)

            # Convert combined and combinedRef into a dictionary which has a mass.
            # data structure: peptideSequence: [monoisotopic mass, [referenceLocation]]
            massDict = combMass(combined, combinedRef)

            # Apply mods to the dictionary values and update the dictionary
            massDict = applyMods(massDict, modList, maxMod)

            # Add m/z ratios for each selected charge to the massDict.
            massDict = chargeIonMass(massDict, chargeFlags)

            # Get the combinedRef data in range form, instead of individuals (0,1,2) -> (0-2)
            massDict = editRefMassDict(massDict)

            # if memFlag is True, the user has selected not to conduct any mgf comparison. Thus, we simply want to
            # extract the peptide sequences from massDict and add them to the tooWriteQueue. modCountDict is set
            # to False as it is irrelevant given no mgf comparison has occured.
            # ** move this to after linsetQueue and break genMassDict after completed, as we do not need to create
            # massDict at all if mgf comparison isn't occuring.
            if mgfFlag:
                #allPeptides = getAllPep(massDict)
                allPeptides = massDict.keys()
                allPeptidesDict = {}
                for peptide in allPeptides:
                    allPeptidesDict[peptide] = protId
                genMassDict.toWriteQueue.put((allPeptidesDict,False))

            # If there is an mgf file AND there is a charge selected, call functions to compare the potential spliced
            # proteins in massDict to the spectra in the MGF file.
            elif mgfData is not None and True in chargeFlags:
                # matched peptides is all the peptides which passed the MGF comparison. modCountDict contains data
                # on the number of times certain mods were included in the output which passed. This information
                # is otherwise lost as all modded peptides which pass MGF comparison are returned to their unmodded
                # form before being added to matchedPeptides.
                matchedPeptides, modCountDict = generateMGFList(protId, mgfData, massDict, modList)
                genMassDict.toWriteQueue.put((matchedPeptides, modCountDict))

            # If csv is selected, write to csv file
            if csvFlag:
                logging.info("Writing locked :(")
                lock.acquire()
                writeToCsv(massDict, protId, finalPath, chargeFlags)
                lock.release()
                logging.info("Writing released!")

        # put PROC_FINISHED flag to toWriteQueue, so that it can update Fasta.completedProcs.
        genMassDict.toWriteQueue.put(PROC_FINISHED)

    # if an error is raised, it is caught by this exception and its details are printed to console.
    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in protein: ' + protId + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e

# set default maxDistance to be absurdly high incase of trans
def outputCreate(spliceType, protein, mined, maxed, overlapFlag, maxDistance=10000000):
    """
    Called from genMassDict(). Recieves a protein and creates the peptides using the relevant input spliceType and in accrodance with the
    min length, max length, overlapFlag and maxDistance input by the user.

    :param spliceType: either the CIS or LINEAR flag, which denotes which type of splicing to conduct.
    :param protein: the input protein sequence.
    :param mined: the minimum number of amino acids in a spliced peptide.
    :param maxed: the maximum number of amino acids in a spliced peptide.
    :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
    permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
    cis spliced peptides.
    :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
    that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.

    :return combined: all the spliced peptides created from the input protein which fit the other criteria.
    :return combinedRef: a list containing where the amino-acids which were used to create a given spliced peptide
    where located within the input protein. Eg: [[0], [0,1], [0,2], [1], [1,2].....]
    :return linSet: if cis splicing is being run, this set will contain all the linear spliced peptides that
    could be produced from the input protein in accordance with the input criteria. If linear splicing is being run,
    this will simply be empty.
    """

    # Splits eg: ['A', 'AB', 'AD', 'B', 'BD']
    # SplitRef eg: [[0], [0,1], [0,2], [1], [1,2]]
    # Produces splits and splitRef arrays which are passed through combined
    splits, splitRef = splitDictPeptide(spliceType, protein, mined, maxed)
    combined, combinedRef = None, None

    if spliceType == CIS:
        # combined eg: ['ABC', 'BCA', 'ACD', 'DCA']
        # combinedRef eg: [[0,1,2], [1,0,2], [0,2,3], [3,2,0]]
        # pass splits through combinedOverlapPeptide() if cis splicing has been selected
        combined, combinedRef, linSet = combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance)

    elif spliceType == LINEAR:
        # if linear splicing, the splits are simply the linear spliced peptides. Set them to equal combined/combinedRef
        # for clear visibility of this. We also want linSet to contain an empty set.
        combined, combinedRef = splits, splitRef
        linSet = set()

    return combined, combinedRef, linSet

def splitDictPeptide(spliceType, protein, mined, maxed):

    """
    Called from outputCreate(). This function creates all the potential cleavages of the input protein which meet
    the input mined, maxed and spliceType parameters. If spliceType == CIS, we want all the cleavages from length 1
    to length maxed-1 inclusive. If sliceType == LINEAR, we want all the cleavage from length mined to length maxed
    inclusive.

    :param spliceType: either CIS or LINEAR, denotes which type of splicing is being run.
    :param protein: the input protein sequence from which the cleavages are generated.
    :param mined: the minimum number of amino acids in a spliced peptide.
    :param maxed: the maximum number of amino acids in a spliced peptide.

    :return splits: a list storing all possible cleavages created from the input sequence in accordance with the
    input criteria. If running linear splicing, this list simply holds all linear spliced peptides.
    :return splitRef: a list containing where each split was located within the the input protein.
    Eg: [[0], [0,1], [0,2], [1], [1,2].....]
    """

    # Makes it easier to integrate with earlier iteration where linearFlag was being passed as an external flag
    # instead of spliceType
    linearFlag = spliceType == LINEAR
    length = len(protein)

    # splits will hold all possible splits that can occur
    splits = []
    # splitRef will hold a direct reference to the characters making up each split string: for starting peptide ABC,
    # the split AC = [0,2]
    splitRef = []

    # embedded for loops build all possible splits
    for i in range(0, length):

        character = protein[i]
        toAdd = ""
        # add and append first character and add and append reference number which indexes this character

        toAdd += character
        ref = list([i+1])
        temp = list(ref)  # use list because otherwise shared memory overwrites

        # linear flag to ensure min is correct for cis and trans
        if linearFlag and minSize(toAdd, mined):
            # Don't need to continue this run as first amino acid is unknown X
            if not toAdd in monoAminoMass.keys():
                continue
            else:
                splits.append(toAdd)
                splitRef.append(temp)

        elif not linearFlag:
            if toAdd in monoAminoMass.keys():
                splits.append(toAdd)
                splitRef.append(temp)

        # iterates through every character after current and adds it to the most recent string if max size
        # requirement is satisfied
        for j in range(i + 1, length):
            toAdd += protein[j]
            if linearFlag:
                ref.append(j+1)
                if maxSize(toAdd, maxed):
                    if minSize(toAdd, mined):
                        # Check that the split to be added contains only valid amino acids. Break inner loop if so.
                        if not aminoCheck(toAdd):
                            break
                        splits.append(toAdd)
                        temp = list(ref)
                        splitRef.append(temp)
                else:
                    break

            else:
                if maxSize(toAdd, maxed-1):
                    # Check that the split to be added contains only valid amino acids. Break inner loop if so.
                    if not aminoCheck(toAdd):
                        print(toAdd)
                        break
                    splits.append(toAdd)
                    ref.append(j+1)
                    temp = list(ref)
                    splitRef.append(temp)
                else:
                    break

    return splits, splitRef

def aminoCheck(split):
    """
    Called by splitDictPeptide() and splitTransPeptide(), this function takes a cleavage from a peptide sequence
    and checks all characters in it are valid amino acids.

    :param split: the cleavage which is to be checked for validity.
    :return: True if the split contains only amino acid characters, False if not.
    """
    for char in split:
        if char in monoAminoMass.keys():
            continue
        else:
            return False
    return True

def combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance):

    """
    Called by outputCreate(). This function takes a list of splits and the corresponding location data within the
    input protein to create all the cis spliced peptides which fit within the other input criteria.
    :param splits: a list storing all possible cleavages created from the input sequence in accordance with the
    input criteria.
    :param splitRef: a list containing where each split was located within the the input protein.
    Eg: [[0], [0,1], [0,2], [1], [1,2].....]
    :param mined: the minimum number of amino acids in a spliced peptide.
    :param maxed: the maximum number of amino acids in a spliced peptide.
    :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
    permitted. Shared amino acids originate from the same amino-acid in the same peptide. This is only relevant to
    cis spliced peptides.
    :param maxDistance: in cis splicing, the maximum distance between any two amino acids in two cleaved peptides
    that are to be recombined to form a cis spliced peptide. If None, the max distance is infinite.

    :return combined: all the cis spliced peptides created from the the splits which fit the input criteria.
    :return combinedRef: a list containing where the amino-acids which were used to create a given spliced peptide
    where located within the input protein. Eg: [[0], [0,1], [0,2], [1], [1,2].....]
    :return linSet: all the linear spliced peptides that could be produced from the input splits in accordance
    with the input criteria.
    """
    # initialise linSet to store the linear spliced peptides which can be created from this set of splits.
    linSet = set()
    # initialise massDict to hold the possible combinations from the input splits. Use a dictionary so that duplicates
    # are automatically removed.
    massDict = {}
    # initialise the output arrays combModless and comModlessRef.
    combModless = []
    combModlessRef = []
    # iterate through all of the splits and build up combinations which meet min/max/overlap criteria
    for i in range(0, len(splits)):

        # toAdd variables hold temporary combinations for insertion in final matrix if it meets criteria
        toAddForward = ""
        toAddReverse = ""

        # iterate through all splits which appear after the current split so that all combinations of splits can
        # be created.
        for j in range(i, len(splits)):
            # create forward combination of split i and j
            toAddForward += splits[i]
            toAddForward += splits[j]
            addForwardRef = splitRef[i] + splitRef[j]
            # create backwards combinations of split i and j
            toAddReverse += splits[j]
            toAddReverse += splits[i]
            addReverseRef = splitRef[j] + splitRef[i]

            # max, min and max distance checks combined into one function for clarity. If it doesn't pass, continue to
            # the next combination of splits.
            if combineCheck(toAddForward, mined, maxed, splitRef[i], splitRef[j], maxDistance):
                # if overlap flag is True, we do need to check if the two splits overlap each other, and
                # not allow combination if they do.
                if overlapFlag:
                    # returns False is the two splitRef lists overlap.
                    if overlapComp(splitRef[i], splitRef[j]):
                        # add toAddReverse to massDict
                        massDict[toAddReverse] = addReverseRef
                        # check if toAddForward is linear and add to linearSet if so
                        if linCisPepCheck(addForwardRef, False):
                            linSet.add(toAddForward)
                        # if toAddForward is not linear, add it to massDict.
                        else:
                            massDict[toAddForward] = addForwardRef

                # if overlap flag is False, we are happy to combine the splits even if they overlap.
                else:
                    # add toAddReverse to massDict
                    massDict[toAddReverse] = addReverseRef
                    # check if toAddForward is linear and add to linearSet if so
                    if linCisPepCheck(addForwardRef, False):
                        linSet.add(toAddForward)
                    else:
                        massDict[toAddForward] = addForwardRef
            # if combineCheck() returns False, check that the massDistance value isn't what caused it to fail.
            # if it was maxDistance, we can break from the inner for loop, as the following splits are only going
            # to be further away.
            elif not maxDistCheck(splitRef[i], splitRef[j], maxDistance):
                break

            # reset toAddForward/Reverse at the end of each iteration through the inner for loop, ready for the
            # next combination of splits.
            toAddForward = ""
            toAddReverse = ""

    # iterate through all items massDict and add the peptide to combModless and reference to combModlessRef if the
    # peptide is not in linSet.
    for peptide, ref in massDict.items():
        if peptide in linSet:
            continue
        else:
            combModless.append(peptide)
            combModlessRef.append(ref)

    return combModless, combModlessRef, linSet


def writer(queue, outputPath, linCisQueue, pepToProtFlag, protToPepFlag, procCompleted, concatFlag, transFlag=False):
    """
    This function is called as the worker to the writerProcess, which is generated for each individual splice type
    from either Fasta.cisAndLinearOutput() or Fasta.transOutput(). It accesses the output data put the toWriteQueue
    by each individual processes and sorts it into a dictionary called seenPeptides. It also accesses the linCisPeptides
    put to the linCisQueue by the individual process and stores this data in a set so that these peptides can be
    removed from seenPeptides before they are output.
    The process also deals with writing the output fasta file once all processes have finished, or once all the
    memory threshold is reached.

    :param queue: the queue (called toWriteQueue in the pool worker functions) to which the processes put their
    output peptide dictionaries once they have completed.
    :param outputPath: the file path name that the output fasta file is to be written to.
    :param linCisQueue: the queue to which all the processes put their linear (if running cis) or linear and cis (if
    running trans) peptides to so that they can be collected and deleted from the output peptides.
    :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
    they originated in listed underneath.
    :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
    which they produced listed underneath.
    :param procCompleted: the queue which each processes puts to once it has finished calculating all the input
    proteins which were allocated to it. It is used to control the progress bar and to stall process generation
    so that a backlog of unfinished processes in the writer queue doesn't clog memory.
    :param concatFlag: if True, an additional file will be output containing the output peptides concatenated to
    into approximately 6000 sequences.
    :param transFlag: True if the writer process is dealing with a trans splicing output. This allows the trans origin
    data, which is distinct from cis and linear, to be formatted sufficiently for output.
    :return:
    """

    # initialise relevant variables
    seenPeptides = {}
    linCisSet = set()
    saveHandle = str(outputPath)
    modCountDict = Counter()
    fileCount = 0
    memThreshFlag = False

    # run on repeat to get from the queue continuously.
    while True:
        # get from cisLinQueue and from matchedPeptide Queue
        matchedTuple = queue.get()
        if not linCisQueue.empty():
            linCisSet = linCisSet | linCisQueue.get()

        # if stop is sent to the matchedPeptide Queue, everything has been output, so we exit the while loop and
        # begin configuring the output in seenPeptides to be written to file.
        if matchedTuple == STOPFLAG:
            logging.info("Everything computed, stop message has been sent")
            break

        # If MEMFLAG is sent, we know the max memory has been hit during process generation. Thus we want to write
        # everything in seenPeptides to file and then communicate back to the main thread that it can continue
        # generating processes.
        if matchedTuple == MEMFLAG:
            # update memThreshFlag so that linCisPeptides can be removed from the file at the end
            memThreshFlag = True
            # remove linear/cis peptides from seenPeptides:
            print('memflag recieved, writing temp file')
            commonPeptides = linCisSet.intersection(seenPeptides.keys())
            for peptide in commonPeptides:
                del seenPeptides[peptide]
            print('deleted common peptides')
            # write to ouptut file.
            fileCount += 1
            writeOutputFiles(seenPeptides, protToPepFlag, pepToProtFlag, transFlag, outputPath, fileCount)
            print('written to output')
            seenPeptides = {}
            # put back to procCompleted queue to tell process generating function that it can restart the queue.
            procCompleted.put(1)
            continue

        # PROC_FINISHED flag gets put via writer queue everytime a process is finished. When recieved, add one to
        # the procCompleted queue so that counter in the main thread can be indexed.
        if matchedTuple == PROC_FINISHED:
            procCompleted.put(1)
            continue
        # if trans, each queue.get that reaches here (ie isn't a flag of some-sort) corresponds to a process
        # thus, we need to add one to procCompleted queue also.
        if transFlag:
            procCompleted.put(1)

        # if metchedTuple is not MEMFLAG, STOPFLAG or PROC_FINISHED, it is a genuine output and we continue as
        # normal to add it to seenPeptides.
        matchedPeptides = matchedTuple[0]
        # if matchedTuple[1] isn't False, it contains relevant data and we need to update modCountDict.
        if matchedTuple[1]:
            modCountDict += matchedTuple[1]

        # Add the matchedPeptides from the given process to seenPeptides.
        for key, value in matchedPeptides.items():
            # convert origins to a list, the delimitter of ; is only relevant for trans.
            origins = value.split(';')
            # update seenPeptides appropriately according to if the peptide already exists and if any of the origin
            # data already exists.
            if key not in seenPeptides.keys():
                seenPeptides[key] = origins
            else:
                for origin in origins:
                    if origin not in seenPeptides[key]:
                        seenPeptides[key].append(origin)

    # once the while loop as been broken, check if seenPeptides contains any elements. This is only to ensure the
    # memory threshold wasn't hit on the last process, and thus all seenPeptides would already have been written
    # to file.
    if seenPeptides:
        # Was not over memory threshold but last items need to be dealt with.
        commonPeptides = linCisSet.intersection(seenPeptides.keys())
        for peptide in commonPeptides:
            del seenPeptides[peptide]

        # write all seenPeptides to output.
        fileCount += 1
        writeOutputFiles(seenPeptides, protToPepFlag, pepToProtFlag, transFlag, outputPath, fileCount)

    # if memThreshFlag is True, we have written a tempFile. We must therefore go back and check that all peptides
    # in linCisSet haven't made their way into the these outputs.
    if memThreshFlag:
        # check that there is anything in linCisSet to compare to (for linear splicing it will be empty)
        if linCisSet:
            # we now want to read each output file into memory, delete peptides in linCisQueue, and write it
            # to file again.
            remFinalCisLin(linCisSet, saveHandle, fileCount)

    if concatFlag:
        concatOutput(saveHandle, fileCount)

    # if modCountDict contains meaningful data, we need to add it to the info file.
    if modCountDict:
        # need to know if related to cis/lin/trans. Replace the relevant portion of the saveHandle with info to
        # get the save handle of the info file.
        if '_Linear' in saveHandle:
            title = 'LINEAR MODIFICATION COUNT' + '\n'
            infoPath = saveHandle.replace("_Linear", "_Info")
            infoPath = infoPath[0:-6] + '.txt'
        elif '_Cis' in saveHandle:
            title = 'CIS MODIFICATION COUNT' + '\n'
            infoPath = saveHandle.replace("_Cis", "_Info")
            infoPath = infoPath[0:-6] + '.txt'
        else:
            title = 'TRANS MODIFICATION COUNT' + '\n'
            infoPath = saveHandle.replace("_Trans", "_Info")
            infoPath = infoPath[0:-6] + '.txt'

        # write the modCountDict data to the info file.
        file = open(infoPath, 'a')
        file.write('\n' + title)
        for key, value in modCountDict.items():
            file.write(key + ': ' + str(value) + '\n')

def getAllPep(massDict):
    """
    Takes all peptides in massDict, converts them to their unmodified equivalent and then adds them to the list
    allPeptides.
    :param massDict: a dictionary containing spliced peptides as keys and a value containing mono-isotopic mass,
    m/z ratios for different charge states and the origin of that peptide within the input protein.
    :return allPeptides: a list of all peptides which were generated, without any modifications.
    """
    allPeptides = set()
    for key, value in massDict.items():
        if not key.isalpha():
            alphaKey = modToPeptide(key)
        else:
            alphaKey = key
        allPeptides.add(alphaKey)
    return allPeptides

def memory_usage_psutil():
    """
    Calculates and returns the percentage of computer's RAM being used.
    :return mem.percent: the percentage of the computer's RAM being used.
    """
    # return the memory usage in percentage like top
    mem = psutil.virtual_memory()

    return mem.percent

def remFinalCisLin(linCisSet, saveHandle, fileCount):
    """
    Called by writer(), this function is run at the end of all output files being written if the memory threshold
    was hit and thus there is more than 1 file. It iterates through each output file excluding the last one, and
    deletes peptides that are also in the linCisSet. Such files were written before the linCisSet was complete, and
    thus not all had been deleted in the original writing.
    :param linCisSet: the final set of linear peptides (for cis splice) or linear and cis peptides (for trans splice)
    which are to be deleted from the output files.
    :param saveHandle: the standard saveHandle of all outputs which is edited to add the output file number.
    :param fileCount: the number of files created in total.
    :return:
    """
    # iterate through each file, excluding the most recently written file.
    for i in range(1, fileCount):
        # declare the path of the file for this iteration
        finalPath = str(saveHandle)[0:-17] + '_' + str(i) + '_' + str(saveHandle)[-17:]
        # open the file to read its contents into the new dictionary.
        with open(finalPath, 'rU') as handle:
            # initialise the dictionary to store the peptides for each file.
            outputDict = {}
            # two counters: one to count the peptide we are up to, one to count the peptide we are adding
            oldPepCount = 0
            newPepCount = 0
            # iterate through each record
            for record in SeqIO.parse(handle, 'fasta'):
                oldPepCount += 1
                peptide = str(record.seq)
                name = str(record.name)
                # if this peptide is in linCisSet, continue to the next peptide. Otherwise add it to outputDict.
                if peptide in linCisSet:
                    continue
                else:
                    newPepCount += 1
                    finalId = name.replace(str(oldPepCount), str(newPepCount), 1)
                    outputDict[peptide] = finalId

        # open the same file path and replace its contents with the
        with open(finalPath, 'w') as handle:
            for peptide, finalId in outputDict.items():
                SeqIO.write(SeqRecord(Seq(peptide), id=finalId, description=""), handle, "fasta")

def concatOutput(saveHandle, fileCount):
    # iterate through each file, excluding the most recently written file.
    for i in range(1, fileCount+1):
        # declare the path of the file for this iteration
        finalPath = str(saveHandle)[0:-17] + '_' + str(i) + '_' + str(saveHandle)[-17:-6]
        # remove subseqs from that file.
        removeSubsetSeq(finalPath)


def writeOutputFiles(finalSeenPeptides, protToPepFlag, pepToProtFlag, transFlag, outputPath, fileCount):
    """
    Called from writer() every time an output file is to be written. This file simply writes finalSeenPeptides to
    an output fasta file, and edits the name of this file to include what number the file is. If pepToProtFlag is True,
    it will also write a csv file with the input proteins as headings, and all peptides to be created from it listed
    below. If protToPepFlag is True, a csv file with peptides as headers and all proteins the peptide originated from
    listed below.

    :param finalSeenPeptides: a dictionary with peptide sequences as keys and a list of the proteins the peptide
    originated from as values.
    :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
    they originated in listed underneath.
    :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
    which they produced listed underneath.
    :param transFlag: True if the writer process is dealing with a trans splicing output. This allows the trans origin
    data, which is distinct from cis and linear, to be formatted sufficiently for output.
    :param outputPath: the generic output path that this data is to be written to. It must first be updated to include
    the file number.
    :param fileCount: the number of output files that will have been written once this one has been completed.
    :return:
    """
    finalPath = Path(str(outputPath)[0:-17] + '_' + str(fileCount) + '_' + str(outputPath)[-17:])
    with open(finalPath, 'w') as output_handle:
        # generate backwardSeenPeptides if protToPep is selected
        backwardsSeenPeptides = {}
        if protToPepFlag:
            # convert seen peptides to backwardsSeenPeptides
            for key, value in finalSeenPeptides.items():
                # check if we are printing trans entries so that we can configure trans data
                # before it goes into backwardsSeenPeptides
                if transFlag:
                    origins = editTransOrigins(value)
                else:
                    origins = value

                # Come back to make this less ugly and more efficient
                for entry in origins:
                    if entry not in backwardsSeenPeptides.keys():
                        backwardsSeenPeptides[entry] = [key]
                    else:
                        backwardsSeenPeptides[entry].append(key)
            writeProtToPep(backwardsSeenPeptides, 'ProtToPep', finalPath)

        if pepToProtFlag:
            writeProtToPep(finalSeenPeptides, 'PepToProt', finalPath)

        logging.info("Writing to fasta")
        SeqIO.write(createSeqObj(finalSeenPeptides), output_handle, "fasta")

def writeProtToPep(seenPeptides, groupedBy, outputPath):
    """
    Called by writeOutputFiles(), this function produces a csv file with either the proteins as headings followed by
    the peptides each protein produces underneath, or with the peptides as headings and the proteins each peptide
    could have come from listed underneath. Which format is used depends on the groupedBy flag, and the format of
    seenPeptides.

    :param seenPeptides: a dictionary containing peptide sequences and details of where each peptide was generated
    from. Which of these two variables is the key and which is the value depends on what the value of groupedBy is.
    :param groupedBy: will contain 'ProtToPep' if the output file is to have proteins as headings, with the
    corresponding produced peptides listed underneath. If so, seenPeptides.keys() will contain each origin protein
    and seenPeptides.values() will contain a list of peptides.
    This param will contain 'PepToProt' if the output file is to have the peptides as headings, with the corresponding
    origin peptides listed underneath. In this scenario, seenPeptides.keys() will contain peptides, and
    seenPeptides.values() will be a list of origin proteins.
    :param outputPath: the generic output location, which will be edited to include the groupedBy flag.
    :return:
    """
    newPath = outputPath.parent / (outputPath.name + groupedBy + '.csv')
    with open(newPath, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        if groupedBy is 'ProtToPep':
            header = 'Protein'
        else:
            header = 'Peptide'
        writer.writerow([header])
        for key, value in seenPeptides.items():
            infoRow = [key]
            writer.writerow(infoRow)
            for peptide in value:
                writer.writerow([peptide])
            writer.writerow([])

def editTransOrigins(origins):
    """
    Called by writeOutputFiles(), this function alters the trans origin data so that only proteins that contribute a
    cleavage of at least 6 amino acids are include in the ProtToPep csv.

    :param origins: a list of origins. It will have the syntax:
    ['Overlap/P04439(13-18)', 'Overlap/P04439(13-19)', 'P04439(13-19)/P30456']
    :return list(set(newOrigins)): a list unique the protein/reference strings which refer to cleavages that contributed
    at least 6 amino acids to the creation of a given trans peptide. The list will have the following format:
    [P04439(13-18), P04439(13-19)]
    """
    newOrigins = []
    for entry in origins:
        prots = entry.split('/')
        for prot in prots:
            if prot[-1] == ')':
                newOrigins.append(prot)
    return list(set(newOrigins))

def createSeqObj(matchedPeptides):
    """
    Given the dictionary of matchedPeptides, converts all of them into SeqRecord objects and adds them to an iterable
    via the yield statement. This iterable is passed back to the SeqIO.write() method which writes the created
    SeqRecord to fasta file.

    :param matchedPeptides: a dictionary with peptide sequences as keys and a list of origin proteins as values.
    :return:
    """
    # declare counter variable
    count = 1

    # iterate through each entry in matchedPeptides
    for sequence, value in matchedPeptides.items():
        # initialise the final peptide name that is to be written to fasta.
        finalId = "ipd|pep"+str(count)+';'
        # iterate through each origin protein and add it to finalId
        for protein in value:
            finalId+=protein+';'

        # return an iterable where each element is a record to written to fasta file.
        yield SeqRecord(Seq(sequence), id=finalId, description="")

        # index the count, ready for the next entry.
        count += 1

def applyMods(combineModlessDict, modList, maxMod):

    """
    Called from genMassDict() and transProcess() to add modified peptides massDict. This function iterates through
    each modification and each amino acid per modification, accesses the finalModTable and calls genericMod() for
    each amino acid to be modified. Note that it applies additional modifications to already modified peptides to
    ensure all potential modifications are produced. It then appends the output of generic mod to the input massDict
    and returns massDict back to the calling function.

    :param combineModlessDict: the dictionary created by genMassDict() and trancProcess() which stores the
    created splice peptides along with location data and mass and charge data. At this stage of computation, massDict
    has the following syntax: massDict[PEPTIDE] = [monoisotopic mass, [referenceLocation]].
    :param modList: a list of the modifications input by the user to be applied to all generated peptides.
    :param maxMod: the maximum number of modifications that can be made per peptide.
    :return combineModlessDict: the input dictionary with the additional modified peptides appended to the end.
    """

    modNo = 0
    for mod in modList:
        # Keep track of which modification is taking place
        modNo += 1

        # Don't need to worry about it if no modification!
        if mod != 'None':
            # Get the list of modifications taking place
            aminoList = finalModTable[mod]
           
            # Go through each character in the modification one by one
            for i in range(0, len(aminoList) - 1):

                char = aminoList[i]
                massChange = aminoList[-1]
                # get the dictionary of mods and their mass
                modDict = genericMod(combineModlessDict, char, massChange, str(modNo), maxMod)
                # Add it to the current list!
                combineModlessDict.update(modDict)
    return combineModlessDict

def genericMod(combineModlessDict, character, massChange, modNo, maxMod):
    """
    Called from applyMods(), this function creates a dictionary containing all the modified peptides which are formed
    by changing the input character. It also applies the relevant mass change to the monoisotopic mass of the
    unmodified peptide to calculate the modified peptide mass. The modified peptides include all permutations of the
    modification within the max mod per peptide set by the user:
    Eg: modifying M in MAMAM will create m1AMAM, MAm1AM, MAMAm1 etc for all permutation.

    :param combineModlessDict: the current version massDict containing all peptides, including peptides with
    modifications that have already been run.
    :param character: the amino acid to be modified.
    :param massChange: the change in mass which results this modification being applied to a single amino acid.
    :param modNo: used to differentiate between different types of modification. The first modification input by
    the user will be 1, the second 2 and the last 3. A modified amino acid is represented by the lower case letter
    followed by the modNo. Eg: MAMAm1 means that the last M has been modified by the first modification.
    :param maxMod: the maximum number of modified peptides allowed per peptide, as input by the user.

    :return modDict: a dictionary of the same structure as massDict which holds all the new modified peptides. These
    are returned to applyMods and added to combineModlessDict ready to either be returned to genMassDict() or
    transProcess(), or put back into this function for the next modification to be applied.
    """

    # set maxMod to 100 if 'None' has been input by the user.
    if maxMod == 'None':
        maxMod = 100
    else:
        maxMod = int(maxMod)

    # initialise the modDict
    modDict = {}

    # Go through each peptide and mod it if necessary
    for string in combineModlessDict.keys():
        # Only need to mod it if it exists (ie : A in ABC)
        if character in string:
            # extract the currentMass so that massChange can be added to it
            currentMass = combineModlessDict[string][0]
            # count how many of the current mod have already been applied to this peptide. This will occur if
            # more than one amino acid is subject to the mod, and is important to know when calculating the mass.
            modsInOrig = string.count(modNo)
            # extract the location data of the peptide which will identical for the modified peptide.
            currentRef = combineModlessDict[string][1]

            # count the number of times the amino acid to be modified occurs in the peptide.
            numOccur = string.count(character)
            # declare the seqToIter variable, which will be added to after separate iterations to ensure all
            # permutation of the modification are applied to this peptide.
            seqToIter = [string]

            # Algorithm for generating all permutations with the mods
            for j in range(numOccur, 0, -1):
                newMods = []
                for seq in seqToIter:
                    # print(seq)
                    noOfMods = seq.count(modNo) + 1
                    if noOfMods > maxMod:
                        continue
                    newMass = currentMass + massChange * (noOfMods - modsInOrig)
                    temp = nth_replace(seq, character, character.lower() + modNo, j)
                    newMods.append(temp)
                    # print(newMods)
                    newValue = [newMass, currentRef]
                    # if trans is running, the original protein tuple must be updated with the modified peptide
                    try:
                        newValue.append(combineModlessDict[string][2])
                    except IndexError:
                        pass
                    modDict[temp] = newValue
                seqToIter += newMods
    return modDict


def linCisPepCheck(refs, transOrigins):
    """
    Called by combineTransPeptide() and combineOverlapPeptide(). This function checks if a given peptide location
    corresponds to a linear (when cis and trans splicing is running) or cis (only trans is running) peptide. This is
    done so that when trans is being run linear and cis peptides are recognised in the combination process and added
    to a set to be deleted from the final output. The same applies for linear peptides in the cis output.

    :param refs: a list of integers corresponding to where a certain peptide was located within the input
    protein/proteins
    :param transOrigins: if trans is being run, this contains a list of index pairs denoting the start and end
    position of each individual protein within the concatenated sequence. Structure: [[0,150], [151,205] ... ].
    If cis is being run, it is set to False. Thus, this variable also serves to ensure only linear peptides are
    recognised when cis splicing is being run.

    :return bool: returns True if the peptide is found to be cis/linear spliced, returns False if it is not.
    """
    # if transOrigins is not False, we know we are checking from a trans process and we need to check it the pep is Cis or Linear.
    if transOrigins != False:
        # return the protein and index that the peptide is from.
        prot1, index1 = findInitProt(refs[0] - 1, transOrigins)
        prot2, index2 = findInitProt(refs[-1] - 1, transOrigins)
        # if the two proteins are the same move to the next check
        if prot1 == prot2:
            # if the set is the same length as the list there is no overlap, and thus it is eiter cis/lin splice.
            if len(set(refs)) == len(refs):
                return True
        return False
    # if transOrigins = False, we are checking from a cis process and need to check if the splice is linear.
    else:
        prevRef = refs[0]
        # iterate through each ref starting from the second one.
        for i in range(1,len(refs)):
            # if at any point a ref is not one more than the previous one, return False.
            if refs[i] == prevRef + 1:
                prevRef = refs[i]
            else:
                return False
        # if we iterate all the way through without returning False, the peptide is Linear.
        return True

def chargeIonMass(massDict, chargeFlags):

    """
    Called by genMassDict() and transProcess(). This function takes the massDict after modification have been applied
    and converts the monoisotopic mass of each peptide to m/z ratios for each charge states selected by the user.
    These m/z ratios are stored in a dictionary which is added to the third element of the values of massDict.

    :param massDict: the dictionary storing the peptide sequence, monoisotopic mass and location data. Has the form:
    massDict[PEPTIDE] = [monoisotopic mass, [referenceLocation]]
    :param chargeFlags: the charge states that are to be considered as input by the user. Stored as a list of booleans,
    where the first element corresponds to if +1 is to be considered, the second to +2 and so on.
    Eg: [True, False, True, False, True]

    :return chargeMassDict: the updated dictionary with charge states and m/z ratios included. Has the form:
    massDict[PEPTIDE] = [monoisotopic mass, [referenceLocation], {charge: m/z ratio}]
    """

    chargeMassDict = {}

    # iterate through each item in massDict
    for key, value in massDict.items():
        # initialise dictionary for storing charge and corresponding m/z for each peptide.
        chargeAssoc = {}
        # iterete through each bool in chargeFlags. Use an index as this will directly correspond to the charge state.
        for z in range(0, len(chargeFlags)):
            # if the bool is True we need to add the m/z ratio for the given charge.
            if chargeFlags[z]:
                # mass charge function returns the m/z ratio given monoisotopic mass and charge.
                chargeMass = massCharge(value[0], z+1)  # +1 to z for actual charge value

                # if no mgf data simply add the chargeMass to the chargeAssoc dictionary.
                if mgfData is None:
                    chargeAssoc[z+1] = chargeMass

                # Make sure the chargemass is less than the maximum possible charge mass in the mgf
                elif chargeMass <= mgfData.chargeMaxDict[z+1]:
                    chargeAssoc[z+1] = chargeMass

        # if ChargeAssoc has been added to, we update the value and add both the key and newVal to an updated dict
        if chargeAssoc:
            # Add it to the 2 as the rest of code acceses it at index 2
            newVal = value
            newVal.insert(2, chargeAssoc)
            chargeMassDict[key] = newVal

    return chargeMassDict

def massCharge(predictedMass, z):
    """
    Called by chargeIonMass() this function calculates the m/z ratio which corresponds to the input monoisotopic mass
    and charge state.

    :param predictedMass: the monoisotopic mass which is to be converted to an m/z ratio.
    :param z: an integer representing the charge state, where i corresponds to +i.

    :return chargeMass: the calculated m/z ratio.
    """
    chargeMass = (predictedMass + (z * 1.00794))/z
    return chargeMass


def maxDistCheck(ref1, ref2, maxDistance):
    """
    called by combineOverlapPeptide(). Given the locations of two cleavages within the input protein, this function
    calculates if the two references are within the maxDistance or not. If they are within the maxDistance, the
    function returns True and maxDistCheck will continue with combining the cleavages to produce a cis spliced peptide.

    :param ref1: a list of integers corresponding to the location of a given cleavage/split.
    :param ref2: a list of integers corresponding to the location of a second cleavage/split which is potentially
    to be combined with the first.
    :param maxDistance: an integer containing the maximum distance that two cleavages can be a part if they are to be
    combined to produce a cis spliced peptide. the distance between the two is considered to be the distance between
    the two closest peptides. It will be set 'None' if the maximum distance is infinte.

    :return bool: True if the peptides are within the maxDistance, False if they are too far apart.
    """
    if maxDistance == 'None':
        return True
    # max distance defined as the number of peptides between to peptide strands
    valid = ref2[0] - ref1[-1] - 1
    if valid > maxDistance:
        return False
    return True


def maxSize(split, maxed):

    """
    Called by combinedCheck to ensure that a given cleavage or peptide is smaller than a given maxSize.

    :param split: the cleavage or peptide that is to have its size checked against max size.
    :param maxed: the max size that the cleavage or peptide is allowed to be.

    :return bool: False if the split is longer than maxed, True if it is shorter.
    """

    if len(split) > maxed:
        return False
    return True


def minSize(split, mined):

    """
    Called by combineCheck() to ensure that a given cleavage or peptide is larger than a given minSize.

    :param split: the cleavage or peptide that is to have its size checked against max size.
    :param mined: the min size that the cleavage or peptide is allowed to be.

    :return bool: False if the split is shorter than mined, True if it is longer.
    """

    if len(split) < mined:
        return False
    return True


def combineCheck(split, mined, maxed, ref1, ref2, maxDistance='None'):
    """
    Called by splitDictPeptide(), splitTransPeptide(), combineOverlapPeptide() and combineTransPeptide() to check that
    a peptide meets the min length, max length and max distance criteria input by the user. It returns True if all
    of these checks pass, and False if at least one is False.

    :param split: the cleavage or peptide that is to have its size checked against max size.
    :param mined: the min size that the cleavage or peptide is allowed to be.
    :param maxed: the max size that the cleavage or peptide is allowed to be.
    :param ref1: a list of integers corresponding to the location of a given cleavage/split.
    :param ref2: a list of integers corresponding to the location of a second cleavage/split.
    :param maxDistance: an integer containing the maximum distance that two cleavages can be a part if they are to be
    combined to produce a cis spliced peptide. the distance between the two is considered to be the distance between
    the two closest peptides. It will be set 'None' if the maximum distance is infinte.

    :return bool: True if maxSize(), minSize() and maxDistCheck() functions all return True, False if any of them
    return False.
    """
    booleanCheck = maxSize(split, maxed) and minSize(split, mined) and maxDistCheck(ref1, ref2, maxDistance)
    return booleanCheck

def overlapComp(ref1, ref2):
    """
    Called by combineOverlapPeptide(), this function checks if two lists contain any identical elements. The input data
    in this case will be lists containing cleavage location indexes. Thus, this function checks if the cleavages
    corresponding to input location references share any amino acids in the origin protein.

    :param ref1: a list of integers corresponding to the location of a given cleavage/split.
    :param ref2: a list of integers corresponding to the location of a second cleavage/split which is to potentially
    be combined with the first.

    :return bool: True if there is no overlap between refs, False if there is overlap.
    """
    S1 = set(ref1)
    S2 = set(ref2)
    if len(S1.intersection(S2)) == 0:
        return True
    return False


def addSequenceList(input_file, multiFileFlag):
    """
    Called by transProcess(). This function takes a path to a fasta file and stores all input records in a dictionary.
    The dictionary has the protein name as a key and protein sequence as a value. The name is editted so that it
    includes the file name if more than one file has been uploaded
    :param input_file: a file path which contains a fasta file of proteins.
    :param multiFileFlag: set to True if the user as uploaded more than one file, and in turn information on which
    file a protein comes from has to be added to the name.

    :return sequenceDictionary: the dictionary containing the protein names as keys and protein sequences as values for
    all proteins in the input fasta file.
    """

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    sequenceDictionary = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        name = name.split('|')[1]
        name = name.split(';')[0]
        if multiFileFlag:
            fileName = input_file.split('/')[-1]
            name = fileName.split('.')[0] + '_' + name

        sequenceDictionary[name] = sequence
    return sequenceDictionary

def combMass(combine, combineRef, origProtTups = None):
    """
    Called by genMassDict() and transProcess() to convert a list of peptides and their location references to a
    dictionary. It also calculates the monoisotopic mass of each peptide and adds this to the dictionary. It firstly
    checks that the mass is not greater than the max monoisotopic mass in the input MGF file, and if it is, this
    peptide is discarded. The output dictionary has the format:
    massDict[PEPTIDE] = [monoisotopic mass, [referenceLocation]]

    :param combine: a list of peptides.
    :param combineRef:  list of the locations of corresponding to where the peptides in combine were situated within
    the input protein. Has the form:[[1,2,3,4], [2,3,4,5]... ]
    :param origProtTups: if trans is being run, the origin data is stored in a series of tuples instead of the
    references contained in combineRef. Each peptide has a pair of tuples, with each tuple storing the origin protein
    and location within that protein of a cleavage used to create the peptide.
    Eg: [[(prot1, "(1-6)"), (prot2,"")], [(prot3, "(1-6)"), (prot4, "(10-16)")]....]

    :return massDict: a dictionary with peptide sequences as keys, and a list containing the monoisotopic mass and
    location data as values. Eg: massDict[PEPTIDE] = [monoisotopic mass, [referenceLocation], originTups]. Note
    that the value will not have a third element "originTups" when cis and linear are being run.
    """

    massDict = {}

    # if mgfData exists, set maxMass to equal it.
    try:
        maxMass = mgfData.maxMass
        # print('in maxmass')
    # if it doesn't exist, let maxMass == 1000000 (outrageously large)
    except:
        maxMass = 1000000

    # iterate through each peptide/location
    for i in range(0, len(combine)):
        # initialise totalMass
        totalMass = 0
        # iterate through each amino-acid, extract its mass from the monoAminoMass dictionary and add it to totalMass.
        for j in range(0, len(combine[i])):
            totalMass += monoAminoMass[combine[i][j]]
        # add the mass of water to get final monoisotopic mass.
        totalMass += H20_MASS
        # check if totalMass is greater than the max mass in the MGF file. If so continue as it will not match to
        # any spectra and thus won't appear in the output.
        if totalMass > maxMass:
            print(combine[i])
            continue
        # if originProtTups == None, we are dealing with a cis/linear splice and do not need to consider originProtTups.
        if origProtTups == None:
            massRefPair = [totalMass, combineRef[i]]
            massDict[combine[i]] = massRefPair
        # when trans is being run,  we need to configure a third element in the value containing the origin protein
        # information in the form of a pair of tuples.
        else:
            # check if this specific peptide has already been added to the massDict (not that this is only possible
            # trans splicing) and add the value to the dictionary in appropriate fashion.
            if combine[i] in massDict.keys():
                    massDict[combine[i]][2].append(origProtTups[i][0])
                    massDict[combine[i]][2].append(origProtTups[i][1])

            else:
                massRefPair = [totalMass, combineRef[i], origProtTups[i]]
                massDict[combine[i]] = massRefPair
    return massDict


def changeRefToDash(ref):
    """
    Called by editRefMassDict, this function takes a location reference in the form [1,2,3,4,5,6] and converts it to a
    range [1-6]

    :param ref: the location reference to be converted to a range: [1,2,3,4,5,6].\
    :return newRef: the location reference in the form of a range: [1-6]
    """
    newRef = []
    for i in range(0,len(ref)):
        refVal = ref[i]
        # check if last element reached and thus input is linear
        if i == len(ref)-1:
            newLinRef = str(ref[0]) + ' - ' + str(ref[-1])
            newRef = [newLinRef]
            return newRef
        # otherwise, check if the next element is still sequential, and if so continue for loop
        if refVal + 1 == ref[i+1]:
            continue
        else:
            if i == 0:
                newStrRef1 = str(ref[0])
            else:
                newStrRef1 = str(ref[0]) + "-" + str(ref[i])
            if i + 1 == len(ref) - 1:
                newStrRef2 = str(ref[-1])
            else:
                newStrRef2 = str(ref[i+1]) + "-" + str(ref[-1])

            newRef = [newStrRef1, newStrRef2]
            return newRef


def editRefMassDict(massDict):
    """
    Called by genMassDict() and transProcess(), this function converts all the location lists contained in the
    value of massDict to a range. That is, from the form [1,2,3,4,5] to [1-5]

    :param massDict: the dictionary containing peptide sequences as keys and a list of location data, monoisotpic mass
    and m/z ratios as the value. Has the form:
    massDict[PEPTIDE] = [monoisotopic mass, [referenceLocation], {charge: m/z ratio}]

    :return:
    """
    for key, value in massDict.items():
        refNew = changeRefToDash(value[1])
        value[1] = refNew
    return massDict

def nth_replace(string, old, new, n=1, option='only nth'):

    # https://stackoverflow.com/questions/35091557/replace-nth-occurrence-of-substring-in-string
    """
    This function replaces occurrences of string 'old' with string 'new'.
    There are three types of replacement of string 'old':
    1) 'only nth' replaces only nth occurrence (default).
    2) 'all left' replaces nth occurrence and all occurrences to the left.
    3) 'all right' replaces nth occurrence and all occurrences to the right.

    :param string: the string from which you want to replace a substring.
    :param old: the substring you want replaced.
    :param new: the substring you want to replace the old substring with.
    :param n: which occurence of the old string you want replaced, will work in conjunction with the only nth flag.
    :param option: if you want to replace only the nth occurence, all to the right of the nth ocurrence or all to the
    left of the nth occurence.
    :return:
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


def writeToCsv(massDict, header, finalPath, chargeFlags):
    """
    Called from transProcess() and genMassDict(), this function writes the data stored in massDict to a csv file
    after all spliced peptides have been generated. It writes the data from each protein sequence separately, with the
    column headers being sequence, mass, position, m/z +1, m/z +2 ... m/z +5.

    :param massDict: the mass dictionary that is to have its data written to file.
    :param header: the name of the protein which this massDict corresponds to.
    :param finalPath: the path the csv file is to be written to.
    :param chargeFlags: which charge states have been selected by the user for consideration.
    :return:
    """
    chargeHeaders = getChargeIndex(chargeFlags)

    with open(finalPath, 'a', newline='') as csv_file:

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
    """
    Called from writeToCsv(), this function returns a list of the charges that have been selected by the user in the
    form [1,3,5] etc.

    :param chargeFlags: a list of bools, where bool at index i denotes if the user has selected to consider charge
    state i+1. That is, [True, False, False, False, True] means the user has chosen to consider +1 and +5 charge.
    :return chargeHeaers: a list of which charges are being considered. Eg: [1, 5]
    """
    chargeHeaders = [i for i, e in enumerate(chargeFlags) if e]
    return chargeHeaders


def getFinalPath(outputPath):
    """
    Called from transProcess() and genMassDict() if the user has selected to create the writeToCsv output. This function
    simply replaces the .fasta extention of outputPath with .csv.

    :param outputPath: the original output path with the .fasta extention.
    :return: the new path with .csv extention as a Path() object.
    """
    outputPathSmall = str(outputPath)[0:-6]
    newPath = str(outputPathSmall) + '.csv'
    return Path(newPath)

def processLockInit(lockVar, toWriteQueue, mgfObj, childTable, linSetQueue):

    """
    Called by self.cisAndLinearOutput() before the multiprocessing pool is created to set up a global lock for a child
    processes and to give all processes access to important queues and shared variables.

    :param lockVar: the multiprocessing.Lock() variable used for to control the access of process to certain tasks.
    :param toWriteQueue: the queue which all processes need to be able to access to enable them to push their output
    to the writer() function.
    :param mgfObj: the object which stores all the data needed from the mgf file, which all processes require access
    to. It was too large to initialise for each process, thus we can initialise it once by giving all processes
    global access.
    :param childTable: a list of all modifications, which has been updated to include custom modifications.
    :param linSetQueue: the queue which all processes require access to so they can push all generated lin peptides
    to the writer() function for deletion from the final cis output.
    :return:
    """

    global lock
    lock = lockVar
    global mgfData
    mgfData = mgfObj
    global finalModTable
    finalModTable = childTable
    genMassDict.toWriteQueue = toWriteQueue
    genMassDict.linSetQueue = linSetQueue

def processLockTrans(lockVar, toWriteQueue, allSplits, allSplitRef, mgfObj, childTable, linCisQueue):

    """
    Called by self.transOutput() before the multiprocessing pool is created to set up a global lock for a child
    processes and to give all processes access to important queues and shared variables.

    :param lockVar: the multiprocessing.Lock() variable used for to control the access of process to certain tasks.
    :param toWriteQueue: the queue which all processes need to be able to access to enable them to push their output
    to the writer() function.
    :param allSplits: a list of all the generate linear cleavages from all input peptides which will be used to
    generate all trans spliced peptide.
    :param allSplitRef: a list of references to the locations of allSplits peptides within the sequence resulting
    from concetenating all input proteins.
    :param mgfObj: the object which stores all the data needed from the mgf file, which all processes require access
    to. It was too large to initialise for each process, thus we can initialise it once by giving all processes
    global access.
    :param childTable: a list of all modifications, which has been updated to include custom modifications.
    :param linSetQueue: the queue which all processes require access to so they can push all generated lin/cis peptides
    to the writer() function for deletion from the final output.
    :return:
    """
    global lock
    lock = lockVar
    global mgfData
    mgfData = mgfObj
    global finalModTable
    finalModTable = childTable
    global splits
    splits = allSplits
    global splitRef
    splitRef = allSplitRef
    transProcess.toWriteQueue = toWriteQueue
    transProcess.linCisQueue = linCisQueue

