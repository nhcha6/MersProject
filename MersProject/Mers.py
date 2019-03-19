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

# define our the spliceType flags
TRANS = "Trans"
LINEAR = "Linear"
CIS = "Cis"

# MEMORY_THRESHOLD is the percentage of RAM being used at which all computation is paused, the current output data is
# written to file, after which the output recommences.
MEMORY_THRESHOLD = 85
# ** obsolete now that temp files are not being used.
MEMORY_THRESHOLD_COMBINE = 90
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

        self.inputFile = [inputFile]
        self.procGenCounter = 0
        self.pepCompleted = multiprocessing.Queue()
        self.completedProcs = 0
        self.totalProcs = 0

    def generateOutput(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                       protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfObj, mgfFlag):

        """
        This function is called from MersGUI.py to initiate the output process once the user input has been recieved
        and confirmed suitable.

        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide.
        :param transFlag: if True, trans splicing should be completed as part of the output.
        :param cisFlag: if True, cis splicing should be completed as part of the output.
        :param linearFlag: if True, linear splicing should be completed as part of the output.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted. ** We may look to consider removing this functionality.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath
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
            self.transOutput(self.inputFile, TRANS, mined, maxed, modList, maxMod,outputPath[TRANS], chargeFlags, mgfObj,
                        modTable, mgfFlag, csvFlag, pepToProtFlag,protToPepFlag)

        if cisFlag:
            self.cisAndLinearOutput(self.inputFile, CIS, mined, maxed,overlapFlag, csvFlag, pepToProtFlag,
                                    protToPepFlag, modList, maxMod, maxDistance, outputPath[CIS], chargeFlags, mgfObj,
                                    modTable, mgfFlag)

        if linearFlag:
            self.cisAndLinearOutput(self.inputFile, LINEAR, mined, maxed, overlapFlag, csvFlag, pepToProtFlag,
                                    protToPepFlag, modList, maxMod, maxDistance, outputPath[LINEAR], chargeFlags,
                                    mgfObj, modTable, mgfFlag)


    def transOutput(self, inputFile, spliceType, mined, maxed, modList, maxMod, outputPath, chargeFlags, mgfObj, modTable, mgfFlag, csvFlag,
                    pepToProtFlag, protToPepFlag):
        """
        This function is called Fasta.generateOutput() if transFlag == True. It controls the creation of all processes
        required to conduct trans splicing.

        :param inputFile: the fasta files containing proteins input by the user.
        :param spliceType: simply holds the TRANS flag. **not needed, as this function is only called by trans splicing.
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
        :param modTable: the modTable dictionary contianing all the modification data. ** does this need to be an input
        variable, it is globally defined??
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted. ** We may look to consider removing this functionality.
        :param pepToProtFlag: if True, a csv file is written with the output peptides as a heading, and the proteins
        they originated in listed underneath.
        :param protToPepFlag: if True, a csv file is written with the input proteins as a heading, and the peptides
        which they produced listed underneath
        :return:
        """

        # initialise the final path
        finalPath = None

        # Open the csv file if the csv file is selected **do we still need csv
        if csvFlag:
            finalPath = getFinalPath(outputPath, spliceType)
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
        finalPeptide, protIndexList, protList = combinePeptides(seqDict)

        # create all the splits and their reference with respect to finalPeptide
        splits, splitRef = splitTransPeptide(spliceType, finalPeptide, mined, maxed, protIndexList)

        # define the length of splits for later calculations
        splitLen = len(splits)

        # define the number of workers to be input into the multiprocessing.Pool() based on the number of CPUs in the
        # computer.
        num_workers = multiprocessing.cpu_count()

        # Used to lock write access to file when writing from separate processes ** do we still need this.
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
                                                                     protToPepFlag, self.pepCompleted, True))
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
        logging.info("All " + spliceType + " !joined")

    def cisAndLinearOutput(self, inputFile, spliceType, mined, maxed, overlapFlag, csvFlag, pepToProtFlag, protToPepFlag,
                           modList, maxMod, maxDistance, outputPath, chargeFlags, mgfObj, childTable, mgfFlag):

        """
        This function is called Fasta.generateOutput() if transFlag == True. It controls the creation of all processes
        required to conduct trans splicing.

        :param inputFile: the fasta files containing proteins input by the user.
        :param spliceType: simply holds the TRANS flag. **not needed, as this function is only called by trans splicing.
        :param mined: the minimum length of an output peptide.
        :param maxed: the maximum length of an ouptut peptide.
        :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
        permitted. Shared amino acids originate from the same amino-acid in the same peptide.
        :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
        or precursor mass comparison is conducted. ** We may look to consider removing this functionality.
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
        :param childTable: the modTable dictionary contianing all the modification data. ** does this need to be an input
        variable, it is globally defined??
        :param mgfFlag: if True, no the user has selected that no mgf comparison be conducted and the raw splice data
        is to be ouptut to Fasta.

        :return:
        """

        # declare final path.
        finalPath = None

        # Open the csv file if the csv file is selected
        if csvFlag:
            finalPath = getFinalPath(outputPath, spliceType)
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
                                    initargs=(lockVar, toWriteQueue, mgfObj, childTable, linSetQueue))
        writerProcess = multiprocessing.Process(target=writer,
                                                args=(toWriteQueue, outputPath, linSetQueue, pepToProtFlag,
                                                      protToPepFlag, self.pepCompleted))
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

                    # add the filename to the seqId if more than one file has been included
                    seqId = seqId.split('|')[1]
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
                                                        initargs=(lockVar, toWriteQueue, mgfObj, childTable,
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
    Takes a dictionary with protein names as keys and proteins sequences as values and concatenates the sequences
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

def splitTransPeptide(spliceType, peptide, mined, maxed, protIndexList):

    """
    Creates all the possible splits (peptide cleavages) from the input protein sequence. This differs from
    splitDictPeptide due to functionality required specifically by trans. The input protein to this function
    is all the inividual input proteins concatenated together, creating potential splits across the border between
    two proteins which cannot be formed given the input data. The protIndexList data is thus included to ensure that
    such cleavages are ignored from the fonal splits output.

    :param spliceType: will always be the TRANS flag ** need to remove.
    :param peptide: the input proteins sequence, which will be all individual proteins concatenated together.
    ** rename to proteinSeq.
    :param mined: the min length of a potential output peptide.
    :param maxed: the max length of a potential output peptide.
    :param protIndexList: a list of location pairs denoting the start and end point of original, individual proteins.
    :return splits: all the potential cleavages which could be recombined to create a trans splice peptide.
    :return splitRef: where each split originated in the input protein.
    """

    # Makes it easier to integrate with earlier iteration where linearFlag was being passed as an external flag
    # instead of spliceType
    linearFlag = spliceType == LINEAR
    length = len(peptide)

    # splits will hold all possible splits that can occur
    splits = []
    # splitRef will hold a direct reference to the characters making up each split string: for starting peptide ABC,
    # the split AC = [0,2]
    splitRef = []

    # embedded for loops build all possible splits
    for i in range(0, length):

        character = peptide[i]
        toAdd = ""

        # figure out which protein the splits starts in, and the max index the splits can reach before it becomes
        # a part of a second peptide
        initProt, protInd = findInitProt(i, protIndexList)

        # add and append first character and add and append reference number which indexes this character
        toAdd += character
        ref = list([i+1])
        temp = list(ref)  # use list because otherwise shared memory overwrites

        # linear flag to ensure min is correct for cis and trans
        if linearFlag and minSize(toAdd, mined):

            # Don't need to continue this run as first amino acid is unknown X
            if 'X' in toAdd or 'U' in toAdd:
                continue
            else:
                splits.append(toAdd)
                splitRef.append(temp)

        elif not linearFlag and 'X' not in toAdd and 'U' not in toAdd:
            splits.append(toAdd)
            splitRef.append(temp)

        # iterates through every character after current and adds it to the most recent string if max size
        # requirement is satisfied
        for j in range(i + 1, length):
            if j > initProt[1]:
                break
            toAdd += peptide[j]
            if linearFlag:
                ref.append(j+1)
                if maxSize(toAdd, maxed):
                    if minSize(toAdd, mined):

                        # All future splits will contain X if an X is found in the current run, hence break
                        if 'X' in toAdd or 'U' in toAdd:
                            break
                        splits.append(toAdd)
                        temp = list(ref)
                        splitRef.append(temp)
                else:
                    break

            else:
                if maxSize(toAdd, maxed-1):
                    # All future splits will contain X if an X is found in the current run, hence break
                    if 'X' in toAdd or 'U' in toAdd:
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
    Controls the flow of computing the trans output of a given subset of splits. It calls functions to create
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
    or precursor mass comparison is conducted. ** We may look to consider removing this functionality.
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
            allPeptides = getAllPep(massDict)
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
            #fulfillPpmReq(mgfObj, massDict)
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
    Takes splits index from the multiprocessing pool and computes all the possible peptides which can be created
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
    Iterates through each combineRef and returns the origin protein names and the index within each protein for
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
    Takes an index to an amino acid within the concatenated input protein, and locates which range within protIndexList
    this reference belongs to. It then return this range of values (relating to the range the origin protein occupies
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
    Called by Fasta.cisAndLinearOutput(), this is the worker function which controls the computation and output of
    a given process for cis and linear splicing. It calls functions to create the cis or linear peptides, calculates
    relevant masses and m/z ratios, applies modifications and compares this data to an MGF file is desire.

    :param spliceType: holds either the CIS or TRANS flag so that the appropriate splicing can be computed.
    :param protDict:
    :param mined: the minimum length of an output peptide.
    :param maxed: the maximum length of an ouptut peptide.
    :param overlapFlag: if True, this flag denotes that combination of splits containing shared amino acids is not
    permitted. Shared amino acids originate from the same amino-acid in the same peptide.
    :param csvFlag: if True, the data contained in massDict should be printed to a csv file before any b-y Ion
    or precursor mass comparison is conducted. ** We may look to consider removing this functionality.
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
            # ** decide on relevance of this.
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
def outputCreate(spliceType, peptide, mined, maxed, overlapFlag, maxDistance=10000000):

    # Splits eg: ['A', 'AB', 'AD', 'B', 'BD']
    # SplitRef eg: [[0], [0,1], [0,2], [1], [1,2]]
    # Produces splits and splitRef arrays which are passed through combined
    splits, splitRef = splitDictPeptide(spliceType, peptide, mined, maxed)
    combined, combinedRef = None, None

    if spliceType == CIS:
        # combined eg: ['ABC', 'BCA', 'ACD', 'DCA']
        # combinedRef eg: [[0,1,2], [1,0,2], [0,2,3], [3,2,0]]
        # pass splits through combined overlap peptide and then delete all duplicates
        combined, combinedRef, linSet = combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance)

    elif spliceType == LINEAR:
        # Explicit change for high visibility regarding what's happening
        combined, combinedRef = splits, splitRef
        linSet = set()

    return combined, combinedRef, linSet

def splitDictPeptide(spliceType, peptide, mined, maxed):

    """
    Inputs: peptide string, max length of split peptide.
    Outputs: all possible splits that could be formed that are smaller in length than the maxed input
    """
    # Makes it easier to integrate with earlier iteration where linearFlag was being passed as an external flag
    # instead of spliceType
    linearFlag = spliceType == LINEAR
    length = len(peptide)

    # splits will hold all possible splits that can occur
    splits = []
    # splitRef will hold a direct reference to the characters making up each split string: for starting peptide ABC,
    # the split AC = [0,2]
    splitRef = []

    # embedded for loops build all possible splits
    for i in range(0, length):

        character = peptide[i]
        toAdd = ""
        # add and append first character and add and append reference number which indexes this character

        toAdd += character
        ref = list([i+1])
        temp = list(ref)  # use list because otherwise shared memory overwrites

        # linear flag to ensure min is correct for cis and trans
        if linearFlag and minSize(toAdd, mined):

            # Don't need to continue this run as first amino acid is unknown X
            if 'X' in toAdd or 'U' in toAdd:
                continue
            else:
                splits.append(toAdd)
                splitRef.append(temp)

        elif not linearFlag and 'X' not in toAdd and 'U' not in toAdd:
            splits.append(toAdd)
            splitRef.append(temp)

        # iterates through every character after current and adds it to the most recent string if max size
        # requirement is satisfied
        for j in range(i + 1, length):
            toAdd += peptide[j]
            if linearFlag:
                ref.append(j+1)
                if maxSize(toAdd, maxed):
                    if minSize(toAdd, mined):

                        # All future splits will contain X if an X is found in the current run, hence break
                        if 'X' in toAdd or 'U' in toAdd:
                            break
                        splits.append(toAdd)
                        temp = list(ref)
                        splitRef.append(temp)
                else:
                    break

            else:
                if maxSize(toAdd, maxed-1):
                    # All future splits will contain X if an X is found in the current run, hence break
                    if 'X' in toAdd or 'U' in toAdd:
                        break
                    splits.append(toAdd)
                    ref.append(j+1)
                    temp = list(ref)
                    splitRef.append(temp)
                else:
                    break

    return splits, splitRef


def combineOverlapPeptide(splits, splitRef, mined, maxed, overlapFlag, maxDistance):

    """
    Input: splits: list of splits, splitRef: list of the character indexes for splits, mined/maxed: min and max
    size requirements, overlapFlag: boolean value true if overlapping combinations are undesired.
    Output: all combinations of possible splits which meets criteria
    """
    # initialise linSet
    linSet = set()
    # initialise combinations array to hold the possible combinations from the input splits
    massDict = {}
    combModless = []
    combModlessRef = []
    # iterate through all of the splits and build up combinations which meet min/max/overlap criteria
    for i in range(0, len(splits)):

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
            if combineCheck(toAddForward, mined, maxed, splitRef[i], splitRef[j], maxDistance):
                # V. messy, need a way to get better visual
                if overlapFlag:
                    if overlapComp(splitRef[i], splitRef[j]):
                        massDict[toAddReverse] = addReverseRef
                        #check if toAdd forward is linear and add to linearSet if so
                        if linCisPepCheck(addForwardRef, False):
                            linSet.add(toAddForward)
                        else:
                            massDict[toAddForward] = addForwardRef


                else:
                    massDict[toAddReverse] = addReverseRef
                    # check if toAddForward is linear and add to linearSet if so
                    if linCisPepCheck(addForwardRef, False):
                        linSet.add(toAddForward)
                    else:
                        massDict[toAddForward] = addForwardRef

            elif not maxDistCheck(splitRef[i], splitRef[j], maxDistance):
                break

            toAddForward = ""
            toAddReverse = ""

    for peptide, ref in massDict.items():
        if peptide in linSet:
            continue
        else:
            combModless.append(peptide)
            combModlessRef.append(ref)

    return combModless, combModlessRef, linSet

def memoryCheck(maxMem):
    process = psutil.Process(os.getpid())
    #print(process.memory_info().rss)
    if process.memory_info().rss > maxMem:
        return True
    else:
        return False

def memoryCheck2():
    memUsed = psutil.virtual_memory()[2]
    #print(memUsed)
    if memUsed > 60:
        return True
    else:
        return False

def getAllPep(massDict):

    allPeptides = set()
    for key, value in massDict.items():
        if not key.isalpha():
            alphaKey = modToPeptide(key)
        else:
            alphaKey = key
        allPeptides.add(alphaKey)
    return allPeptides

def memory_usage_psutil():
    # return the memory usage in percentage like top
    mem = psutil.virtual_memory()

    return mem.percent

def writer(queue, outputPath, linCisQueue, pepToProtFlag, protToPepFlag, procCompleted, transFlag = False):

    seenPeptides = {}
    linCisSet = set()
    saveHandle = str(outputPath)
    modCountDict = Counter()
    fileCount = 0

    while True:
        # get from cisLinQueue and from matchedPeptide Queue
        matchedTuple = queue.get()

        if not linCisQueue.empty():
            linCisSet = linCisSet | linCisQueue.get()

        # if stop is sent to the matchedPeptide Queue, everything has been output,
        # so we exit the while loop.
        if matchedTuple == STOPFLAG:
            logging.info("Everything computed, stop message has been sent")
            break

        # If mem is sent, we know the max memory has been hit during process generation.
        if matchedTuple == MEMFLAG:
            # remove linear/cis peptides from seenPeptides:
            print('memflag recieved, writing temp file')
            commonPeptides = linCisSet.intersection(seenPeptides.keys())
            for peptide in commonPeptides:
                del seenPeptides[peptide]
            print('deleted common peptides')
            # write to ouptut
            fileCount += 1
            writeOutputFiles(seenPeptides, protToPepFlag, pepToProtFlag, transFlag, outputPath, fileCount)
            print('written to output')
            seenPeptides = {}
            # put back to procCompleted queue to tell process generating function that it can restart the queue.
            procCompleted.put(1)
            continue

        # flag gets put via writer queue everytime a process is finished
        if matchedTuple == PROC_FINISHED:
            procCompleted.put(1)
            continue

        # if trans, each queue.get that reaches here (ie isn't a flag of some-sort) corresponds to a process
        # thus, we need to index
        if transFlag:
            procCompleted.put(1)

        # if metchedTuple is not MEMFLAG or STOPFLAG, it is a genuine output and we continue as
        # normal to add it to seenPeptides.
        matchedPeptides = matchedTuple[0]
        if matchedTuple[1]:
            modCountDict += matchedTuple[1]

        # each queue.get() returns the matchedPeptides from an individual process.
        # Add  the matchedPeptides from the given process to seenPeptides.
        for key, value in matchedPeptides.items():
            origins = value.split(';')
            if key not in seenPeptides.keys():
                seenPeptides[key] = origins
            else:
                for origin in origins:
                    if origin not in seenPeptides[key]:
                        seenPeptides[key].append(origin)

    if seenPeptides:
        # Was not over memory threshold but last items need to be dealt with.
        commonPeptides = linCisSet.intersection(seenPeptides.keys())
        for peptide in commonPeptides:
            del seenPeptides[peptide]

        # write to ouptut
        fileCount += 1
        writeOutputFiles(seenPeptides, protToPepFlag, pepToProtFlag, transFlag, outputPath, fileCount)

    if modCountDict:
        # need to know if related to cis/lin/trans. Replace the relevant
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
        #print(infoPath)
        file = open(infoPath, 'a')
        file.write('\n' + title)
        for key, value in modCountDict.items():
            file.write(key + ': ' + str(value) + '\n')
        #print(modCountDict)

def writeOutputFiles(finalSeenPeptides, protToPepFlag, pepToProtFlag, transFlag, outputPath, fileCount):
    finalPath = str(outputPath)[0:-17] + '_' + str(fileCount) + '_' + str(outputPath)[-17:]
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
            writeProtToPep(backwardsSeenPeptides, 'ProtToPep', outputPath)

        if pepToProtFlag:
            writeProtToPep(finalSeenPeptides, 'PepToProt', outputPath)

        logging.info("Writing to fasta")
        SeqIO.write(createSeqObj(finalSeenPeptides), output_handle, "fasta")

def combineAllTempFasta(linCisSet, outputTempFiles, outputPath):
    counter = 0
    while not outputTempFiles.empty():

        # Get the two files at the top of the tempFiles queue for combination.
        # Note that there will never be one temp file in the queue when the
        # while loop is being checked, so you will always be able to get two
        # temp files from the queue if it passes the not empty check.
        fileOne = outputTempFiles.get()
        fileTwo = outputTempFiles.get()

        # if this reduces the queue to empty, break the loop. We do this to avoid merging
        # the last two temp files, adding the result to the queue and then passing a queue with
        # only one temp file in it into the while loop.
        if outputTempFiles.empty():
            break

        # when there are still more temp files in the queue, extract seenPeptides from the
        # current two temp files, write them to a new tempFile and add it to the temp file Queue.
        seenPeptides, counter = combineTempFile(fileOne, fileTwo, linCisSet, counter, outputPath)
        tempName = writeTempFasta(seenPeptides)
        outputTempFiles.put(tempName)

    # once the while loop breaks, return the finalSeenPetides from the remaining two tempFiles.
    finalSeenPeptides, counter = combineTempFile(fileOne, fileTwo, linCisSet, counter, outputPath)

    # Return the last combination of two files remaining
    return finalSeenPeptides


def combineTempFile(fileOne, fileTwo, linCisSet, counter, outputPath):
    seenPeptides = {}
    with open(fileOne, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):

            peptide = str(record.seq)
            protein = str(record.name)

            origins = protein.split(';')
            if peptide not in seenPeptides.keys():
                seenPeptides[peptide] = origins
            else:
                for origin in origins:
                    if origin not in seenPeptides[peptide]:
                        seenPeptides[peptide].append(origin)

    with open(fileTwo, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):

            peptide = str(record.seq)
            protein = str(record.name)

            origins = protein.split(';')
            if peptide not in seenPeptides.keys():
                seenPeptides[peptide] = origins
            else:
                for origin in origins:
                    if origin not in seenPeptides[peptide]:
                        seenPeptides[peptide].append(origin)

            if memory_usage_psutil() > MEMORY_THRESHOLD_COMBINE:
                if seenPeptides:
                    commonPeptides = linCisSet.intersection(seenPeptides.keys())
                    for peptide in commonPeptides:
                        del seenPeptides[peptide]

                outputPath = outputPath.split('.fasta')[0]
                finalPath = outputPath + "-file" + str(counter) + '.fasta'
                with open(finalPath, 'w') as output_handle:
                    SeqIO.write(createSeqObj(seenPeptides), output_handle, "fasta")
                seenPeptides = {}
                counter += 1

    # Delete temp files as they are used up
    os.remove(fileOne)
    os.remove(fileTwo)
    # In case final record
    if seenPeptides:
        commonPeptides = linCisSet.intersection(seenPeptides.keys())
        for peptide in commonPeptides:
            del seenPeptides[peptide]

    return seenPeptides, counter


def writeTempFasta(seenPeptides):
    temp = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fasta", delete=False)
    for key, value in seenPeptides.items():
        temp.writelines(">")
        counter = 0
        for protein in value:
            if counter != 0:
                temp.writelines(';')
            temp.writelines(str(protein))
            counter += 1
        temp.writelines("\n")
        temp.writelines(str(key))
        temp.writelines("\n")
    return temp.name

def writeProtToPep(seenPeptides, groupedBy, outputPath):
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
    newOrigins = []
    for entry in origins:
        prots = entry.split('/')
        for prot in prots:
            if prot[-1] == ')':
                newOrigins.append(prot)
    return list(set(newOrigins))


def fulfillPpmReq(mgfObj, massDict):
    """
    Assumption there are charges. Get the peptides that match, and writes them to the output fasta file
    """

    matchedPeptides = generateMGFList(mgfObj, massDict)

    lock.acquire()
    logging.info("Writing to fasta")
    with open("OutputMaster.fasta", "a") as output_handle:
        SeqIO.write(createSeqObj(matchedPeptides), output_handle, "fasta")

    lock.release()
    logging.info("Writing complete")


def createSeqObj(matchedPeptides):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []
    try:
        for sequence, value in matchedPeptides.items():

            finalId = "ipd|pep"+str(count)+';'

            for protein in value:
                finalId+=protein+';'

            yield SeqRecord(Seq(sequence), id=finalId, description="")

            count += 1
    except AttributeError:
        print(str(matchedPeptides))


    # return seqRecords

def applyMods(combineModlessDict, modList, maxMod):

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
    From the modless dictionary of possible combinations, this function returns a
    dictionary containing all the modifications that arise from changing character. The
    key of the output is simply the modified peptide, and the value is the mass which
    results as set by massChange
    """

    if maxMod == 'None':
        maxMod = 100
    else:
        maxMod = int(maxMod)

    # A, B, C  convert to ai, bi, ci where i is the modNo
    modDict = {}

    # Go through each combination and mod it if necessary
    for string in combineModlessDict.keys():
        currentMass = combineModlessDict[string][0]
        modsInOrig = string.count(modNo)
        currentRef = combineModlessDict[string][1]

        # Only need to mod it if it exists (ie : A in ABC)
        if character in string:

            numOccur = string.count(character)
            seqToIter = [string]
            # Generate all permutations with the mods
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
    chargeFlags: [True, False, True, False, True]
    """

    chargeMassDict = {}

    for key, value in massDict.items():
        chargeAssoc = {}
        for z in range(0, len(chargeFlags)):
            if chargeFlags[z]:
                chargeMass = massCharge(value[0], z+1)  # +1 for actual value
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
    chargeMass = (predictedMass + (z * 1.00794))/z
    return chargeMass


def writeToCsv(massDict, header, finalPath, chargeFlags):

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
    chargeHeaders = [i for i, e in enumerate(chargeFlags) if e]
    return chargeHeaders


def maxDistCheck(ref1, ref2, maxDistance):
    if maxDistance == 'None':
        return True
    # max distance defined as the number of peptides between to peptide strands
    valid = ref2[0] - ref1[-1] - 1
    if valid > maxDistance:
        return False
    return True


def maxSize(split, maxed):

    """
    ensures length of split is smaller than or equal to max
    """

    if len(split) > maxed:
        return False
    return True


def minSize(split, mined):

    """
    ensures length of split is greater than min
    """

    if len(split) < mined:
        return False
    return True


def combineCheck(split, mined, maxed, ref1, ref2, maxDistance='None'):
    booleanCheck = maxSize(split, maxed) and minSize(split, mined) and maxDistCheck(ref1, ref2, maxDistance)
    return booleanCheck


def linearCheck(toAdd, combinedLinearSet):
    # Look to do this better, not comparing to combineLinearSet, instead checking that the splitsRefs aren't linearly
    # ordered: [1, 2, 3] and [4, 5, 6] are obviously linearly ordered.
    if combinedLinearSet is None:
        return True
    if toAdd in combinedLinearSet:
        return False
    return True


def overlapComp(ref1, ref2):

    """
    checks if there is an intersection between two strings. Likely input it the splitRef data.
    Outputs True if no intersection
    overlapComp needs to delete paired reference number if being applied to the splits output
    """
    S1 = set(ref1)
    S2 = set(ref2)
    if len(S1.intersection(S2)) == 0:
        return True
    return False


def addSequenceList(input_file, multiFileFlag):

    """
    input_file is the file path to the fasta file.
    The function reads the fasta file into a dictionary with the format {proteinRef: protein}
    """

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    sequenceDictionary = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        name = name.split('|')[1]
        if multiFileFlag:
            fileName = input_file.split('/')[-1]
            name = fileName.split('.')[0] + '_' + name
            print(name)

        sequenceDictionary[name] = sequence
    return sequenceDictionary

def removeDupsQuick(seq, seqRef):

    seen = set()
    seen_add = seen.add
    initial = []
    initialRef = []
    # initial = [x for x in seq if not (x in seen or seen_add(x))]
    for i in range(0, len(seq)):
        if not (seq[i] in seen or seen_add(seq[i])):
            initial.append(seq[i])
            initialRef.append(seqRef[i])

    return initial, initialRef


def combMass(combine, combineRef, origProtTups = None):
    massDict = {}
    try:
        maxMass = mgfData.maxMass
        # print('in maxmass')
    except:
        maxMass = 1000000
    for i in range(0, len(combine)):
        totalMass = 0
        for j in range(0, len(combine[i])):
            totalMass += monoAminoMass[combine[i][j]]
        totalMass += H20_MASS
        if totalMass > maxMass:
            print(combine[i])
            continue
        if origProtTups == None:
            massRefPair = [totalMass, combineRef[i]]
            massDict[combine[i]] = massRefPair
        # when trans is being run, peptide which have already been added to massDict need their original peptide location
        # information updated so it is not lost in the running process.
        else:
            if combine[i] in massDict.keys():
                    massDict[combine[i]][2].append(origProtTups[i][0])
                    massDict[combine[i]][2].append(origProtTups[i][1])

            else:
                massRefPair = [totalMass, combineRef[i], origProtTups[i]]
                massDict[combine[i]] = massRefPair
    return massDict


def changeRefToDash(ref):
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
    for key, value in massDict.items():
        refNew = changeRefToDash(value[1])
        value[1] = refNew
    return massDict


def getFinalPath(outputPath, spliceType):
    outputPathSmall = str(outputPath)[0:-6]
    newPath = str(outputPathSmall) + '-' + spliceType + '.csv'
    return Path(newPath)

def nth_replace(string, old, new, n=1, option='only nth'):

    # https://stackoverflow.com/questions/35091557/replace-nth-occurrence-of-substring-in-string
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


def processLockInit(lockVar, toWriteQueue, mgfObj, childTable, linSetQueue):

    """
    Designed to set up a global lock for a child processes (child per protein)
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
    Designed to set up a global lock for a child processes (child per protein)
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


