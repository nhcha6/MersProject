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


TRANS = "Trans"
LINEAR = "Linear"
CIS = "Cis"

MEMORY_THRESHOLD = 20
MEMORY_THRESHOLD_COMBINE = 90
NUM_PROC_TOTAL = 1000
MAX_PROC_ALIVE = 40
MEMFLAG = 'mem'
STOPFLAG = 'stop'
PROC_FINISHED = 'Process Finished'
TEMP_WRITTEN = "temp written"

logging.basicConfig(level=logging.DEBUG, format='%(message)s')
# logging.disable(logging.INFO)

mgfData = None

# massDict syntax: PEPTIDE: [monoisotopic mass, [referenceLocation], {charge: m/z}]
# for example: {'DQQ': [389.15466499999997, ['121', '117-118'], {2: 195.58527249999997}],
# 'QDQ': [389.15466499999997, ['118-119', '117'], {2: 195.58527249999997}]......


class Fasta:

    """
    Class that represents the input from a fasta file
    """

    def __init__(self, inputFile):

        self.inputFile = [inputFile]
        self.allProcessList = []
        self.procGenCounter = 0
        self.pepCompleted = multiprocessing.Queue()
        self.completedProcs = 0
        self.totalProcs = 0

    def generateOutput(self, mined, maxed, overlapFlag, transFlag, cisFlag, linearFlag, csvFlag, pepToProtFlag,
                       protToPepFlag, modList, maxMod, maxDistance, outputPath, chargeFlags, mgfObj, mgfFlag):

        """
        Function that literally combines everything to generate output
        """
        # calculate the total number of processes that will be generated (NUM_PROC_TOTAL for every splice type running)
        if transFlag:
            self.totalProcs += NUM_PROC_TOTAL
        if cisFlag:
            self.totalProcs += NUM_PROC_TOTAL
        if linearFlag:
            self.totalProcs += NUM_PROC_TOTAL

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

        finalPath = None

        # Open the csv file if the csv file is selected
        if csvFlag:
            finalPath = getFinalPath(outputPath, spliceType)
            open(finalPath, 'w')

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

        # if len(seqDict) <= 1:
        #     logging.info('Only 1 protein, therefore trans not relevant')
        #     return

        finalPeptide, protIndexList, protList = combinePeptides(seqDict)

        splits, splitRef = splitTransPeptide(spliceType, finalPeptide, mined, maxed, protIndexList)

        splitLen = len(splits)

        # configure mutliprocessing functionality
        num_workers = multiprocessing.cpu_count()

        # Used to lock write access to file
        lockVar = multiprocessing.Lock()

        toWriteQueue = multiprocessing.Queue()
        linCisQueue = multiprocessing.Queue()

        pool = multiprocessing.Pool(processes=num_workers, initializer=processLockTrans, initargs=(lockVar, toWriteQueue,
                                                                                                   splits,splitRef, mgfObj,
                                                                                                   modTable, linCisQueue))

        writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, linCisQueue, pepToProtFlag,
                                                                     protToPepFlag, self.pepCompleted, True))
        writerProcess.start()

        # Create a process for pairs of splits, pairing element 0 with -1, 1 with -2 and so on.
        splitsIndex = []
        procSize = math.ceil(splitLen / (NUM_PROC_TOTAL*2))

        # if procSize = 1, NUM_PROC_TOTAL may be much larger than the number of procs actually generated. Thus, to
        # ensure the progress bar still works readjust self.totalProcs to be splitLen/2 instead.
        if procSize == 1:
            self.totalProcs = self.totalProcs - NUM_PROC_TOTAL + math.ceil(splitLen/2)

        #maxMem = psutil.virtual_memory()[1] / 2

        for i in range(0, math.ceil(splitLen / 2), procSize):
            if i + procSize > math.floor(splitLen / 2):
                for j in range(i, splitLen - i):
                    splitsIndex.append(j)
            else:
                for j in range(i, i + procSize):
                    splitsIndex.append(j)
                    splitsIndex.append(splitLen - 1 - j)

            # wait until the number processes started but not completed is less than 40 to begin a new process.
            self.procGenCounter += 1
            while True:
                if self.procGenCounter - self.completedProcs < MAX_PROC_ALIVE:
                    break

            pool.apply_async(transProcess, args=(splitsIndex, mined, maxed, modList, maxMod,
                                                 finalPath, chargeFlags, mgfFlag, csvFlag, protIndexList, protList))
            splitsIndex = []

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

        pool.close()
        pool.join()

        toWriteQueue.put(STOPFLAG)
        writerProcess.join()
        logging.info("All " + spliceType + " !joined")

    def cisAndLinearOutput(self, inputFile, spliceType, mined, maxed, overlapFlag, csvFlag, pepToProtFlag, protToPepFlag,
                           modList, maxMod, maxDistance, outputPath, chargeFlags, mgfObj, childTable, mgfFlag):

        """
        Process that is in charge for dealing with cis and linear. Creates sub processes for every protein to compute
        their respective output
        """

        finalPath = None

        # Open the csv file if the csv file is selected
        if csvFlag:
            finalPath = getFinalPath(outputPath, spliceType)
            open(finalPath, 'w')

        num_workers = multiprocessing.cpu_count()

        # Used to lock write access to file
        lockVar = multiprocessing.Lock()

        toWriteQueue = multiprocessing.Queue()
        linSetQueue = multiprocessing.Queue()
        pool = multiprocessing.Pool(processes=num_workers, initializer=processLockInit,
                                    initargs=(lockVar, toWriteQueue, mgfObj, childTable, linSetQueue))
        writerProcess = multiprocessing.Process(target=writer,
                                                args=(toWriteQueue, outputPath, linSetQueue, pepToProtFlag,
                                                      protToPepFlag, self.pepCompleted))
        writerProcess.start()

        # maxMem = psutil.virtual_memory()[1] / 2

        # calculate total size of input fasta
        totalProt = 0
        for file in inputFile:
            with open(file, "rU") as handle:
                for entry in SeqIO.parse(handle, 'fasta'):
                    totalProt += 1

        pepPerProc = math.ceil(totalProt / NUM_PROC_TOTAL)
        print("Process Size: " + str(pepPerProc))
        # if pepPerProc == 1, it is likely there are less proteins than NUM_PROC_TOTAL. Therefore, we adjust
        # the self.totalProcs counter (used in the progress bar).
        if pepPerProc == 1:
            self.totalProcs = self.totalProcs - NUM_PROC_TOTAL + totalProt


        # initialise counter variables
        counter = 0
        protDict = {}
        for file in inputFile:
            with open(file, "rU") as handle:
                for record in SeqIO.parse(handle, 'fasta'):
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

                    if counter % pepPerProc == 0:
                        self.procGenCounter += 1
                        # if more than the max number of processes has been generated, wait for a process
                        # to finish before another one is started
                        while True:
                            if self.procGenCounter - self.completedProcs < MAX_PROC_ALIVE:
                                break
                        # logging.info(spliceType + " process started for: " + seq[0:5])
                        # Start the processes for each protein with the targe function being genMassDict
                        pool.apply_async(genMassDict, args=(spliceType, protDict, mined, maxed, overlapFlag,
                                                            csvFlag, modList, maxMod, maxDistance, finalPath,
                                                            chargeFlags, mgfFlag))
                        protDict = {}

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
                            pool = multiprocessing.Pool(processes=num_workers, initializer=processLockInit,
                                                        initargs=(lockVar, toWriteQueue, mgfObj, childTable,
                                                                  linSetQueue))

        # add left over to pool if there is any
        if protDict:
            self.procGenCounter += 1
            pool.apply_async(genMassDict, args=(spliceType, protDict, mined, maxed, overlapFlag,
                                                csvFlag, modList, maxMod, maxDistance, finalPath,
                                                chargeFlags, mgfFlag))

        pool.close()
        pool.join()

        toWriteQueue.put(STOPFLAG)
        writerProcess.join()
        logging.info("All " + spliceType + " !joined")

# takes splits index from the multiprocessing pool and adds to writer the output. Splits and SplitRef are global
# variables within the pool.
def transProcess(splitsIndex, mined, maxed, modList, maxMod, finalPath,
                 chargeFlags, mgfFlag, csvFlag, protIndexList, protList):


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

# Only works if we presume Cis proteins aren't being created in the trans process.
def findOrigProt(combinedRef, protIndexList, protList):
    proteinTups = []
    for i in range(0, len(combinedRef)):
        protRef1 = ""
        protRef2 = ""
        ref = combinedRef[i]
        protIndex1, protIter1 = findInitProt(ref[0] - 1, protIndexList)
        #print(protIndex1)
        prot1 = protList[protIter1]

        # special check if peptide is ovelap spliced
        if len(set(ref)) != len(ref):
            protRef1 += ('(' + str(min(ref) - protIndex1[0]))
            protRef1 += ('-' + str(max(ref) - protIndex1[0]) + ')')
            proteinTups.append([(prot1, protRef1),('Overlap',"")])
            continue

        for j in range(1,len(ref)):
            #print(j)
            if ref[j] - 1 > protIndex1[1] or ref[j] - 1 < protIndex1[0]:
                #check to see if the first split is at least 6 amino acids in length.
                # if so append the location of the split within the peptide to prot1
                if j > 5:
                    protRef1 += ('(' + str(ref[0] - protIndex1[0]))
                    protRef1 += ('-' + str(ref[j-1] - protIndex1[0]) + ')')

                protIndex2, protIter2 = findInitProt(ref[j] - 1, protIndexList)
                prot2 = protList[protIter2]
                # same as above, check if second split is at least 6 amino acids long
                if len(ref) - j > 5:
                    protRef2 += ('(' + str(ref[j] - protIndex2[0]))
                    protRef2 += ('-' + str(ref[-1] - protIndex2[0]) + ')')

                proteinTups.append([(prot1,protRef1),(prot2,protRef2)])
                # combinedRef[i].insert(j, prot2)
                # combinedRef[i].insert(0,prot1)
                break
    return proteinTups


def findInitProt(index, protIndexList):
    length = protIndexList[-1][-1]
    #print(length)
    # plus 1 needed for when Protein length is perfectly divisible by protein index length
    aveLen = math.ceil(length/len(protIndexList)) + 1
    #print(aveLen)
    protIter = math.floor(index/aveLen)
    #print(protIndexList[protIter][0])
    if protIter == len(protIndexList):
        protIter -= 1
        #print(protIter)
    while True:
        lower = protIndexList[protIter][0]
        upper = protIndexList[protIter][1]
        if lower <= index:
            if upper >= index:
                #print(protIndexList[protIter])
                return protIndexList[protIter], protIter
            else:
                protIter += 1
        else:
            protIter -= 1

def splitTransPeptide(spliceType, peptide, mined, maxed, protIndexList):

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


def combineTransPeptide(splits, splitRef, mined, maxed, splitsIndex, protIndexList):

    """
    Input: splits: list of splits, splitRef: list of the character indexes for splits, mined/maxed: min and max
    size requirements, overlapFlag: boolean value true if overlapping combinations are undesired.
    Output: all combinations of possible splits which meets criteria
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

def genMassDict(spliceType, protDict, mined, maxed, overlapFlag, csvFlag, modList, maxMod,
                maxDistance, finalPath, chargeFlags, mgfFlag):

    """
    Compute the peptides for the given protein
    """
    try:
        for protId, protSeq in protDict.items():
            start = time.time()

            # Get the initial peptides and their positions, and the set of linear peptides produced for this protein
            combined, combinedRef, linSet = outputCreate(spliceType, protSeq, mined, maxed, overlapFlag, maxDistance)

            # add this set of linear proteins to the linProt queue
            genMassDict.linSetQueue.put(linSet)

            # Convert it into a dictionary that has a mass
            massDict = combMass(combined, combinedRef)

            # Apply mods to the dictionary values and update the dictionary
            massDict = applyMods(massDict, modList, maxMod)

            # Add the charge information along with their masses
            massDict = chargeIonMass(massDict, chargeFlags)

            # Get the positions in range form, instead of individuals (0,1,2) -> (0-2)
            massDict = editRefMassDict(massDict)

            if mgfFlag:
                #allPeptides = getAllPep(massDict)
                allPeptides = massDict.keys()
                allPeptidesDict = {}
                for peptide in allPeptides:
                    allPeptidesDict[peptide] = protId
                genMassDict.toWriteQueue.put((allPeptidesDict,False))

            # If there is an mgf file AND there is a charge selected
            elif mgfData is not None and True in chargeFlags:
                #fulfillPpmReq(mgfObj, massDict)
                matchedPeptides, modCountDict = generateMGFList(protId, mgfData, massDict, modList)
                genMassDict.toWriteQueue.put((matchedPeptides, modCountDict))

            # If csv is selected, write to csv file
            if csvFlag:
                logging.info("Writing locked :(")
                lock.acquire()

                writeToCsv(massDict, protId, finalPath, chargeFlags)
                lock.release()
                logging.info("Writing released!")

        # put PROC_FINISHED flag to toWriteQueue
        genMassDict.toWriteQueue.put(PROC_FINISHED)

    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in protein: ' + protId + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e

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
    memThreshFlag = False

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
            # update memThreshFlag so that linCisPeptides can be removed from the file at the end
            memThreshFlag = True
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

    # if memThreshFlag is True, we have written a tempFile. We must therefore go back and check that all peptides
    # in linCisSet haven't made their way into the these outputs.
    if memThreshFlag:
        # check that there is anything in linCisSet to compare to (for linear splicing it will be empty)
        if linCisSet:
            # we now want to read each output file into memory, delete peptides in linCisQueue, and write it
            # to file again.
            remFinalCisLin(linCisSet, saveHandle, fileCount)

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

def remFinalCisLin(linCisSet, saveHandle, fileCount):
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


def combinePeptides(seqDict):

    """
    combines an array of strings into one string. Used for ultimately segments from multiple peptides
    """

    dictlist = []
    protIndexList = []
    protList = []
    ind = 0
    for key, value in seqDict.items():
        dictlist.append(value)
        protIndexList.append([ind,ind + len(value) - 1])
        protList.append(key)
        ind += len(value)

    # print(protIndexList)

    finalPeptide = ''.join(dictlist)
    return finalPeptide, protIndexList, protList


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


