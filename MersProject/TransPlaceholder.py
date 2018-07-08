from Mers import *
def transOutput(finalPeptide, mined, maxed, overlapFlag, modList, outputPath, chargeFlags, linearFlag=False):
    print('output create trans')
    finalPath = str(outputPath) + '/Trans.csv'
    open(finalPath, 'w', newline='')
    splits, splitRef = splitDictPeptide(finalPeptide, mined, maxed, linearFlag)

    # combined eg: ['ABC', 'BCA', 'ACD', 'DCA']
    # combinedRef eg: [[0,1,2], [1,0,2], [0,2,3], [3,2,0]]
    # pass splits through combined overlap peptide and then delete all duplicates

    createTransProcess(splits, splitRef, mined, maxed, overlapFlag, modList, outputPath, chargeFlags)


def createTransProcess(splits, splitsRef, mined, maxed, overlapFlag, modList, outputPath, chargeFlags):

    print('create trans thread')
    combModless = []
    combModlessRef = []
    length = len(splits)


    # iterate through all of the splits creating a process for each split
    combProcesses = []
    startComb = time.time()
    print(len(splits))
    manager = multiprocessing.Manager()
    finalMassDict = manager.dict()
    for i in range(0, length):
        subsetSplits = splits[i:-1]
        subsetSplitsRef = splitsRef[i:-1]
        combinationProcess = multiprocessing.Process(target=specificTransProcess,
                                                     args=(subsetSplits, subsetSplitsRef, mined, maxed,
                                                           overlapFlag, modList, outputPath, chargeFlags, finalMassDict))
        combProcesses.append(combinationProcess)
        combinationProcess.start()



    for process in combProcesses:
        process.join()

    writeToCsv(finalMassDict, 'a', 'Trans', outputPath, 'Trans', chargeFlags)
    endComb = time.time()
    print('combination time: ' + str(endComb - startComb))


def specificTransProcess(subsetSplits, subSplitsRef, mined, maxed, overlapFlag, modList, outputPath, chargeFlags,
                         finalMassDict):
    print("Created trans process!!")
    subsetComb, subsetCombRef = combinePeptideTrans(subsetSplits, subSplitsRef, mined, maxed, overlapFlag)
    subsetComb, subsetCombRef = removeDupsQuick(subsetComb, subsetCombRef)
    massDict = combMass(subsetComb, subsetCombRef)
    massDict = applyMods(massDict, modList)
    chargeIonMass(massDict, chargeFlags)
    finalMassDict.update(massDict)
    print("Printed trans process to csv!")