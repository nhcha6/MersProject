



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

def combinePeptideTrans(splits, splitRef, mined, maxed, overlapFlag):

    splitComb = []
    splitCombRef = []

    # toAdd variables hold temporary combinations for insertion in final matrix if it meets criteria
    toAddForward = ""
    # addForwardRef = []
    toAddReverse = ""
    # addReverseRef = []
    for j in range(0, len(splits)):
        # create forward combination of i and j
        toAddForward += splits[0]
        toAddForward += splits[j]
        addForwardRef = splitRef[0] + splitRef[j]
        toAddReverse += splits[j]
        toAddReverse += splits[0]
        addReverseRef = splitRef[j] + splitRef[0]

        # max, min and max distance checks combined into one function for clarity for clarity
        if combineCheck(toAddForward, mined, maxed, splitRef[0], splitRef[j]):
            if overlapFlag:
                if overlapComp(splitRef[0], splitRef[j]):
                    splitComb.append(toAddForward)
                    splitComb.append(toAddReverse)
                    splitCombRef.append(addForwardRef)
                    splitCombRef.append(addReverseRef)
            else:
                splitComb.append(toAddForward)
                splitComb.append(toAddReverse)
                splitCombRef.append(addForwardRef)
                splitCombRef.append(addReverseRef)

        toAddForward = ""
        toAddReverse = ""
        # addForwardRef = []
        # addReverseRef = []
    return splitComb, splitCombRef

