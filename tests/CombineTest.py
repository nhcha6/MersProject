import unittest
from Mers import *
from CombineTestData import *

class combinePeptideKnownResult(unittest.TestCase):

    # overlapFlag = False, maxDistance = 'None'
    def testCombineMin2Max7(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 2, 7, False, 'None')
            self.assertEqual(combine, min2max7Combine[i], "Combine incorrect for min 2 max 7")
            self.assertEqual(combineRef, min2max7CombineRef[i], "CombineRef incorrect for min 2 max 7")

    # overlapFlag = False, maxDistance = 'None'
    def testCombineMin1Max8(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 1, 8, False, 'None')
            self.assertEqual(combine, min1max8Combine[i], "Combine incorrect for min 1 max 8")
            self.assertEqual(combineRef, min1max8CombineRef[i], "CombineRef incorrect for min 1 max 8")

    # overlapFlag = False, maxDistance = 'None'
    def testCombineMin3Max6(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 3, 6, False, 'None')
            self.assertEqual(combine, min3max6Combine[i], "Combine incorrect for min 3 max 6")
            self.assertEqual(combineRef, min3max6CombineRef[i], "CombineRef incorrect for min 3 max 6")

    # overlapFlag = False, maxDistance = 'None'
    def testCombineMin4Max5(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 4, 5, False, 'None')
            self.assertEqual(combine, min4max5Combine[i], "Combine incorrect for min 4 max 5")
            self.assertEqual(combineRef, min4max5CombineRef[i], "CombineRef incorrect for min 4 max 5")

    # overlapFlag = False, maxDistance = 'None'
    def testCombineMin5Max4(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 5, 4, False, 'None')
            self.assertEqual(combine, min5max4Combine[i], "Combine incorrect for min 5 max 4")
            self.assertEqual(combineRef, min5max4CombineRef[i], "CombineRef incorrect for min 5 max 4")

    # min = 3 max = 6 overlapFlag = False
    def testCombineMaxDist4(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 3, 6, False, 4)
            self.assertEqual(combine, maxDist4Combine[i], "Combine incorrect for maxDist = 4")
            self.assertEqual(combineRef, maxDist4CombineRef[i], "CombineRef incorrect for maxDist = 4")

    # min = 4 max = 5 max overlapFlag = False
    def testCombineMaxDist3(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 4, 5, False, 3)
            self.assertEqual(combine, maxDist3Combine[i], "Combine incorrect for maxDist = 3")
            self.assertEqual(combineRef, maxDist3CombineRef[i], "CombineRef incorrect for maxDist = 3")

    # min = 4 max = 5 max overlapFlag = False
    def testCombineMaxDist5(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 4, 5, False, 5)
            self.assertEqual(combine, maxDist5Combine[i], "Combine incorrect for maxDist = 5")
            self.assertEqual(combineRef, maxDist5CombineRef[i], "CombineRef incorrect for maxDist = 5")

    # min = 3 max = 6 maxDist = 'None'
    def testCombineOverlapTrue(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 3, 6, True, 'None')
            self.assertEqual(combine, overlapTrueCombine[i], "Combine incorrect for min 3 max 6")
            self.assertEqual(combineRef, overlapTrueCombineRef[i], "CombineRef incorrect for overlap = True")

    # min = 4 max = 5 maxDistance = 5
    def testCombineOverlapTrueMaxDist5(self):
        for i in range(0,len(splits)):
            combine, combineRef = combineOverlapPeptide(splits[i], splitRef[i], 4, 5, True, 5)
            self.assertEqual(combine, overlapTrueMaxDist5Combine[i], "Combine incorrect for overlap = True Maxdist = 5")
            self.assertEqual(combineRef, overlapTrueMaxDist5CombineRef[i], "CombineRef incorrect for overlap = True Maxdist = 5")

if __name__ == "__main__":
    unittest.main()