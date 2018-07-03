# Unit Testing Files
import unittest
from Mers import *
from SplitTestData import *


# class for testing splits output with a false negative flag
class splitsKnownValuesFalse(unittest.TestCase):

    # Minimum value does not make a difference, the primary difference is made by max
    testMin1Max7 = {genericPeptide: [min1max12Splits, min1max12Ref]}

    testMin2Max2 = {genericPeptide: [min2max2Splits, min2max2Ref]}
    testMin2Max3 = {genericPeptide: [min2max3Splits, min2max3Ref]}
    testMin2Max4 = {genericPeptide: [min2max4Splits, min2max4Ref]}
    testMin2Max5 = {genericPeptide: [min2max5Splits, min2max5Ref]}
    testMin2Max6 = {genericPeptide: [min2max6Splits, min2max6Ref]}

    def testSplitsMin1Max12(self):
        for key, value in self.testMin1Max7.items():
            result, ref = splitDictPeptide(key, 1, 8, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)

    def testSplitsMin2Max2(self):
        for key, value in self.testMin2Max2.items():
            result, ref = splitDictPeptide(key, 2, 2, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)

    def testSplitsMin2Max3(self):
        for key, value in self.testMin2Max3.items():
            result, ref = splitDictPeptide(key, 2, 3, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)

    def testSplitsMin2Max4(self):
        for key, value in self.testMin2Max4.items():
            result, ref = splitDictPeptide(key, 2, 4, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)

    def testSplitsMin2Max5(self):
        for key, value in self.testMin2Max5.items():
            result, ref = splitDictPeptide(key, 2, 5, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)

    def testSplitsMin2Max6(self):
        for key, value in self.testMin2Max6.items():
            result, ref = splitDictPeptide(key, 2, 6, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)
    '''
    def testSplitsMin3Max7(self):
    
    def testSplitsMin2Max8(self):

    def testSplitsMin2Max9(self):

    def testSplitsMin2Max10(self):
    
    def testSplitsMin2Max11(self):

    '''
"""
# class for testing splits output with a true linear flag
class splitsKnownValuesTrue(unittest.TestCase):

     genericPeptideTrue = 'ABCDEFG'
     genericSplitsTrue = ['A','AB','ABC','ABCD','ABCDE','ABCDEF',
                    'ABCDEFG','B', 'BC', 'BCD', 'BCDE',
                    'BCDEF', 'BCDEFG', 'C', 'CD', 'CDE', 'CDEF',
                    'CDEFG', 'D', 'DE', 'DEF', 'DEFG', 'E', 'EF',
                    'EFG', 'F', 'FG', 'G']
     genericRefTrue = [[0], [0,1], [0,1,2], [0,1,2,3], [0,1,2,3,4],
                     [0,1,2,3,4,5], [0,1,2,3,4,5,6], [1], [1,2], [1,2,3], [1,2,3,4],
                     [1,2,3,4,5], [1,2,3,4,5,6], [2], [2,3], [2,3,4],
                     [2,3,4,5], [2,3,4,5,6], [3], [3,4],
                     [3,4,5], [3,4,5,6], [4],
                     [4,5], [4,5,6], [5], [5,6], [6]]

     min3max5SplitsTrue = ['ABC','ABCD', 'ABCDE', 'BCD', 'BCDE', 'BCDEF', 'CDE', 'CDEF',
                           'CDEFG', 'DEF', 'DEFG', 'EFG']
     min3max5RefTrue = [[0,1,2], [0,1,2,3], [0,1,2,3,4], [1,2,3], [1,2,3,4], [1,2,3,4,5], [2,3,4],
                        [2,3,4,5], [2,3,4,5,6], [3,4,5], [3,4,5,6], [4,5,6]]


     def testSplitsGenericCase(self):
        result, ref = splitDictPeptide(self.genericPeptideTrue, 1, 12, True)
        self.assertEqual(self.genericSplitsTrue, result)
        self.assertEqual(self.genericRefTrue, ref)

     def testSplitsMinMax(self):
         result, ref = splitDictPeptide(self.genericPeptideTrue, 3, 5, True)
         self.assertEqual(self.min3max5SplitsTrue, result)
         self.assertEqual(self.min3max5RefTrue, ref)


"""
if __name__ == "__main__":
    unittest.main()