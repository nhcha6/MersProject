# Unit Testing Files
import unittest
from Mers import *
from SplitTestData import *


# class for testing splits output with a false negative flag
class splitsKnownValuesFalse(unittest.TestCase):

    def testSplitsMin1Max12(self):
        for key, value in testMin1Max12.items():
            result, ref = splitDictPeptide(key, 1, 12, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)


    def testSplitsMin2Max2(self):

    def testSplitsMin2Max3(self):

    def testSplitsMin2Max4(self):

    def testSplitsMin2Max5(self):

    def testSplitsMin2Max6(self):

    def testSplitsMin2Max7(self):

    def testSplitsMin2Max8(self):

    def testSplitsMin2Max9(self):

    def testSplitsMin2Max12(self):
    
    def testSplitsMin3Max4(self):

    def testSplitsMin3Max5(self):
        for key, value in testMin3Max5.items():
            result, ref = splitDictPeptide(key, 3, 5, False)
            self.assertEqual(value[0], result)
            self.assertEqual(value[1], ref)
    
    def testSplitsMin3Max6(self):
    
    def testSplitsMin3Max7(self):
    
    def testSplitsMin3Max8(self):
    
    def testSplitsMin3Max9(self):
    
    def testSplitsMin3Max12(self):
    
    def testSplitsMin4Max5(self):
    
    def testSplitsMin4Max6(self):

    def testSplitsMin4Max7(self):

    def testSplitsMin4Max8(self):

    def testSplitsMin4Max9(self):

    def testSplitsMin4Max10(self):

    def testSplitsMin4Max11(self):

    def testSplitsMin4Max12(self):

    def testSplitsMin5Max6(self):

    def testSplitsMin5Max7(self):

    def testSplitsMin5Max8(self):

    def testSplitsMin5Max9(self):

    def testSplitsMin5Max10(self):

    def testSplitsMin5Max11(self):

    def testSplitsMin5Max12(self):

    def testSplitsMin6Max7(self):

    def testSplitsMin6Max8(self):

    def testSplitsMin6Max9(self):

    def testSplitsMin6Max10(self):

    def testSplitsMin6Max11(self):

    def testSplitsMin6Max12(self):

    def testSplitsMin7Max8(self):

    def testSplitsMin7Max9(self):

    def testSplitsMin7Max10(self):

    def testSplitsMin7Max11(self):

    def testSplitsMin7Max12(self):

    def testSplitsMin8Max9(self):

    def testSplitsMin8Max10(self):

    def testSplitsMin8Max11(self):

    def testSplitsMin8Max12(self):

    def testSplitsMin9Max10(self):

    def testSplitsMin9Max12(self):

    def testSplitsMin10Max11(self):

    def testSplitsMin10Max12(self):

    def testSplitsMin11Max12(self):

    def testSplitsMin12Max12(self):


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