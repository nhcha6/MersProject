genericPeptide = 'ABCDEFG'
genericSplits = ['A', 'AB', 'ABC', 'ABCD', 'ABCDE', 'ABCDEF',
                 'ABCDEFG', 'B', 'BC', 'BCD', 'BCDE',
                 'BCDEF', 'BCDEFG', 'C', 'CD', 'CDE', 'CDEF',
                 'CDEFG', 'D', 'DE', 'DEF', 'DEFG', 'E', 'EF',
                 'EFG', 'F', 'FG', 'G']
genericRef = [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4],
              [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5, 6], [1], [1, 2], [1, 2, 3], [1, 2, 3, 4],
              [1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6], [2], [2, 3], [2, 3, 4],
              [2, 3, 4, 5], [2, 3, 4, 5, 6], [3], [3, 4],
              [3, 4, 5], [3, 4, 5, 6], [4],
              [4, 5], [4, 5, 6], [5], [5, 6], [6]]

min3max5Splits = ['A', 'AB', 'ABC', 'ABCD', 'B', 'BC', 'BCD', 'BCDE',
                  'C', 'CD', 'CDE', 'CDEF',
                  'D', 'DE', 'DEF', 'DEFG', 'E', 'EF',
                  'EFG', 'F', 'FG', 'G']
min3max5Ref = [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3],
               [1], [1, 2], [1, 2, 3], [1, 2, 3, 4],
               [2], [2, 3], [2, 3, 4],
               [2, 3, 4, 5], [3], [3, 4],
               [3, 4, 5], [3, 4, 5, 6], [4],
               [4, 5], [4, 5, 6], [5], [5, 6], [6]]

# Unit Testing Files
import unittest
from Mers import *

splitPeptideTest = []
splitTest = []
splitRefTests = []

splitPeptideTest.append(genericPeptide)
splitTest.append(genericSplits)
splitRefTests.append(genericRef)


# class for testing splits output with a false negative flag
class splitsKnownValues(unittest.TestCase):



     def testSplitsGenericCase(self):
         
        result, ref = splitDictPeptide(self.genericPeptide, 1, 12, False)
        self.assertEqual(self.genericSplits, result)
        self.assertEqual(self.genericRef, ref)

     def testSplitsMinMax(self):
         result, ref = splitDictPeptide(self.genericPeptide, 3, 5, False)
         self.assertEqual(self.min3max5Splits, result)
         self.assertEqual(self.min3max5Ref, ref)

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



if __name__ == "__main__":
    unittest.main()