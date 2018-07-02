genericPeptide = 'ABCDEFG'

min1max12Splits = ['A', 'AB', 'ABC', 'ABCD', 'ABCDE', 'ABCDEF',
                 'ABCDEFG', 'B', 'BC', 'BCD', 'BCDE',
                 'BCDEF', 'BCDEFG', 'C', 'CD', 'CDE', 'CDEF',
                 'CDEFG', 'D', 'DE', 'DEF', 'DEFG', 'E', 'EF',
                 'EFG', 'F', 'FG', 'G']
min1max12Ref = [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4],
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

splitPeptideTest = [genericPeptide]
minedTest = [1, 3]
maxedTest = [12, 5]
splitTests = [[min1max12Splits, min3max5Splits]]
splitRefTests = [[min1max12Ref, min3max5Ref]]

# class for testing splits output with a false negative flag
class splitsKnownValues(unittest.TestCase):

     def testSplitsGenericCase(self):
         for i in range(0, len(splitPeptideTest)):
             for j in range(0, len(minedTest)):
                 result, ref = splitDictPeptide(splitPeptideTest[i], minedTest[j], maxedTest[j], False)
                 print(splitTests[i][j])
                 self.assertEqual(splitTests[i][j], result)
                 self.assertEqual(splitRefTests[i][j], ref)




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