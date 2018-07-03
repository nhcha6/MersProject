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

testMin1Max12 = {genericPeptide: [min1max12Splits, min1max12Ref]}
testMin3Max5 = {genericPeptide: [min3max5Splits, min3max5Ref]}