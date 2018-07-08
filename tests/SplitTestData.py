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

min2max2Splits = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
min2max2Ref = [[0], [1], [2], [3], [4], [5], [6]]

min2max3Splits = ['A', 'AB', 'B', 'BC', 'C', 'CD',
                  'D', 'DE', 'E', 'EF', 'F', 'FG', 'G']
min2max3Ref = [[0], [0, 1], [1], [1, 2], [2], [2, 3],
               [3], [3, 4], [4], [4, 5], [5], [5, 6], [6]]

min2max4Splits =  ['A', 'AB', 'ABC', 'B', 'BC', 'BCD', 'C',
                   'CD', 'CDE', 'D', 'DE', 'DEF', 'E', 'EF',
                   'EFG', 'F', 'FG', 'G']
min2max4Ref = [[0], [0, 1], [0, 1, 2], [1], [1, 2], [1, 2, 3],
               [2], [2, 3], [2, 3, 4], [3], [3, 4],
               [3, 4, 5], [4], [4, 5], [4, 5, 6],
               [5], [5, 6], [6]]

min2max5Splits =  ['A', 'AB', 'ABC', 'ABCD', 'B', 'BC', 'BCD', 'BCDE',
                   'C', 'CD', 'CDE', 'CDEF', 'D', 'DE', 'DEF', 'DEFG',
                   'E', 'EF', 'EFG', 'F', 'FG', 'G']

min2max5Splits = ['A', 'AB', 'ABC', 'ABCD', 'B', 'BC', 'BCD', 'BCDE',
                  'C', 'CD', 'CDE', 'CDEF',
                  'D', 'DE', 'DEF', 'DEFG', 'E', 'EF',
                  'EFG', 'F', 'FG', 'G']
min2max5Ref = [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3],
               [1], [1, 2], [1, 2, 3], [1, 2, 3, 4],
               [2], [2, 3], [2, 3, 4],
               [2, 3, 4, 5], [3], [3, 4],
               [3, 4, 5], [3, 4, 5, 6], [4],
               [4, 5], [4, 5, 6], [5], [5, 6], [6]]

min2max6Splits = ['A', 'AB', 'ABC', 'ABCD', 'ABCDE',
                  'B', 'BC', 'BCD', 'BCDE','BCDEF', 'C',
                  'CD', 'CDE', 'CDEF','CDEFG', 'D', 'DE',
                  'DEF', 'DEFG', 'E', 'EF','EFG', 'F', 'FG', 'G']
min2max6Ref = [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4],
              [1], [1, 2], [1, 2, 3], [1, 2, 3, 4],
              [1, 2, 3, 4, 5], [2], [2, 3], [2, 3, 4],
              [2, 3, 4, 5], [2, 3, 4, 5, 6], [3], [3, 4],
              [3, 4, 5], [3, 4, 5, 6], [4],
              [4, 5], [4, 5, 6], [5], [5, 6], [6]]

min2max7Splits = ['A', 'AB', 'ABC', 'ABCD', 'ABCDE', 'ABCDEF',
                 'ABCDEFG', 'B', 'BC', 'BCD', 'BCDE',
                 'BCDEF', 'BCDEFG', 'C', 'CD', 'CDE', 'CDEF',
                 'CDEFG', 'D', 'DE', 'DEF', 'DEFG', 'E', 'EF',
                 'EFG', 'F', 'FG', 'G']



