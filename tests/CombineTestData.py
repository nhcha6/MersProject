splits = [['A', 'AB', 'ABC', 'B', 'BCD', 'BCDE']]
splitRef = [[[1], [1, 2], [1, 2, 3], [4], [4, 5, 6], [4, 5, 6, 7]]]

min2max7Combine = [['AAB', 'ABA', 'AABC', 'ABCA', 'AB', 'BA', 'ABCD', 'BCDA', 'ABCDE', 'BCDEA',
                    'ABABC', 'ABCAB', 'ABB', 'BAB', 'ABBCD', 'BCDAB', 'ABBCDE', 'BCDEAB',
                    'ABCB', 'BABC', 'ABCBCD', 'BCDABC', 'ABCBCDE', 'BCDEABC',
                    'BBCD', 'BCDB', 'BBCDE', 'BCDEB', 'BCDBCDE', 'BCDEBCD']]
min2max7CombineRef = [[[1, 1, 2], [1, 2, 1], [1, 1, 2, 3], [1, 2, 3, 1], [1, 4], [4, 1], [1, 4, 5, 6], [4, 5, 6, 1],
     [1, 4, 5, 6, 7], [4, 5, 6, 7, 1], [1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4], [4, 1, 2], [1, 2, 4, 5, 6],
     [4, 5, 6, 1, 2], [1, 2, 4, 5, 6, 7], [4, 5, 6, 7, 1, 2], [1, 2, 3, 4], [4, 1, 2, 3], [1, 2, 3, 4, 5, 6],
     [4, 5, 6, 1, 2, 3], [1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 1, 2, 3], [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7],
     [4, 5, 6, 7, 4], [4, 5, 6, 4, 5, 6, 7], [4, 5, 6, 7, 4, 5, 6]]]

min1max8Combine = [['AAB', 'ABA', 'AABC', 'ABCA', 'AB', 'BA', 'ABCD', 'BCDA', 'ABCDE', 'BCDEA',
                    'ABABC', 'ABCAB', 'ABB', 'BAB', 'ABBCD', 'BCDAB', 'ABBCDE', 'BCDEAB',
                    'ABCB', 'BABC', 'ABCBCD', 'BCDABC', 'ABCBCDE', 'BCDEABC',
                    'BBCD', 'BCDB', 'BBCDE', 'BCDEB', 'BCDBCDE', 'BCDEBCD']]
min1max8CombineRef = [[[1, 1, 2], [1, 2, 1], [1, 1, 2, 3], [1, 2, 3, 1], [1, 4], [4, 1], [1, 4, 5, 6], [4, 5, 6, 1],
     [1, 4, 5, 6, 7], [4, 5, 6, 7, 1], [1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4], [4, 1, 2], [1, 2, 4, 5, 6],
     [4, 5, 6, 1, 2], [1, 2, 4, 5, 6, 7], [4, 5, 6, 7, 1, 2], [1, 2, 3, 4], [4, 1, 2, 3], [1, 2, 3, 4, 5, 6],
     [4, 5, 6, 1, 2, 3], [1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 1, 2, 3], [4, 4, 5, 6], [4, 5, 6, 4],
     [4, 4, 5, 6, 7], [4, 5, 6, 7, 4], [4, 5, 6, 4, 5, 6, 7], [4, 5, 6, 7, 4, 5, 6]]]

min3max6Combine = [['AAB', 'ABA', 'AABC', 'ABCA', 'ABCD', 'BCDA', 'ABCDE', 'BCDEA',
                    'ABABC', 'ABCAB', 'ABB', 'BAB', 'ABBCD', 'BCDAB', 'ABBCDE', 'BCDEAB',
                    'ABCB', 'BABC', 'ABCBCD', 'BCDABC',
                    'BBCD', 'BCDB', 'BBCDE', 'BCDEB']]
min3max6CombineRef = [[[1, 1, 2], [1, 2, 1], [1, 1, 2, 3], [1, 2, 3, 1], [1, 4, 5, 6], [4, 5, 6, 1], [1, 4, 5, 6, 7],
     [4, 5, 6, 7, 1],[1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4], [4, 1, 2], [1, 2, 4, 5, 6], [4, 5, 6, 1, 2],
     [1, 2, 4, 5, 6, 7],[4, 5, 6, 7, 1, 2],[1, 2, 3, 4], [4, 1, 2, 3], [1, 2, 3, 4, 5, 6], [4, 5, 6, 1, 2, 3],
     [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7], [4, 5, 6, 7, 4]]]

min3max6Combine = [['AAB', 'ABA', 'AABC', 'ABCA', 'ABCD', 'BCDA', 'ABCDE', 'BCDEA',
                    'ABABC', 'ABCAB', 'ABB', 'BAB', 'ABBCD', 'BCDAB', 'ABBCDE', 'BCDEAB',
                    'ABCB', 'BABC', 'ABCBCD', 'BCDABC',
                    'BBCD', 'BCDB', 'BBCDE', 'BCDEB']]
min3max6CombineRef = [[[1, 1, 2], [1, 2, 1], [1, 1, 2, 3], [1, 2, 3, 1], [1, 4, 5, 6], [4, 5, 6, 1], [1, 4, 5, 6, 7],
     [4, 5, 6, 7, 1],[1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4], [4, 1, 2], [1, 2, 4, 5, 6], [4, 5, 6, 1, 2],
     [1, 2, 4, 5, 6, 7],[4, 5, 6, 7, 1, 2],[1, 2, 3, 4], [4, 1, 2, 3], [1, 2, 3, 4, 5, 6], [4, 5, 6, 1, 2, 3],
     [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7], [4, 5, 6, 7, 4]]]

min4max5Combine = [['AABC', 'ABCA', 'ABCD', 'BCDA', 'ABCDE', 'BCDEA',
                    'ABABC', 'ABCAB', 'ABBCD', 'BCDAB',
                    'ABCB', 'BABC', 'BBCD', 'BCDB', 'BBCDE', 'BCDEB']]
min4max5CombineRef = [[[1, 1, 2, 3], [1, 2, 3, 1], [1, 4, 5, 6], [4, 5, 6, 1], [1, 4, 5, 6, 7],
     [4, 5, 6, 7, 1], [1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4, 5, 6], [4, 5, 6, 1, 2],
     [1, 2, 3, 4], [4, 1, 2, 3], [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7], [4, 5, 6, 7, 4]]]

min5max4Combine = [[]]
min5max4CombineRef = [[]]

maxDist4Combine = [['AAB', 'ABA', 'AABC', 'ABCA',
                    'ABABC', 'ABCAB', 'ABB', 'BAB',
                    'ABCB', 'BABC',
                    'BBCD', 'BCDB', 'BBCDE', 'BCDEB']]
maxDist4CombineRef = [[[1, 1, 2], [1, 2, 1], [1, 1, 2, 3], [1, 2, 3, 1],
     [1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4], [4, 1, 2],
     [1, 2, 3, 4], [4, 1, 2, 3],
     [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7], [4, 5, 6, 7, 4]]]

maxDist3Combine = [['AABC', 'ABCA',
                    'ABABC', 'ABCAB',
                    'ABCB', 'BABC', 'BBCD', 'BCDB', 'BBCDE', 'BCDEB']]
maxDist3CombineRef = [[[1, 1, 2, 3], [1, 2, 3, 1],
     [1, 2, 1, 2, 3], [1, 2, 3, 1, 2],
     [1, 2, 3, 4], [4, 1, 2, 3], [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7], [4, 5, 6, 7, 4]]]

maxDist5Combine = [['AABC', 'ABCA', 'ABCD', 'BCDA',
                    'ABABC', 'ABCAB', 'ABBCD', 'BCDAB',
                    'ABCB', 'BABC', 'BBCD', 'BCDB', 'BBCDE', 'BCDEB']]
maxDist5CombineRef = [[[1, 1, 2, 3], [1, 2, 3, 1], [1, 4, 5, 6], [4, 5, 6, 1],
     [1, 2, 1, 2, 3], [1, 2, 3, 1, 2], [1, 2, 4, 5, 6], [4, 5, 6, 1, 2],
     [1, 2, 3, 4], [4, 1, 2, 3], [4, 4, 5, 6], [4, 5, 6, 4], [4, 4, 5, 6, 7], [4, 5, 6, 7, 4]]]

overlapTrueCombine = [['ABCD', 'BCDA', 'ABCDE', 'BCDEA',
                    'ABB', 'BAB', 'ABBCD', 'BCDAB', 'ABBCDE', 'BCDEAB',
                    'ABCB', 'BABC', 'ABCBCD', 'BCDABC']]
overlapTrueCombineRef = [[[1, 4, 5, 6], [4, 5, 6, 1], [1, 4, 5, 6, 7],
     [4, 5, 6, 7, 1], [1, 2, 4], [4, 1, 2], [1, 2, 4, 5, 6], [4, 5, 6, 1, 2],
     [1, 2, 4, 5, 6, 7],[4, 5, 6, 7, 1, 2],[1, 2, 3, 4], [4, 1, 2, 3], [1, 2, 3, 4, 5, 6], [4, 5, 6, 1, 2, 3]]]

overlapTrueMaxDist5Combine = [['ABCD', 'BCDA', 'ABBCD', 'BCDAB', 'ABCB', 'BABC']]
overlapTrueMaxDist5CombineRef = [[[1, 4, 5, 6], [4, 5, 6, 1], [1, 2, 4, 5, 6], [4, 5, 6, 1, 2],
                                  [1, 2, 3, 4], [4, 1, 2, 3]]]