# Adapted from https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
# removes duplicates given a list
def removeDups(seq):
    seen = set()
    seen_add = seen.add
    initial = [x for x in seq if not (x in seen or seen_add(x))]
    temp = []
    #Remove permutations
    for i in range(0, len(initial)):
        for j in range(i + 1, len(initial)):
            if (combRemove(initial[i], initial[j])):
                temp.append(initial[j])

    final = [x for x in initial if x not in temp]
    return final

"""
Returns true if two strings are a permutation of each other -> Called in removeDups
"""


def combRemove(ref1, ref2):
    return (len(ref1) == len(ref2) and sorted(ref1) == sorted(ref2))