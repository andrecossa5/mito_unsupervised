"""
Implement the logic for searching iteratively with a stop criterion.
"""

import time
import numpy as np
import pandas as pd
from itertools import chain


##


# Utils
def splitter(L):
    """
    Split a list L in chunks of variable n of elements
    """
    splitted = []
    i = 0
    go = True
    
    while go:
        n = np.random.choice([1,2,3])
        splitted.append(L[i:i+n])
        i += n
        if i >= len(L):
            go = False
    
    return splitted

# L = list(range(100))
# splitter(L, 15)[0]


##


def tester(g, n=1):
    """
    Simple test
    """

    test = np.sum([ x == 1 for x in g ]) >= n
    partition = g if test else []
    to_final = g if not test else []

    return test, partition, to_final

# L = np.random.choice([0,0,0,1,1], 20).tolist()
# L_splitted = splitter(L)
# [ tester(x) for x in L_splitted ]  


##


############################ Prova recursive
# def do_something(item):
#     print(item, np.mean(item))
# 
# def traverse_nested_list(nested_list):
#     for i, item in enumerate(nested_list):
#         if not all([ isinstance(x, list) for x in item ]):
#             do_something(item)
#         else:
#             traverse_nested_list(item) 
# 
# 
# # Usage
# nested_list = [[[1, 2], [3, 4]], [[[5], [6]], [7]], [[1, 3]], [[4], [5, 6, 7]], [[[[2]]]]]
# traverse_nested_list(nested_list)
############################ 


##


def main(L, n=2): 

    partitions = [L]
    final_partitions = []
    go = True
    i = 0

    while go:

        print(f'Iteration {i} partitions: {partitions}')
        splitted_partitions = [ splitter(x) for x in partitions ]
        # n_partitions = len(splitted_partitions) * len(splitted_partitions[0])
        # print(f'Begin iteration {i}: {n_partitions} partitions, {n} elements each.')

        print(f'Iteration {i} splitted partitions: {splitted_partitions}')
        tests, partitions, to_final = zip(
            *[ tester(x, n=n) for x in splitted_partitions[0] ]
        )

        splitted_partitions[0][0]

        partitions = [ x for x in partitions if x if len(x)>0 ]
        final_partitions.extend([ x for x in to_final if x if len(x)>0 ])
        
        print(f'Finished Iteration {i}: found {np.sum(tests)} partitions to re-split.')
        print('\n')

        go = any(tests)
        i += 1

    # final_partitions.extend([ x for x in splitted_partitions[0] ])
    unassigned = [ x for x in splitted_partitions[0] ]

    return final_partitions, unassigned


# Test
L = np.random.choice([0,0,0,1,1], 20).tolist()
final_partitions, unassigned = main(L, 2)



len(final_partitions)
len(L)

len(list(chain.from_iterable(final_partitions)))
len(L)


len(set(L) - set(chain.from_iterable(final_partitions)))




