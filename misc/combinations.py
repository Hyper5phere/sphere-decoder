# test script for the combinations algorithm

def comb_wrapper(symbset, dim):
    a = dim*[symbset[0]]
    comblist = []
    combinations(comblist, symbset, a, dim-1)
    return comblist

def combinations(comblist, symbols, combination, dim):
    if combination not in comblist:
        comblist.append(combination)

    if dim >= 0:
        for s in symbols:
            tmp = combination.copy()
            tmp[dim] = s
            combinations(comblist, symbols, tmp, dim-1)

symbset = [-3, -1, 1, 3]
code_length = 20

comblist = comb_wrapper(symbset, code_length)
print("All possible combinations for code vector:")
for i, c in enumerate(comblist):
    print(str(c) + "\t", end="")
    if (i+1) % code_length == 0:
        print("")