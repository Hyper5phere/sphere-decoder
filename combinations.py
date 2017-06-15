# test script for the combinations algorithm
# (does not work)

def combinations(symbols, combination, dim):
    curr = []
    if dim == -1:
        return [combination]
    else:
        for s in symbols:
            combination[dim] = s
            curr += combinations(symbols, combination, dim-1)
        return curr

print(repr(combinations([1,2,3], [1,1,1], 2)))