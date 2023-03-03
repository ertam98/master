import sys

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def combine_lists(*lists):
    """Adapted from https://stackoverflow.com/a/1675380
    Takes lists and returns a set with all elements."""
    s = set(lists[0])
    for mylist in lists[1:]:
        s.update(mylist)
    return s

def geometric_mean(mylist):
    n = len(mylist)
    result = 1
    for element in mylist:
        result *= element

    return result**(1/n)

def myhash(sub, val, tol):
    temp = sorted(zip(sub, val))
    sub = tuple(element[0] for element in temp)
    val = tuple(round(tol*element[1]) for element in temp)
    
    return hash((sub, val))