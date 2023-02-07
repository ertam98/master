import sys

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def combine_lists(list1, list2):
    """Code from https://stackoverflow.com/a/1675380
    Takes two lists and returns a set with all elements."""
    s = set(list1)
    s.update(list2)
    return s