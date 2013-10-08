import random
import pickle

def dict_sample(d):
    """
    Sample a key from a dictionary using the values as probabilities (unnormalized)
    """
    cutoff = random.random()
    normalizer = float(sum(d.values()))
    #print "Normalizer: ", normalizer

    current = 0
    for i in d:
        assert(d[i] > 0)
        current += float(d[i]) / normalizer
        if current >= cutoff:
            #print "Chose", i
            return i
    print "Didn't choose anything: ", cutoff, current

def arg_max(dictionary):
    best_val = max(dictionary.values())
    return [x for x in dictionary if dictionary[x] == best_val][0]

def levenshtein(a,b, add_value=1, delete_value=1, substitute_value=1):
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n

    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+add_value, current[j-1]+delete_value
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + substitute_value
            current[j] = min(add, delete, change)

    return current[n]

def reverse_dict(initial_dict):
    assert(type(initial_dict) == type({}))
    d = {}
    for key in initial_dict:
        d[initial_dict[key]] = key
    return d

def dictionary_sum(a, b):
    """
    Sums two dictionaries and returns the result.  Makes no promises about not
    modifying arguments.
    """
    for ii in b:
        if ii in a:
            a[ii] += b[ii]
        else:
            a[ii] = b[ii]
    return a

def count_line(counts):
    """
    Given a dictionary of words (indexed by number), return an LDA line.
    """
    s = str(len(counts))
    for ii in counts:
        s += " %i:%i" % (ii, counts[ii])
    s += "\n"
    return s

def flatten(iterable):
    """
    Given a list of lists (of lists ...), iterate through all of the elements.
    """

    it = iter(iterable)
    for e in it:
        if isinstance(e, (list, tuple)):
            for f in flatten(e):
                yield f
        else:
            yield e

def poll_iterator(iter_list, shuffle=True):
    """
    Given a list of iterable things, poll each of them one at a time
    until they run out, then continue polling from the others.
    """
    if shuffle:
        random.shuffle(iter_list)

    alive = [True] * len(iter_list)
    iters = [x.__iter__() for x in iter_list]
    num_iters = len(iters)
    while max(alive):
        for ii in xrange(num_iters):
            if alive[ii]:
                try:
                    yield iters[ii].next()
                except StopIteration:
                    alive[ii] = False
    return

def write_pickle(obj, filename, protocol=-1):
    """
    I can never remember the syntax
    """

    o = open(filename, 'wb')
    pickle.dump(obj, o, protocol)
    o.close()

def read_pickle(filename, protocol=-1):
    infile = open(filename, 'rb')
    val = pickle.load(infile)
    infile.close()
    return val
