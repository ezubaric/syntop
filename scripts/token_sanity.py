
from collections import defaultdict
import sys

if __name__ == "__main__":
    input_file = sys.argv[1]
    topics = set()
    sums = defaultdict(dict)

    for ii in open(input_file):
        if not "\t" in ii:
            continue
        key, val = ii.split()
        fields = key.split("_")
        key = fields[0]

        print(fields)
        if len(fields) > 3:
            key += fields[-1]
            fields = map(int, fields[1:-1])
        else:
            fields = map(int, fields[1:])
        print(key, fields)

        if fields[1] != -1:
            sums[key][fields[0]] = sums[key].get(fields[0], 0.0) + float(val)
            if not key.startswith("gamma"):
                topics.add(fields[0])

    for ii in topics:
        for jj in sums:
            print jj, ii, sums[jj].get(ii, 0.0)
