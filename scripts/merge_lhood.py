from glob import glob
import sys

if __name__ == "__main__":
    d = {}

    path = sys.argv[1]

    for ii in glob("%s/*.lhood" % path):
        print ii
        if ii.endswith("current.lhood"):
            continue

        iter = int(ii.rsplit("/", 1)[-1].split(".", 1)[0])

        d[iter] = float(open(ii).read().split()[1])

    o = open("%s/current.lhood" % path, 'w')
    for ii in sorted(d):
        o.write("%i\t%f\n" % (ii, d[ii]))
