
if __name__ == "__main__":
    import sys

    K = int(sys.argv[1])	# number of topics
    D = int(sys.argv[2])	# number of documents
    V = int(sys.argv[3])	# vocabulary size
    I = int(sys.argv[4])	# iteration number
    dataset = sys.argv[5]	# dataset dir
    suffix = sys.argv[6]	# suffix

    tmp_nu_dir = dataset + "_params_" + suffix + "/tmp-nu-"
    tmp_tau_dir = dataset + "_params_" + suffix + "/tmp-tau-"
    tmp_gamma_dir = dataset + "_params_" + suffix + "/tmp-gamma.txt"

    output = open(dataset + "_params_" + suffix + "/" + str(I) + ".nu", 'w')
    for k in xrange(K + 1):
        f = open(tmp_nu_dir + str(k) + ".txt", 'r')
        output_string = range(K)
        for line in f:
            tokens = line.split()
            indexs = tokens[0].split("_")
            output_string[int(indexs[2])] = tokens[1]
        output.write("\t".join(output_string) + "\n")
    print "output file", str(I) + ".nu to file system"

    output = open(dataset + "_params_" + suffix + "/" + str(I) + ".tau", 'w')
    for k in xrange(K):
        f = open(tmp_tau_dir + str(k) + ".txt", 'r')
        output_string = range(V)
        for line in f:
            tokens = line.split()
            indexs = tokens[0].split("_")
            output_string[int(indexs[2])] = tokens[1]
        output.write("\t".join(output_string) + "\n")
    print "output file", str(I) + ".tau to file system"

    output = open(dataset + "_params_" + suffix + "/" + str(I) + ".gamma", 'w')
    f = open(tmp_gamma_dir, 'r')
    output_gamma = {}
    for d in range(D):
        output_gamma[d] = range(K)
    for line in f:
        tokens = line.split()
        indexs = tokens[0].split("_")
        output_gamma[int(indexs[1])][int(indexs[2])] = tokens[1]
    for d in range(D):
        output.write("\t".join(output_gamma[d]) + "\n")
    print "output file", str(I) + ".gamma to file system"
