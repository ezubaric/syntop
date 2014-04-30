if __name__ == "__main__":
    import sys

    tmp_tau_dir = sys.argv[1];
    tmp_word_dir = sys.argv[2];
    output_dir = sys.argv[3];

    '''
    input = open(tmp_word_dir, 'r');
    word_list = [];
    for line in input:
        word_list.append(line.strip());
    '''

    input = open(tmp_word_dir, 'r');
    word_list = {};
    for line in input:
        line = line.strip();
        content = line.split();
        #word_list[int(content[2])] = content[1];
        #word_list[len(word_list)] = content[2];
        #add by Thiago
        word_list[len(word_list)] = content[0];

    output = open(output_dir, 'w')
    f = open(tmp_tau_dir, 'r');
    for line in f:
        output.write("==============================================\n");

        k=0;
        tokens = line.split();

        assert len(tokens)==len(word_list)
        fd = [(float(tokens[t]), word_list[t]) for t in xrange(len(tokens))]
        fd.sort()
        fd.reverse()

        for (fd_value, fd_key) in fd:
            output.write(fd_key + "\t" + str(fd_value) + "\n");
        k+=1;
