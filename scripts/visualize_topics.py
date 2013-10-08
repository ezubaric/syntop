import sys
from string import strip
from time import sleep
from math import exp

def read_beta(beta_file, topics):
    try:
        beta = open(beta_file).read().split()
    except IOError:
        print "Beta not found"
        for i in xrange(len(topics)):
            topics[i].add_beta(1.0)
        return
    for i in xrange(len(beta)):
        topics[i].add_beta(float(beta[i]))

def read_gamma(gamma_file, topics):
    line = 0
    num_topics = len(topics)
    try:
        input = open(gamma_file, "r")
        for line in input:
            vals = []
            tokens = line.split()
            for token in tokens:
                vals.append(token)
        vals = map(float, vals)
    except IOError:
        print "Gamma not found"
        for i in xrange(num_topics):
            topics[i].add_gamma(1.0 / float(num_topics))
        return 

    scale = max(vals)
    scale += 0.001

    for i in xrange(len(vals)):
        z = i % num_topics
        topics[z].add_gamma(vals[i] / scale)
        

def read_nu(nu_file, cutoff_multiplier, topics):
    # There is an extra topic for the start state
    num_topics = len(topics) - 1
    cutoff = cutoff_multiplier / float(num_topics)

    try:
        infile = open(nu_file)
    except IOError:
        return 
    # We include the extra start state here
    for i in xrange(num_topics + 1):
        line = []
        tokens = infile.readline().split()
        sum = 0.0
        for j in xrange(num_topics):
            val = float(tokens[j])
            sum += val
            line.append([val, j])
        line = map(lambda x: [x[1], x[0]], line)
        line = filter(lambda x: x[1]/sum > cutoff, line)
        topics[i].set_transitions([sum, line])

def read_tau(tau_file, wordlist, num_words, topics):
    infile = open(tau_file)
    # There's one extra "topic" for the starting state
    for z in xrange(len(topics)-1):
        vals = []
        tokens = infile.readline().split()
        # Since the file is one entry per line, we need to read a
        # 'row' like this
        for token in tokens:
            vals.append(float(token))
        word_strength = map(lambda x,y,z: [float(x),y + " %#.4g" % exp(float(x)),z], vals, wordlist, range(len(wordlist)))
        if len(wordlist) > 10000:
            word_strength = filter(lambda x: x[2] > 50 and x[2] < 10000, word_strength)
        elif len(wordlist) > 1000:
            word_strength = filter(lambda x: x[2] > 100, word_strength)
        word_strength = map(lambda x: x[:2], word_strength)
        word_strength.sort(lambda x,y: cmp(y,x))

        topics[z].set_words(word_strength[:num_words])

class Topic:

    def __init__(self, id):
        self.id_ = id
        self.gamma_ = 0
        self.beta_ = 0
        self.tau_ = {}
        self.nu_ = []

    def set_words(self, words):
        self.tau_ = words

    def set_transitions(self, links):
        self.nu_sum_ = links[0]                              
        self.nu_ = links[1]

    def add_gamma(self, x):
        self.gamma_ = max(self.gamma_, x)

    def add_beta(self, x):
        self.beta_ = x

    def graphviz(self):
        s = ""
        s += "\"topic" + str(self.id_) + "\" [\n"
        s += "   label = \"Topic " + str(self.id_) + ":" + str(self.beta_)
        for i in self.tau_:
            s += "| " + i[1] #+ ":" + str(i[0]) + " "
        s += "\"\n"
        s += "   shape = \"record\"\n"
        s += "   style = filled\n"
        s += "   fillcolor = grey" + str(int(50 + 50 * self.gamma_)) + "\n"
        s += "];\n"
        for i in self.nu_:
            s += "\"topic" + str(self.id_) + "\" -> \"" 
            scaled = " %#.4g" % (i[1]/self.nu_sum_)
            raw = " %#.4g" % i[1]
            s += "topic" + str(i[0]) + "\" [label = \"" + scaled + " (" + raw + ")\"];\n"
        if self.gamma_ < 0.1 and self.id_ < len(self.tau_):
            s = ""
        return s

def graphviz(topics, outfile_name="out.dot"):
    s = "digraph g {\n"
    s += "graph [\n"
    s += "rankdir = \"LR\"\n"

    s += "];\n"
    for i in topics:
        s += i.graphviz()
    s += "}"
    
    outfile = open(outfile_name, 'w')
    outfile.write(s)

def read_list(filename):
    return map(strip, open(filename).readlines())

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "USAGE: python visualize_topics.py vocab_file file_root num_topics num_words_per_record transition_multiplier_cutoff"
    else:
        vocab_file = sys.argv[1]
        file_root = sys.argv[2]
        num_topics = int(sys.argv[3])

        num_topic_words_printed = int(sys.argv[4])
        transition_multiplier_cutoff = float(sys.argv[5])
        topics = map(lambda x: Topic(x), range(num_topics + 1))
        vocab = read_list(vocab_file)

        repeat_delay = 0
        if len(sys.argv) == 7:
            repeat_delay = int(sys.argv[6])
        first_time = True

        while repeat_delay > 0 or first_time:
            first_time = False
            read_tau(file_root + ".tau", vocab, num_topic_words_printed, topics)
            read_gamma(file_root + ".gamma", topics)
            read_nu(file_root + ".nu", transition_multiplier_cutoff, topics)
            read_beta(file_root + ".beta", topics)
            graphviz(topics, file_root + ".dot")
            if repeat_delay > 0:
                try:
                    print open(file_root + ".lhood").readlines()[-1].strip()
                except IndexError:
                    break
            sleep(repeat_delay)
