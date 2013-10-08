import sys
from string import strip
from random import sample, shuffle
from math import log, exp
from sets import dict_sample

MAX_TOKENS = 1000

INCLUSION_THRESHOLD = 0.00025
# How many words a POS needs before it will be split up by document
SINGLE_SEED_WORD_COUNT = 100
# How many words a POS needs before it will be culled
MIN_WORDS_IN_POS_CATEGORY = 10
# How many words a topic needs to be included
MIN_WORDS_IN_TOPIC = 10
MAX_TRIES = 10

def load_vocab(filename, lower=False):
    words = map(strip, open(filename, 'r').readlines())
    count = 0
    d = {}

    for i in words:
        if lower:
            word = i.lower()
        else:
            word = i
        if i in d:
            print "DUPLICATE:", i

        word = word.split()[1]
        d[word] = count
        d[count] = word
        count += 1
    print "vocab size:", len(words), len(d), len(filter(lambda x: type(x)==type(""), d.keys()))
    return d

def SampleWeightedDictionary(d, num):
    temp_list = []
    for i in d:
        temp_list += [i] * int(500 * d[i])
    return sample(temp_list, num)

def TopicLineFromCounts(voc, counts, smoothing=10.0):
    print "TOPIC LINE", counts
    s = ""
    normalizer = len(voc) / 2 * smoothing + sum(counts.values())
    for i in xrange(len(voc) / 2):
        if i != 0:
            s += "\n"
        val = log(counts.get(i, 0.0) + float(smoothing)) - log(normalizer)
        if counts.has_key(i):
            print voc[i], counts.get(i, 0.0)
        s += str(val)
    size = len(counts.keys())
    return s, size

def TopicLine(doc, voc, pos, smoothing=10):
    new_doc = []
    for i in doc:
        if pos == i[1]:
            print "+", voc[i[0]], i[1],
            new_doc += [i[0]]
    doc = new_doc
    print ""

    #print "MAX:", max(map(lambda x: tf[x]*idf[x], new_doc)), "MIN:", min(map(lambda x: tf[x]*idf[x], new_doc))
    #doc = filter(lambda x: tf[x]*idf[x] > INCLUSION_THRESHOLD, new_doc)

    #print "Words retained:", map(lambda x: voc[x], filter(lambda y: tf[y]*idf[i] < INCLUSION_THRESHOLD, new_doc))
    #print doc

    # We choose the word that has the best tf-idf
    #source_word = doc[0]
    #for i in doc:
        #print source_word, i, voc[i], tf[i], idf[i], tf[i]*idf[i]
    #    if tf[i]*idf[i] > tf[source_word]*idf[source_word]:
    #        source_word = i

    #related_words = linsim[voc[source_word]]

    #doc = filter(lambda x: voc[x] in related_words or x==source_word, doc)
    #print voc[source_word], map(lambda x: voc[x], doc)

    #print doc
    counts = {}
#    source_count = max(linsim[voc[source_word]].values() + [smoothing])
    for i in doc:
        counts[i] = counts.get(voc[i], smoothing) + 1#linsim[voc[source_word]].get(voc[i], source_count)

    #print counts
    return TopicLineFromCounts(voc, counts)

def GenerateTopics(number_topics, voc, pos_counts, docs, output_file):

    #assert(number_topics <= len(docs))
    print number_topics, len(docs), 0 <= number_topics <= len(docs), 0<=number_topics, number_topics <= len(docs)
    shuffled_keys = docs.keys()
    shuffle(shuffled_keys)
    #print len(sampled_docs)
    print "Creating ", len(docs) , " topics of size ", len(voc) / 2

    good_lines = []
    bad_lines = []
    # "First try to get good topics that don't have pont masses"
    #for i in bnc.function_pos:
    #    line, length = TopicLineFromCounts(voc, bnc.get_pos(i))
    #    good_lines.append(line)

    # Convert the pos counts into entropies
    entropies = {}
    single_seeds = []
    for i in pos_counts:
        print "Considering pos", i
        # We're actually not using entropies any more, just raw word counts
        count = len(pos_counts[i])
        if count > SINGLE_SEED_WORD_COUNT:
            entropies[i] = count
        else:
            single_seeds.append(i)

    print "ENTROPIES", entropies

    # Do all the simple parts of speech that don't have that many words
    for pos in single_seeds:
        print "##", pos
        line, length = TopicLineFromCounts(voc, pos_counts[pos])
        good_lines.append(line)

    if len(entropies) != 0:
        # Choose the part of speech of the focus word
        pos = dict_sample(entropies)
        tries = 0

        print shuffled_keys[:10]
        for i in shuffled_keys:
            print "This pos has ", len(pos_counts[pos]), " representatives"
            print "selected doc", i
            tries += 1
            if tries <= MAX_TRIES:
                line, length = TopicLine(docs[i], voc, pos)
            else:
                line, length = TopicLineFromCounts(voc, pos_counts[pos])

            if len(good_lines) == number_topics:
                break
            if length >= MIN_WORDS_IN_TOPIC:
                good_lines.append(line)

                # Choose the part of speech of the focus word
                # for the next document
                pos = dict_sample(entropies)
                tries = 0
            else:
                bad_lines.append(line)

    print "Using ", len(good_lines), " good lines."
    o = open(output_file, 'w')
    o.write("\n".join(good_lines))
    if len(good_lines) < number_topics:
        bad_lines[:(number_topics-len(good_lines))]
        o.write("\n")
        o.write("\n".join(bad_lines))

def LoadDocuments(vocab, input="sentences.dat"):
    corpus = {}
    pos_count = {}
    doc = 0
    cur_doc = []
    token_count = 0

    import glob
    glob.glob(input);

    for file in glob.glob(input):
        data = open(file, 'r');
        for i in data:
            if doc % 50000 == 0 and doc > 0:
            #return docs, pos_count
                print doc, len
            print "test data", i
            fields = i.split("\t")
            if len(fields)==1 and len(cur_doc) > 0:
                print doc, token_count, len(cur_doc)
                corpus[doc] = cur_doc
                doc += 1
                token_count = 0
                cur_doc = []
            elif len(fields)==1:
                continue
            else:
                tokens = i.split()
                tokens.reverse()
                while (len(tokens) > 0):
                    cur_token = int(tokens.pop())
                    if cur_token >= 0:
                        token_count += 1
                        pos = int(tokens.pop())
                        if token_count < MAX_TOKENS:
                            cur_doc += [[cur_token, pos]]
                        if not pos_count.has_key(pos):
                            pos_count[pos] = {}
                        pos_count[pos][cur_token] = pos_count[pos].get(cur_token, 0) + 1
                        tokens.pop()
    return corpus, pos_count

def LoadDocumentWords(vocab, input="sentences.dat"):
    doc = 0
    cur_doc = {}
    pos_count = {}
    docs = []
    roles = []
    for i in open(input, 'r'):
        doc += 1
        if doc % 50000 == 0:
            #return docs, pos_count
            print doc
        fields = i.split()
        if len(fields)==1 and len(cur_doc) > 0:
            docs.append(cur_doc)
            cur_doc = {}
        elif len(fields)==1:
            continue
        else:
            tokens = i.split()
            tokens.reverse()
            while (len(tokens) > 0):
                cur_token = int(tokens.pop())
                if cur_token >= 0:
                    cur_doc[cur_token] = cur_doc.get(cur_token, 0) + 1
                    pos = tokens.pop()
                    if not pos_count.has_key(pos):
                        pos_count[pos] = {}
                    pos_count[pos][cur_token] = pos_count[pos].get(cur_token, 0) + 1
                    pos_count.get(pos, {})
                    tokens.pop()
    if len(cur_doc) > 0:
        docs.append(cur_doc)
    return docs, pos_count

def old_main():
#if __name__=="__main__":
    if len(sys.argv) != 5:
        print "USAGE: python seed.py dat_file vocab num_topics output_file"
    else:
        voc = load_vocab(sys.argv[2])
        #linsim = LinSim(voc, ["/n/fs/topics/data/linsim/all.lsp", "/Users/jbg/topics/data/all.lsp"])

        docs, pos = LoadDocuments(voc, sys.argv[1])
        print "Read in", len(docs), "documents with vocab", len(voc)

        # Remove nearly empty parts of speech
        pos_keys = pos.keys()
        for i in pos_keys:
            print i
            if len(pos[i]) < MIN_WORDS_IN_POS_CATEGORY or sum(pos[i].values()) < 100:
                print "Deleting pos ", i
                del pos[i]
        print ""

        GenerateTopics(int(sys.argv[3]), voc, pos, docs, sys.argv[4])

if __name__=="__main__":
#def new_main():
    # added by Ke Zhai
    import sys

    import random;
    from numpy import *;

    num_topic = int(sys.argv[1]);
    num_doc = int(sys.argv[2]);
    num_voc = int(sys.argv[3]);
    output_prefix = sys.argv[4];

    beta_param = 1.0;

    o = open(output_prefix+"gamma.origin", 'w')
    for i in xrange(num_doc):
        o.write("\t".join([str(1.0) for j in xrange(num_topic)]));
        o.write("\n");

        #a = array([random.random() for j in xrange(num_topic)]);
        #o.write("\t".join(map(str, a)) + "\n")
    print "output seed gamma..."

    o = open(output_prefix+"nu.origin", 'w')
    for i in xrange(num_topic+1):
        o.writelines("\t".join([str(1.0) for j in xrange(num_topic)]));
        o.writelines("\n");
    print "output seed nu..."

    o = open(output_prefix+"tau.origin", 'w')
    for i in xrange(num_topic):
        a = array([random.random() for j in xrange(num_voc)]);
        a = log(a/sum(a));
        o.write("\t".join(map(str, a)) + "\n")
        # o.write("\t".join([str(log(1.0/num_voc)) for j in xrange(num_voc)]) + "\n");
    print "output seed tau..."

    o = open(output_prefix+"beta.origin", 'w')
    uniform = True;
    if not uniform:
        beta_prime = random.beta(1.0, beta_param, num_topic);
        one_minus_beta_prime = 1-beta_prime;
        beta = beta_prime;
        for i in xrange(1, num_topic-1):
            beta[i] = beta[i] * one_minus_beta_prime[0:i].prod()
        beta[num_topic-1] = 1.0-beta[0:num_topic-1].sum()
    else:
        beta = ones(num_topic);
        beta = beta/sum(beta);
    print "beta is ", beta
    o.write("\t".join(map(str, beta)) + "\n")
    print "output seed beta..."
