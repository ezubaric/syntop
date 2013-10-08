import sys
import random

def sample_wr(population, k):
    """
    Chooses k random elements (with replacement) from a population
    """

    n = len(population)
    _random, _int = random.random, int  # speed hack 
    result = [None] * k
    for i in xrange(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return result

class PartOfSpeech:
    """
    Defines a part of speech for generating synthetic data
    """

    def __init__(self, name, elements, sucessors, min_dep, max_dep, ID):
        self.name = name
        self.ID = ID
        self.elements = elements
        self.sucessors = sucessors

        self.max_dep = max_dep
        self.min_dep = min_dep

    def doc_vocab(self, vocab):
        """ Choose a topic and apply the vocab to it """
        return map(lambda x: vocab[x], random.sample(self.elements, 1)[0])

    def choose_children(self):
        num_children = random.sample(range(self.min_dep, self.max_dep + 1), 1)[0]
        return sample_wr(self.sucessors, num_children)

def walk_tree(pos, words, current):
    """ Returns the string representation of a node and its children """

    # Pick the word
    s = str(random.sample(words[current], 1)[0]) + " "

    # Get the role
    s += str(pos[current].ID) + " "

    # Get the children
    children = pos[current].choose_children()
    s += str(len(children)) + " "

    for i in children:
        s += walk_tree(pos, words, i)

    # Walk back to the parent
    s += "-1 "

    return s

def create_document(parts_of_speech, vocab, max_sen=10,
                    start="verb", min_sen=1):
    """ 
    Returns a string corresponding to the randomly created document
    """

    word_classes = {}
    for i in parts_of_speech:
        word_classes[i] = parts_of_speech[i].doc_vocab(vocab)

    print "WORD_CLASSES: ", word_classes

    num_sentences = random.sample(range(min_sen, max_sen), 1)[0]

    s = str(num_sentences) + "\t"
    for i in xrange(num_sentences):
        s += walk_tree(parts_of_speech, word_classes, start).strip() + "\t"
        
    s += "\n"

    return s


def create_parts_of_speech():
    d = {}
    d["verb"] = PartOfSpeech("verb", [["sits", "runs", "falls", "walks"], ["falls", "climbs", "surges", "bucks"], ["ponders", "queries", "discusses"], ["despairs", "dreads", "hates", "mourns", "fears"]], ["noun", "prep"], 1, 3, 0)
    d["noun"] = PartOfSpeech("noun", [["COW", "SHEEP", "PONY"], ["GRAD_STUDENT", "PROFESSOR", "PHD_CANDIDATE"], ["STOCK", "MUTUAL_FUND", "SHARE"]], ["adj", "prep", "det"], 0, 2, 1)
    d["prep"] = PartOfSpeech("prep", [["on", "about", "over", "with"]], ["noun"], 1, 1, 2)
    d["det"] = PartOfSpeech("det", [["this", "that", "the", "a"]], [], 0, 0, 3)
    d["adj"] = PartOfSpeech("adj", [["stupid", "insolent"], ["American", "Russian", "German"], ["happy", "elated", "glad"], ["red", "blue", "white", "purple"], ["evil"]], [], 0, 0, 4)
    return d

def write_roles(parts_of_speech, output_file):
    l = [""]*len(parts_of_speech)
    for i in parts_of_speech:
        l[parts_of_speech[i].ID] = i

    o = open(output_file, 'w')
    for i in l:
        o.write(parts_of_speech[i].name + "\n")
    o.close()
        
def create_vocab(words, output_file):
    o = open(output_file, 'w')

    # Get a list without duplicates
    d = []
    for i in words:
        for j in words[i].elements:
            for k in j:
                if not k in d:
                    d.append(k)

    print "VOCAB: ",
    # Get a dictionary so we can just write the 
    id_map = {}
    for i in xrange(len(d)):
        print d[i], i,
        id_map[d[i]] = i
    print "\n"

    print 'ID_MAP: ', id_map

    for i in d:
        o.write(i + "\n")
    o.close()
    
    return id_map
        
            
if __name__=="__main__":
    output_root = sys.argv[1]
    num_docs = int(sys.argv[2])
    num_sen = int(sys.argv[3])

    pos = create_parts_of_speech()
    vocab = create_vocab(pos, output_root + ".voc")
    write_roles(pos, output_root + ".roles")

    o = open(output_root + ".0.prs", 'w')
    o_doc = open(output_root + ".doc", 'w')
    #o.write(str(num_docs) + "\n")
    for ii in xrange(num_docs):
        o.write("%i\t%s" % (ii, create_document(pos, vocab, num_sen)))
        #o.write("%i\t%s" % (ii, create_document(pos, vocab, num_sen, 'noun', 1)))
        o_doc.write("%i\n" % ii)
