
from glob import glob
from topicmod.lib.util import flags
from topicmod.lib.util.parse_reader import Tree
from collections import defaultdict
import gzip

class Vocab:
    def __init__(self, input_base):
        self._lookup = {}
        self._relations = {}
        self._pos = {}

        with open("%s.word" % input_base) as infile:
            for ww in infile:
                self._lookup[ww.strip()] = len(self._lookup)

        with open("%s.pos" % input_base) as infile:
            for pos in infile:
                self._lookup["OOV-%s" % pos.lower().strip()] = len(self._lookup)
                self._pos[pos.strip().lower()] = len(self._pos)

        with open("%s.rel" % input_base) as infile:
            for pos in infile:
                self._relations[pos.strip().lower()] = len(self._relations)

    def word(self, word, pos):
        if word in self._lookup:
            return self._lookup[word]
        else:
            return self._lookup["OOV-%s" % pos]

    def relation(self, relation):
        return self._relations[relation.lower()]

    def pos(self, pos):
        return self._pos[pos.lower()]

def tree_from_wacky_sentence(sentence, vocab, use_lemma=False, use_lower=True):
    children = defaultdict(list)
    nodes = {}

    for parent, child in [(x.parent, x.offset) for x in sentence]:
        children[parent].append(child)

    head = None
    for ww in sentence:
        if use_lemma and use_lower:
            word = ww.lemma.lower()
        elif use_lemma:
            word = ww.lemma
        elif use_lower:
            word = ww.word.lower()
        else:
            word = ww.word
        word = vocab.word(word, ww.pos)
        relation = vocab.relation(ww.relation.strip())
        pos = vocab.pos(ww.pos)

        nodes[ww.offset] = (ww.offset, word, relation, pos)

        if x.parent == 0:
            head = x.offset

    if head:
        return Tree(head, children, nodes)
    else:
        return None

class WackyLine:

    def __init__(self, line):
        self.word, self.lemma, self.pos, self.offset, \
          self.parent, self.relation = line.lower().split("\t")
        self.offset = int(self.offset)
        self.parent = int(self.parent)

class WackyDocument:

    def __init__(self, sentences, title):
        self._raw = None
        self.author = ""
        self.title = title
        self._sentences = sentences
        self._id = None

    def tokens(self):
        num_sentences = max(self._sentences) + 1
        for ii in xrange(num_sentences):
            for jj in self._sentences[ii]:
                yield jj.word

    def pos_tags(self):
        num_sentences = max(self._sentences) + 1
        for ii in xrange(num_sentences):
            for jj in self._sentences[ii]:
                yield jj.pos

    def lemmas(self, stemmer):
        num_sentences = max(self._sentences) + 1
        for ii in xrange(num_sentences):
            for jj in self._sentences[ii]:
                yield jj.lemma

    def chunked_nps(self, min_count=2):
        for ii in self._sentences:
            tokens = []
            count = 0
            in_np = False

            for jj in self._sentences[ii]:
                if jj.pos.startswith("np"):
                    if not in_np:
                        tokens.append("<")
                    in_np = True
                elif in_np:
                    tokens.append(">")
                    in_np = False
                    count += 1
                if not ">" in jj.word and not "<" in jj.word:
                    tokens.append(jj.word)

            if count >= min_count:
                yield " ".join(tokens)
            else:
                yield ""

    def num_sentences(self):
        return len(self._sentences)

    def relations(self):
        num_sentences = max(self._sentences) + 1
        for ii in xrange(num_sentences):
            for jj in self._sentences[ii]:
                yield jj.rel

def documents(filename):
    content_file = gzip.open(filename, 'r')

    buffer = defaultdict(list)
    current_sentence = 0

    for ii in latin_iter(content_file):
        # print "wacky-sentences?", ii.encode("ascii", "ignore")
        if ii.startswith("<s>"):
            continue
        elif ii.startswith(u'<text id'):
            title = ii.split("wikipedia:")[1].replace('">', '')
            print "title is", title
        elif ii.startswith("</s>"):
            current_sentence += 1
        elif ii.startswith("</text>"):
            wd = WackyDocument(buffer, title)
            print len(buffer), wd.num_sentences()
            if wd.num_sentences() > 20:
                yield wd

            buffer = defaultdict(list)
            title = ""
        else:
            buffer[current_sentence].append(WackyLine(ii))


def latin_iter(input_file):
    line = input_file.readline()
    while line:
        yield line.decode("latin-1")
        line = input_file.readline()
    return

if __name__ == "__main__":
    flags.define_string("pattern", None, "pattern for input files")
    flags.define_string("vocab", None, "vocab base filename (.word, .pos, .rel)")
    flags.define_string("doc_out", None, "where we write document")
    flags.define_string("title_out", None, "where we write the title out")
    flags.define_bool("use_lemma", False, "Use lemma form")
    flags.define_bool("use_lower", True, "Use lower case")
    flags.define_bool("use_relation", True, "Use relationship in final printout (alternative is POS)")
    flags.InitFlags()

    vocab = Vocab(flags.vocab)

    doc_out = open("%s.dat" % flags.doc_out, 'w')
    title_out = open("%s.title" % flags.doc_out, 'w')
    for ii in glob(flags.pattern):
        for doc in documents(ii):
            title_out.write("%s\n" % doc.title.strip())
            sentences = []
            for ii in doc._sentences:
                tree = tree_from_wacky_sentence(doc._sentences[ii], vocab)
                if tree:
                    tree_string = tree.PrintString(flags.use_relation)
                else:
                    tree_string = ""
                if tree_string:
                    sentences.append(tree_string)
            doc_out.write("%i %s\n" % (len(sentences), " ".join(sentences)))
