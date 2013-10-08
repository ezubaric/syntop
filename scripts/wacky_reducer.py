#
# wacky_reducer.py
#
# File to turn protocol buffers into a test-only input file readable by
# mapreduce implementation of syntactic topic model.

from collections import defaultdict

from nltk import FreqDist

from topicmod.corpora.proto.corpus_pb2 import *
from topicmod.corpora.corpus_vocab_wrapper import CorpusVocabWrapper
from topicmod.util import flags
from parse_reader import *

flags.define_int("docs_per_file", 100, "Number of documents per file")
flags.define_int("vocab_size", 5000, "Maximum vocabulary size")
flags.define_bool("remove_stop", False, "remove stopwords")
flags.define_bool("use_lemma", False, "Use lemmas instead of raw tokens")
flags.define_bool("use_relation", False,
                  "Use relation (synset) instead of pos")
flags.define_glob("vocab_source", None, "Where we get initial vocabulary")
flags.define_string("output_path", None,
                    "Where we write the translated corpuss")
flags.define_string("output_filename", "wacky_en_reduced.index",
                    "Filename of index")
flags.define_int("min_length", 100, "Number of characters in document line")

class CorpusTranslator:

  def __init__(self, output_path, use_lemma, docs_per_file):
    self.output_base = output_path
    self.document_list = []
    self.use_lemma = use_lemma

    # A lookup for each language
    self.vocab = defaultdict(dict)
    self.roles = defaultdict(dict)
    self.output_corpus = Corpus()

    self.docs_per_file = docs_per_file
    self.num_files = 0

  def merge_and_truncate_vocab(self, list_of_original_vocabs, vocab_size,
                               remove_stop):
#    print "original_vocabs is", list_of_original_vocabs

    source_vocab = {}
    freq_total = defaultdict(FreqDist)
    for ii in list_of_original_vocabs:
      source_vocab[ii] = CorpusVocabWrapper(ii)
      print "source_vocab is", source_vocab[ii]

    print "USE RELATION", flags.use_relation

    if flags.use_relation:
      roles = source_vocab[ii].synsets
    else:
      roles = source_vocab[ii].pos

    for ii in roles:
      for jj in roles[ii]:
        entry = roles[ii][jj]
        if entry.original in self.roles[ii]:
          assert self.roles[ii][entry.original] == entry.id
        self.roles[ii][entry.original] = entry.id

    print "USE LEMMA", self.use_lemma
    for ii in source_vocab:
      source = source_vocab[ii].tokens
      #print "source is", source
      if self.use_lemma:
        source = source_vocab[ii].lemmas

      for lang in source:
        for kk in source[lang]:
          entry = source[lang][kk]
#          print "entry is", entry
          if entry.stop_word:
            print "entry is stop_word", entry.stop_word, entry.original
          if remove_stop and entry.stop_word:
            continue
          freq_total[lang].inc(entry.original, entry.frequency)

        print "freq_total is", lang, len(freq_total[lang])

    # Now that we have sum of all frequencies, remove low frequency words
    print "Total words observed:", [len(freq_total[x]) for x in freq_total]

    self.vocab = defaultdict(dict)
    for lang in freq_total:
      num_types = 0
      for jj in freq_total[lang]:
        print "word is", lang, jj.encode("ascii", "ignore"), num_types
#        if num_types < 100 or u'nor' in jj:
#          print lang, jj.encode("ascii", "ignore"), num_types
        self.vocab[lang][jj] = num_types
        num_types += 1

        if num_types >= flags.vocab_size:
          print "warning: exceed vocab_size", num_types, flags.vocab_size, self.vocab[lang]
          break

  def translate_document(self, doc_num, doc_filename):

    if self.num_files % self.docs_per_file == 0:
      split = self.num_files // self.docs_per_file
      self.parse_output = open("%s.%i.prs" % (self.output_base, split), 'w')
      self.doc_output = open("%s.%i.doc" % (self.output_base, split), 'w')
    self.num_files += 1

    doc = Document()
    print "Reading document from", doc_filename
    doc.ParseFromString(open(doc_filename, 'r').read())

    t_doc = Document()
    t_doc.id = doc.id

    sentences = []
    for sent in doc.sentences:
      print "sentence starts..."

      head = -1
      valid = True
      nodes = {}
      children = defaultdict(set)

      for word in sent.words:
        token = self.vocab_lookup(doc.language, word.token, word.lemma,
                                  word.pos)
#        print token, self.string_lookup(doc.language, word.token, word.lemma)

#        print "word and token is", word, word.token, token
#        print "Token %s %s" % (token, self.string_lookup(doc.language,
#                                                         word.token,
#                                                         word.lemma))
#        print self.string_lookup(doc.language, word.token, word.lemma)

        nodes[word.offset] = (word.offset, token, word.synset, word.pos)
        #print self.string_lookup(doc.language, word.token, word.lemma), word.pos, word.offset
        
        if word.parent == 0:
          if head == -1:
            head = word.offset
          else:
            valid = False
        else:
          children[word.parent].add(word.offset)

      if head != -1 and valid:
        t = Tree(head, children, nodes)
        sent_string = t.PrintString(flags.use_relation)
#        print "sent_string is |%s|" % sent_string
        if sent_string:
          sentences.append(sent_string.strip())

      print "sentence finishes..."

    line = "%i\t%i\t%s" % (doc_num, len(sentences), "\t".join(sentences))
    line = line.strip()
    #print "|%s|" % line
    if len(line) > flags.min_length:
      self.parse_output.write("%s\n" % line)
      self.doc_output.write("%s\n" % doc.title)
      return True
    else:
      return False

  def string_lookup(self, language, original, lemma):
    term = original
    if self.use_lemma:
      term = lemma
    term = self.term_mapping[language][term].original

    return term

  def vocab_lookup(self, language, original, lemma, pos):
    """
    token = self.vocab_lookup(doc.language, word.token, word.lemma, word.pos)

    Convert a string to an integer.

    As a side-effect, expand vocab to include possible OOV terms (that is, the
    terms that might be added to account for words not explicitly in the
    vocabulary).  These terms are of the form OOV-$pos, where $pos is a part fo
    speech.
    """
#    print "language is", original, lemma, pos

    term = original
    if self.use_lemma:
      term = lemma
    term = self.term_mapping[language][term].original
#    print term.encode("ascii", "ignore"), pos

    if not term in self.vocab[language]:
      print "MISS", term.encode("ascii", "ignore"), "UNK-%i" % pos
      term = "UNK-%i" % pos
      #print term,
      if not term in self.vocab[language]:
        self.vocab[language][term] = max(self.vocab[language].values()) + 1

#    print "self.vocab[language][term] is", original, self.vocab[language][term]
    return self.vocab[language][term]

  def set_vocab(self, corpus_part):
    corpus_part_vocab = CorpusVocabWrapper(corpus_part)
    self.term_mapping = corpus_part_vocab.tokens
    
    for kk in corpus_part_vocab.tokens.keys():
      print "corpus token is"
      for ll in corpus_part_vocab.tokens[kk].keys():
        print "token is", ll, corpus_part_vocab.tokens[kk][ll];

    if self.use_lemma:
      self.term_mapping = corpus_part_vocab.lemmas

  def write_roles(self, filename):
    """
      Write the protocol buffer that includes the list of all documents and the
      complete vocabulary.  Because vocab can expand as documents are
      translated, this should be done after all documents are translated.
    """

    o = codecs.open("%s.roles" % flags.output_path, 'w', 'latin-1')

    for ll in self.roles:
      for ii in self.roles[ll]:
        o.write("%i\t%s\n" % (ll, ii.strip()))

    o.close()

  def write_vocab(self, filename):
    """
      Write the protocol buffer that includes the list of all documents and the
      complete vocabulary.  Because vocab can expand as documents are
      translated, this should be done after all documents are translated.
    """

    o = codecs.open("%s.voc" % flags.output_path, 'w', 'latin-1')

#    print "self.vocab is", self.vocab
#    print "self.vocab is", self.vocab[0].keys()

    for ll in self.vocab:
      invert_vocab = {};
      for ii in self.vocab[ll]:
        #invert_vocab[len(invert_vocab)] = ii;
        invert_vocab[self.vocab[ll][ii]] = ii;
      for ii in xrange(len(invert_vocab)):
        o.write("%i\t%i\t%s\n" % (ll, ii, invert_vocab[ii]));
      #for ii in self.vocab[ll]:
      #  o.write("%i\t%s\t%i\n" % (ll, ii, self.vocab[ll][ii]))

    o.close()

if __name__ == "__main__":
  flags.InitFlags()

  ct = CorpusTranslator(flags.output_path, flags.use_lemma,
                        flags.docs_per_file)

  # Build vocab
#  print "flags.vocab_size is", flags.vocab_size
  ct.merge_and_truncate_vocab(flags.vocab_source, flags.vocab_size,
                              flags.remove_stop)

  # Rewrite corpus
  doc_num = 0
  for ii in flags.vocab_source:
    cp = Corpus()
    cp.ParseFromString(open(ii, 'r').read())
    print "Reading documents from ", ii
    root = ii.rsplit("/", 1)[0]

    ct.set_vocab(ii)

    for jj in cp.doc_filenames:
      new_doc = ct.translate_document(doc_num, "%s/%s" % (root, jj))
      if new_doc:
        doc_num += 1
      print "Done with doc", jj

  ct.write_vocab("%s/%s" % (flags.output_path, flags.output_filename))
  ct.write_roles("%s/%s" % (flags.output_path, flags.output_filename))
