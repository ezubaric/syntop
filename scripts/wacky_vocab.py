
from glob import glob
from topicmod.lib.util import flags
import gzip
from collections import defaultdict
import operator

if __name__ == "__main__":
    flags.define_string("output", None, "Where we write vocab file")
    flags.define_bool("use_lemma", False, "Use the lemma for words")
    flags.define_int("vocab_size", 10000, "Vocab size")
    flags.define_string("pattern", None, "pattern for input files")
    flags.define_bool("use_lower", True, "Transform to lower case")
    flags.InitFlags()

    count = defaultdict(int)
    rel = set()
    pos = set()

    # Get the count of all words
    for ii in glob(flags.pattern):
        for jj in gzip.open(ii, 'rb'):
            if '\t' in jj:
                fields = jj.split('\t')
                pos.add(fields[2])
                rel.add(fields[-1].strip())
                if flags.use_lemma:
                    word = fields[1]
                else:
                    word = fields[0]

                if flags.use_lower:
                    count[word.lower()] += 1
                else:
                    count[word] += 1

    print("%i tokens counted with max count of %i" % (len(count), max(count.values())))

    # Write the vocab
    num_written = 0
    with open("%s.word" % flags.output, 'w') as outfile:
        for ww, cc in sorted(count.iteritems(), key=operator.itemgetter(1), reverse=True):
            num_written += 1

            outfile.write("%s\n" % ww)

            if num_written > flags.vocab_size:
                break

    with open("%s.rel" % flags.output, 'w') as outfile:
        for ii in rel:
            outfile.write("%s\n" % ii)

    with open("%s.pos" % flags.output, 'w') as outfile:
        for ii in pos:
            outfile.write("%s\n" % ii)
