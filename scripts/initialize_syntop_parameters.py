from math import log
from random import random

from topicmod.util import flags
from syntop_parameters_pb2 import *

flags.define_string("vocab", None, "Size of vocabulary")
flags.define_int("num_docs", None, "Numer of documents")
flags.define_int("num_topics", 128, "Number topics")
flags.define_string("model_name", "output/model", "Name of model")

flags.define_bool("finite", False, "Use finite model")
flags.define_bool("ignore_trans", False, "Use only documents")
flags.define_bool("ignore_docs", False, "Use only syntax")
flags.define_bool("shortcut_gsl", False, "Use closed form updates when possible")

flags.define_int("max_doc_iterations", 5, "Number of e-step rounds per-document")
flags.define_int("alpha_doc", 1.0, "DP parameter for doc distributions")
flags.define_int("alpha_trans",  1.0, "DP parameter for transition distributions")
flags.define_int("alpha_top", 1.0, "DP parameter for top-level stick-breaking distribution")
flags.define_int("vocab_sigma", 0.1, "Vocab hyperparametersx")

if __name__ == "__main__":
  flags.InitFlags()

  params = SyntopParameters()

  params.finite = flags.finite
  params.ignore_trans = flags.ignore_trans
  params.ignore_docs = flags.ignore_docs
  params.shortcut_gsl = flags.shortcut_gsl
  
  params.num_topics = flags.num_topics
  params.model_name = flags.model_name
  params.num_docs = flags.num_docs
  params.vocab_size = len(open(flags.vocab).readlines())
  print "Read %i terms from %s" % (params.vocab_size, flags.vocab)

  params.max_doc_iterations = flags.max_doc_iterations
  params.alpha_doc = flags.alpha_doc
  params.alpha_trans = flags.alpha_trans
  params.alpha_top = flags.alpha_top
  params.vocab_sigma = flags.vocab_sigma

  beta = [1.0] * flags.num_topics
  scale = params.alpha_top / (1.0 + params.alpha_top)
  beta_sum = 0.0

  for ii in xrange(flags.num_topics):
    beta[ii] = beta[ii - 1] * scale
    beta_sum += beta[ii]

  for ii in xrange(flags.num_topics):
    beta[ii] /= beta_sum

  beta_file = open("%s.beta" % params.model_name, 'w')
  for ii in beta:
    beta_file.write("%f\n" % ii)

  gamma_file = open("%s.gamma" % params.model_name, 'w')
  for ii in xrange(params.num_docs):
    for jj in xrange(params.num_topics):
      gamma_file.write("%f\n" % params.alpha_doc)

  nu_file = open("%s.nu" % params.model_name, 'w')
  for ii in xrange(params.num_topics):
    for jj in xrange(params.num_topics):
      if ii == jj: # Discourage self-transitions
        nu_file.write("%f\n" % (params.alpha_trans * 0.001))
      else:
        nu_file.write("%f\n" % params.alpha_trans)

  tau_file = open("%s.tau" % params.model_name, 'w')
  for ii in xrange(params.num_topics):
    tau = [10.0] * params.vocab_size
    for jj in xrange(params.vocab_size):
      tau[jj] += (random() - 0.5)
    tau_sum = log(sum(tau))
    for jj in tau:
      tau_file.write("%f\n" % (log(jj) - tau_sum))

  f = open("%s.params" % params.model_name, "wb")
  s = params.SerializeToString()
  f.write(s)
  f.close()
