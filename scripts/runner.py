import sys
import os

from topicmod.util import flags
from syntop_parameters_pb2 import *

flags.define_int("num_iterations", 1, "Number of iterations")
flags.define_string("model_name", "output/model", "Where we find data")

flags.define_string("corpus", None, "The source corpus")

flags.define_bool("hadoop", False, "Do we use hadoop or local batch")
flags.define_bool("doc_step", True, "Do we call the document-centric parts")
flags.define_bool("merge_step", True, "Do we merge doc step results (and compute new topics)")
flags.define_bool("update_global", True, "Do we compute new transition and DP variational parameters")

class Array:
  def __init__(self, name):
    self._rows = {}
    self._name = name

  def __getitem__(self, index):
    if not index in self._rows:
      self._rows[index] = defaultdict(float)
    return self._rows[index]

  def __iter__(self):
    for ii in self._rows:
      yield self._rows[ii]

  def parse(self, key, val):
    assert key.startswith(self._name), "%s/%s label mismatch" % (key, val)
    name, ii, jj = key.split("_")
    self._rows[int(ii)][int(jj)] += val

  def dims(self):
    dim1 = max(self._rows) + 1
    dim2 = max(len(self._rows[x]) for x in self._rows) + 1
    return dim1, dim2

  def add_prior_and_normalize(self, prior):
    d1, d2 = self.dims()
    for ii in xrange(d1):
      for jj in xrange(d2):
        self[ii][jj] = (self[ii][jj] + prior) / self[ii][-1]

  def dump(self, root):
    o = open("%s.%s" % (root, self._name), 'w')

    d1, d2 = self.dims()
    for ii in xrange(d1):
      for jj in xrange(d2):
        o.write("%f\n" % self[ii][jj])

def parse_estep(params):
  """
  Returns the likelihood (of the document part)
  """ 

  array_lookup = {}

  array_lookup["nutrans"] = Array("nu_exp_topic_counts")
  array_lookup["nudoc"] = Array("nu_exp_doc_counts")
  array_lookup["tau"] = Array("tau")
  array_lookup["lhood"] = Array("lhood")

  for ii in open("%s.dump" % flags.model_name):
    key, val = ii.split()
    bucket = key.split("_")[0]

    array_lookup[bucket].parse(key, val)

  array_lookup["tau"].normalize(params.vocab_sigma)

  for ii in array_lookup:
    array_lookup[ii].dump(params.model_name)

if __name__ == "__main__":
  flags.InitFlags()
  
  # Read parameters
  params = SyntopParameters()
  f = open("%s.params" % flags.model_name, "rb")
  params.ParseFromString(f.read())

  for ii in xrange(flags.num_iterations):
    if flags.doc_step:
      if flags.hadoop:
        # Copy the distributed cache
        None
      else:
        for ii in os.popen("./batch_syntop %s %s.params" % (flags.corpus, flags.model_name)):
          print ii

    if flags.merge_step:
      lhood = parse_estep(params)
      print("Iteration %i Doc lhood %f" % (ii, lhood))

    # Update nu and beta
    if flags.update_global:
      None
