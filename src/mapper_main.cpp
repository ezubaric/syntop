/*
 * Copyright 2009 Jordan Boyd-Graber
 */
#include "variational_parameters.h"
#include "document_mapper.h"

#include "topicmod/lib/util/flags.h"

DEFINE_string(model_name, "", "Where we load data from");
DEFINE_string(doc, NULL, "Filename of all documents");
DEFINE_int(vocab_size, -1, "Number of terms");
DEFINE_int(max_line, 1048576, "Max line length");
DEFINE_int(max_sent, 200, "Maximum number of sentences");

int main(int argc, char *argv[]) {
  InitFlags(argc, argv);

  SyntopParameters params = SyntopParameters();
  params.set_model_name(FLAGS_model_name);
  params.set_vocab_size(FLAGS_vocab_size);

  int num_docs = 0;
  ifstream docfile(FLAGS_doc.c_str());
  while (docfile.good()) {
    string line;
    getline(docfile, line);
    num_docs++;
  }
  params.set_num_docs(num_docs);

  scoped_ptr<VariationalParameters>
    var(new VariationalParameters(params, true));
  docfile.seekg(0);
  scoped_ptr<DocumentMapper> mapper(new DocumentMapper(var.get(), &params));
  scoped_ptr<Document> doc;
  StringMap count;

  while (docfile.good()) {
    string line;
    getline(docfile, line);
    doc.reset(ParseDocument(line, FLAGS_max_sent));

    mapper->Emit(doc.get(), &count);
  }
}
