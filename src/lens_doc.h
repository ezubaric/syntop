/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef LENS_DOC_INCLUDED
#define LENS_DOC_INCLUDED

#include <gsl/gsl_vector.h>

#include <boost/scoped_array.hpp>
#include <iomanip>
#include <fstream>

#include "topicmod/lib/util/strings.h"

#include "util.h"
#include "sentence.h"

using lib_util::SplitStringUsing;
using lib_util::ParseLeadingIntValue;
using boost::scoped_array;
using boost::shared_ptr;

class Document {
 public:
  LensNode* get_word_(int word);
  LensNode* get_head_(int sentence);
  explicit Document(int id);
  int get_number_sentences_();
  int get_number_words_();
  void AddSentence(int sentence_num, string line);
  void SetNumberSentences(int number_sentences);
  void FillWordArray();
  int get_id();

  friend class LensModel;
  double get_omega_(int i);
  double get_xi_(int i);
  gsl_vector* get_omega_();
  void ResetPhi(int number_topics);
  /*
  void Visualize(string filename,
                 gsl_vector* gamma,
                 const vector<string>& vocab,
                 const vector<string>& roles,
                 const gsl_matrix* tau);
  */
 protected:
  int number_sentences_;
  int total_words_;
  scoped_array<LensNode*> words_;

  void ResetNuExpInteraction();
  void AddNuExpInteraction(int index, double contribution);
  void ClearPhi();
  int index_;
  shared_ptr<gsl_vector> omega_;
 private:
  Document(const Document& copy);
  scoped_array< Sentence<LensNode> > sentences_;
  int id_;
};

/*
void LensDocument::Visualize(string filename,
                             gsl_vector* gamma,
                             const vector<string>& vocab,
                             const vector<string>& roles,
                             const gsl_matrix* tau) {
  fstream outfile;
  double gamma_sum = gsl_blas_dsum(gamma);
  outfile.open(filename.c_str(), ios::out);

  outfile << "digraph g{" << endl;
  outfile << "graph [" << endl;
  outfile << "label = \"(";
  for(unsigned int i=0; i<gamma->size; i++) {
    if(i!=0) outfile << ", ";
    outfile << i << ":" << gsl_vector_get(gamma, i) / gamma_sum;
  }
  outfile << ")\"" << endl;
  outfile << "rankdir = \"LR\"" << endl;
  outfile << "];" << endl;

  for(int n=0; n<total_words_; n++) {
    LensNode* word = words_[n];

    outfile << "\"" << n << "\" [" << endl;
    outfile << "\tlabel = \"" << vocab[word->term_];

    gsl_vector* phi = word->phi_.get();
    outfile << "| omega = " << fixed << setprecision(3) << get_omega_(n);
    outfile << "| xi = " << fixed << setprecision(3) << get_xi_(n);
    for(unsigned int i=0; i<phi->size; i++) {
      outfile << "| " << i << ":" << fixed << setprecision(3)
              << gsl_vector_get(phi, i) << "("
              << gsl_matrix_get(tau, i, word->term_) << ")";
    }
    // Print out topic proportions

    outfile << "\"" << endl;
    outfile << "\tshape = record" << endl;
    outfile << "];" << endl;

    for(int c=0; c<word->number_children_; c++) {
      int child_index = word->children_positions_[c];
      outfile << "\"" << n << "\" -> \"" << child_index
              << "\" [label = \"" << roles[words_[child_index]->relation_]
              << "\"]" << endl;
    }
  }
  outfile.close();
}
*/
Document* ParseDocument(string line, int max_sent);

#endif
