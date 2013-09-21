/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#include "lens_doc.h"

using std::cerr;
using std::endl;

Document::Document(int id) {
  total_words_ = 0;
  number_sentences_ = 0;
  id_ = id;
}

int Document::get_id() {
  return id_;
}

void Document::FillWordArray() {
  // Count the number of words
  int number_words = 0;
  int word = 0;
  for (int i = 0; i < this->number_sentences_; i++) {
    // this->sentences_[i].head_->PrintTree(); cout << endl;
    // Note that this should be an assignment, not an increment.
    // FillWordArray increments number words starting from the position
    // designated by the third argument
    number_words = this->sentences_[i].head_->
      FillWordArray(NULL,
                    -1,
                    number_words,
                    this->sentences_[i].head_.get());
    // cout << "NW= " << number_words << ", ";
  }
  // cout << "NW=" << number_words << ", TW=" << total_words_ << endl;
  assert(number_words == total_words_);
  this->words_.reset(new LensNode*[total_words_]);
  for (int i = 0; i < this->number_sentences_; i++) {
    word = this->sentences_[i].head_->
      FillWordArray(words_.get(),
                    -1,
                    word,
                    this->sentences_[i].head_.get());
  }
}

void Document::SetNumberSentences(int number_sentences) {
  assert(number_sentences > 0);
  Sentence<LensNode>* sentences = new Sentence<LensNode>[number_sentences];
  this->sentences_.reset(sentences);
  this->number_sentences_ = number_sentences;
}

int Document::get_number_sentences_() {
  return this->number_sentences_;
}

int Document::get_number_words_() {
  return this->total_words_;
}

LensNode* Document::get_word_(int word) {
  return this->words_[word];
}

LensNode* Document::get_head_(int sentence) {
  return this->sentences_[sentence].head_.get();
}


/*
void Document::Visualize(string filename,
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

gsl_vector* Document::get_omega_() {
  return omega_.get();
}

double Document::get_omega_(int i) {
  return gsl_vector_get(omega_.get(), i);
}

/*
 * Note that this must be called after FillWordArray
 */
void Document::ResetPhi(int number_topics) {
  this->omega_.reset(gsl_vector_alloc(this->total_words_), gsl_vector_free);
  gsl_vector_set_all(omega_.get(), 0.0);
  if (this->total_words_ == 0) {
    cerr << "Is this being called without running FillWordArray first?" << endl;
  }
  for (int i = 0; i < this->total_words_; i++) {
    this->words_[i]->phi_.reset(gsl_vector_alloc(number_topics),
                                gsl_vector_free);
    gsl_vector_set_all(this->words_[i]->phi_.get(), 0.0);
    uniform(this->words_[i]->phi_.get());
  }
}

void Document::ClearPhi() {
  this->omega_.reset();
  for (int i = 0; i < this->total_words_; i++) {
    this->words_[i]->phi_.reset();
  }
}

Document* ParseDocument(string line, int max_sent) {
  vector<string> sents;
  // Subtract two because the first fields are the number of sentences
  // and the id of the document
  int num_sents = SplitStringUsing(line, "\t", &sents) - 2;
  // This assertion tests that assumption
  int doc_id = ParseLeadingIntValue(sents[0]);
  assert(num_sents == ParseLeadingIntValue(sents[1]));

  if (max_sent > 0 && num_sents > max_sent) num_sents = max_sent;

  Document* new_doc = new Document(doc_id);
  new_doc->SetNumberSentences(num_sents);
  for (int ii = 0; ii < num_sents; ++ii) {
    // Add one because the first field is the number of sentences
    new_doc->AddSentence(ii, sents[ii + 2]);
  }
  new_doc->FillWordArray();
  return new_doc;
}

void Document::AddSentence(int sentence_num, string line) {
  total_words_ += sentences_[sentence_num].LoadSentence(line);
}
