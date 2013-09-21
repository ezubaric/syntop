/*
 *  sentence.h
 *
 *
 *  Created by Jordan Boyd-Graber on 27.09.07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SENTENCE_INCLUDED
#define SENTENCE_INCLUDED

// The number of numbers needed to adequately define a node; it's the
// number of children, the word, and the relation
#define NODE_DEF_NUMBERS_NEEDED 3

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/scoped_ptr.hpp>

#include "node.h"

using boost::is_any_of;
using boost::split;

// class Model;
// struct SamplingNode: public Node;
// {
//      void Sample(Model* model, bool init);
//      std::string StateString();
// }

template<class T> class Sentence {
 public:
        Sentence();
        int LoadSentence(string input);
        boost::scoped_ptr<T> head_;
 private:
        Sentence(const Sentence& copy);
};

template<class T>
Sentence<T>::Sentence() {}

/*
template<class T> std::string SampleSentence::StateString() {
        return this->head_->StateString();
}*/

/*
template<class T> void SampleSentence::Sample(Model* model, bool init) {
        // We use the const to designate that it doesn't have a parent synset.
        this->head_->Sample(model, NO_PARENT, init);
}
*/


/*
 * This turns a text representation of a parse into a tree in memory.
 */
template<class T>
int Sentence<T>::LoadSentence(string line) {
  vector<T*> parents;
  vector<string> tokens;

  // TODO(jbg): This creates problems when there are two spaces between tokens,
  // and doesn't work like python's split
  split(tokens, line, is_any_of(" \t"));
  int number_words = 0;
  int token_position = 0;
  int last_position = -3;
  int num_tokens=(int)tokens.size();
  T* cur = NULL;
  while (token_position < num_tokens) {
    // This should hopefully ensure that if we add more to the
    // definition of a node, it will be incorporated into this function to
    // read them in.  We must have either read NODE_DEF_NUMBERS_NEEDED *or*
    // we read in a -1, which tells us to backtrack up one level.
    assert(last_position + NODE_DEF_NUMBERS_NEEDED == token_position ||
           last_position + 1 == token_position);
    last_position = token_position;

    int token = atoi(tokens[token_position].c_str());
    // cout << token_position << ":" << token << endl;
    ++token_position;
    if (token >= 0) {
      // It's a word

      int relation = atoi(tokens[token_position++].c_str());
      // cout << "Relation: " << relation << endl;
      int number_children = atoi(tokens[token_position++].c_str());
      if (parents.size()) {
        ++number_words;
        cur = parents.back()->AddChild(token, relation, number_children);
      } else {  // It's the root node
        // The next node to be read should be the first node after the root node
        number_words = 1;
        assert(token_position == NODE_DEF_NUMBERS_NEEDED);
        cur = new T();
        this->head_.reset(cur);

        head_->relation_ = relation;
        head_->term_ = token;
        head_->SetNumberChildren(number_children);
      }
      parents.push_back(cur);
    } else {
      // It's telling us to go back up
      cur = parents.back();
      parents.pop_back();
    }
  }
  assert(parents.empty() && cur == head_.get());
  return number_words;
}

#endif
