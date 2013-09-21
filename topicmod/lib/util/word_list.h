/*
 * word_list.h
 *
 * Copyright 2009, Jordan Boyd-Graber
 *
 * Wrapper class for a list option in flags.h
 */
#ifndef TOPICMOD_LIB_UTIL_WORD_LIST_H_
#define TOPICMOD_LIB_UTIL_WORD_LIST_H_

#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;

namespace lib_util {

class WordList {
 public:
  explicit WordList(const string& s);
  friend std::ostream & operator<<(std::ostream &, WordList const &);
  size_t size() const;

  const string& operator[] (int idx) const;
 private:
    vector<string> words_;
};
}

#endif  // TOPICMOD_LIB_UTIL_WORD_LIST_H_
