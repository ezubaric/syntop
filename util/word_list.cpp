/*
 * word_list.cpp
 *
 * Copyright 2009, Jordan Boyd-Graber
 *
 * Wrapper class for an array that allows it to be used as a flag type in flags.h
 */

#include <string>

#include "topicmod/lib/util/word_list.h"
#include "topicmod/lib/util/strings.h"

namespace lib_util {

WordList::WordList(const string& s) {
  SplitStringUsing(s, ",", &words_);
}

size_t WordList::size() const {
  return words_.size();
}

std::ostream & operator<<(std::ostream & out, WordList const & list ) {
  string s = "[";
  for (int ii = 0; ii < (int)list.size(); ++ii) {
    s += list[ii];
    if (ii < (int)list.size() - 1) s += ", ";
  }
  s += "]";
  out << s;
  return out;
}

const string& WordList::operator[] (int idx) const {
  return words_[idx];
}
}
