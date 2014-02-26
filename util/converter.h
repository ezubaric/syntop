/*
 * converter.h
 *
 * Copyright 2007, Jonathan Chang
 */

#ifndef TOPICMOD_LIB_UTIL_CONVERTER_H_
#define TOPICMOD_LIB_UTIL_CONVERTER_H_

#include <string>
#include <vector>

#include "topicmod/lib/util/strings.h"
#include "topicmod/lib/util/word_list.h"

using std::string;

namespace lib_util {

template<typename T>
struct converter {
  T operator()(const char* arg);
};

template<>
struct converter<string> {
  string operator()(const char* arg) {
    return string(arg);
  }
};

template<>
struct converter< vector<string> > {
  vector<string> operator()(const char* arg) {
    vector<string> result;
    SplitStringUsing(arg, ",", &result);
    return result;
  }
};

template<>
struct converter<int> {
  int operator()(const char* arg) {
    return ParseLeadingIntValue(arg);
  }
};

template<>
struct converter<size_t> {
  int operator()(const char* arg) {
    return ParseLeadingSizeValue(arg);
  }
};

template<>
struct converter<lib_util::WordList> {
  lib_util::WordList operator()(const char* arg) {
    return WordList(arg);
  }
};

template<>
struct converter<double> {
  double operator()(const char* arg) {
    return ParseLeadingDoubleValue(arg);
  }
};


template<>
struct converter<bool> {
  bool operator()(const char* arg) {
    return ParseLeadingBoolValue(arg);
  }
};
}

#endif  // TOPICMOD_LIB_UTIL_CONVERTER_H_
