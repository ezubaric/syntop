/*
 * Copyright 2008, Jonathan Chang
 *
 */

#include <string>
#include <vector>

#include "topicmod/lib/util/strings.h"

namespace lib_util {

const char* kWhitespace = "\t\n\v\f\r ";

string SimpleDtoa(double d) {
  ostringstream s;
  s << d;
  return s.str();
}

string SimpleItoa(int n) {
  ostringstream outs;
  outs << n;
  return outs.str();
}


int SplitStringUsing(const string& s, const char* delim,
                     vector<string>* result) {
  string::size_type pos = 0;
  string::size_type index = 0;
  result->clear();
  while (pos < s.size()) {
    // Position pos to the first non-delimiter character.
    pos = s.find_first_not_of(delim, pos);
    if (pos == string::npos) {
      break;
    }
    // Find the first delimiter.
    index = s.find_first_of(delim, pos);
    if (index == string::npos) {
      index = s.size();
    }
    result->push_back(s.substr(pos, index - pos));

    // Position pos to the next possible non-delimiter.
    pos = index + 1;
  }
  return result->size();
}


double ParseLeadingDoubleValue(const string& s) {
  return strtod(s.c_str(), NULL);
}

double ParseLeadingDoubleValue(const char* s) {
  return strtod(s, NULL);
}

bool ParseLeadingBoolValue(const char* s) {
  if (!stricmp(s, "true")) {
    return true;
  } else if (!stricmp(s, "false")) {
    return false;
  } else {
    assert(false);
    return false;
  }
}

size_t ParseLeadingSizeValue(const char* s) {
  return strtoull(s, NULL, 10);
}

size_t ParseLeadingSizeValue(const string& s) {
  return strtoull(s.c_str(), NULL, 10);
}

int ParseLeadingIntValue(const char* s) {
  return strtol(s, NULL, 10);
}

int ParseLeadingIntValue(const string& s) {
  return strtol(s.c_str(), NULL, 10);
}

void StripWhitespace(string* s) {
  string::size_type index = s->find_first_not_of(kWhitespace);
  if (index == string::npos) {
    *s = "";
    return;
  }
  *s = s->substr(index);
  index = s->find_last_not_of(kWhitespace);
  *s = s->substr(0, index + 1);
}
}
