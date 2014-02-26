/*
 * strings.h
 *
 * Copyright 2007, Jonathan Chang
 *
 * Utils for manipulating strings
 */

#ifndef TOPICMOD_LIB_UTIL_STRINGS_H__
#define TOPICMOD_LIB_UTIL_STRINGS_H__
#include <assert.h>

#include <cstdlib>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

using std::ostringstream;
using std::string;
using std::vector;

#ifndef WIN32
  #define stricmp strcasecmp
  #define strnicmp strncasecmp
#endif

namespace lib_util {

// Returns the number of pieces found
int SplitStringUsing(const string& s, const char* delim,
                     vector<string>* result);
// convert double to string
string SimpleDtoa(double d);

string SimpleItoa(int n);

void StripWhitespace(string* s);

int ParseLeadingIntValue(const string& s);

int ParseLeadingIntValue(const char* s);

size_t ParseLeadingSizeValue(const char* s);

size_t ParseLeadingSizeValue(const string& s);

double ParseLeadingDoubleValue(const char* s);

double ParseLeadingDoubleValue(const string& s);

bool ParseLeadingBoolValue(const char* s);
/*
// This hash function for strings was taken from:
// http://www.partow.net/programming/hashfunctions/
namespace __gnu_cxx {
  template<> class hash<string> {
  public:
    size_t operator()(const string& str) const {
      size_t hash = 0;
      
      for(size_t ii = 0; ii < str.length(); ++ii) {
	hash = str[ii] + (hash << 6) + (hash << 16) - hash;
      }1
      return hash;
    } 
  };
}
*/
}

#endif  // TOPICMOD_LIB_UTIL_STRINGS_H__
