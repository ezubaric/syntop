/*
 * Flag library
 *
 * Copyright 2007, Jonathan Chang
 *
 * flags.h
 */

#ifndef __UTIL_FLAGS_H__
#define __UTIL_FLAGS_H__

#include <vector>
#include <string>

#ifdef WIN32                             // NOLINT
#include "getopt.h"                      // NOLINT
#else                                    // NOLINT
#include <getopt.h>                      // NOLINT
#endif                                   // NOLINT

#include <iostream>                      // NOLINT

#include "topicmod/lib/util/converter.h" // NOLINT

using std::cout;
using std::endl;
using std::vector;
using lib_util::WordList;

#define DEFINE_string(name, default, description) string FLAGS_##name; Flag<string> FLAGSCONTAINER_##name(#name, default, description, &FLAGS_##name) // NOLINT
#define DEFINE_int(name, default, description) int FLAGS_##name; Flag<int> FLAGSCONTAINER_##name(#name, default, description, &FLAGS_##name) // NOLINT
#define DEFINE_list(name, default, description) WordList FLAGS_##name(default); Flag<WordList> FLAGSCONTAINER_##name(#name, FLAGS_##name, description, &FLAGS_##name) // NOLINT
#define DEFINE_size(name, default, description) size_t FLAGS_##name; Flag<size_t> FLAGSCONTAINER_##name(#name, default, description, &FLAGS_##name) // NOLINT
#define DEFINE_double(name, default, description) double FLAGS_##name;  Flag<double> FLAGSCONTAINER_##name(#name, default, description, &FLAGS_##name) // NOLINT
#define DEFINE_bool(name, default, description) bool FLAGS_##name; Flag<bool> FLAGSCONTAINER_##name(#name, default, description, &FLAGS_##name) // NOLINT

class FlagBase;

class FlagSingleton {
 public:
  static void AddFlag(FlagBase* flag);
  static void AddUnparsedArgument(const char* s);

  static void Build();
  static void Print();
  static void Parse(int index, const char* arg);

  static const vector<string>& unparsed_arguments();
  static struct option* longopts();
 private:
  static FlagSingleton& get();

  struct option* options_;
  vector<string> unparsed_arguments_;
  vector<FlagBase*> flags_;
};

class FlagBase {
 public:
  FlagBase(const char* name, const char* description);

  virtual void Print(std::ostream& s) = 0; // NOLINT
  virtual void Parse(const char* arg) = 0;

  virtual ~FlagBase();

  const char* name();
  const char* description();
 private:
  const char* name_;
  const char* description_;
};

template<typename T>
class Flag : public FlagBase {
 public:
  Flag(const char* name, const T& defaultval, const char* description, T* val);
  virtual void Parse(const char* arg);
  virtual void Print(std::ostream& s); // NOLINT

 private:
  T* val_;
  lib_util::converter<T> converter_;
};

template<typename T>
void Flag<T>::Parse(const char* arg) {
  *val_ = converter_(arg);
}

template<typename T>
void Flag<T>::Print(std::ostream& s) { // NOLINT
  const T& val = *val_;
  s << val;
}

template<typename T>
Flag<T>::Flag(const char* name,
           const T& defaultval,
           const char* description,
           T* val) : FlagBase(name, description), val_(val) {
  *val_ = defaultval;
}

void InitFlags(int argc, char *argv[]);
const vector<string>& UnparsedFlagArguments();

#endif  // __UTIL_FLAGS_H__
