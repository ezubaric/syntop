/*
 * flags.cpp
 * Copyright 2007, Jonathan Chang
 * Edited by Jordan Boyd-Graber
 *
 * Adds a flag parser
 */

#include <string>
#include <vector>

#include "topicmod/lib/util/flags.h"

void FlagSingleton::AddFlag(FlagBase* flag) {
    get().flags_.push_back(flag);
}

struct option* FlagSingleton::longopts() {
    return get().options_;
}

void FlagSingleton::AddUnparsedArgument(const char* s) {
  get().unparsed_arguments_.push_back(string(s));
}

const vector<string>& FlagSingleton::unparsed_arguments() {
  return get().unparsed_arguments_;
}

FlagSingleton& FlagSingleton::get() {
  static FlagSingleton* flag_singleton = new FlagSingleton;
  return *flag_singleton;
}

FlagBase::FlagBase(const char* name,
                   const char* description) : name_(name),
                                              description_(description) {
  FlagSingleton::AddFlag(this);
}

const char* FlagBase::name() {
  return name_;
}


const char* FlagBase::description() {
  return description_;
}


FlagBase::~FlagBase() { }


void FlagSingleton::Build() {
  get().options_ = new struct option[get().flags_.size() + 2];

  get().options_[0].name = "help";
  get().options_[0].has_arg = 0;
  get().options_[0].flag = NULL;
  get().options_[0].val = 0;

  for (int ii = 1; ii < (int)get().flags_.size() + 1; ++ii) {
    get().options_[ii].name = (char*)get().flags_[ii - 1]->name();
    get().options_[ii].has_arg = 1;
    get().options_[ii].flag = NULL;
    get().options_[ii].val = 0;
  }

  get().options_[get().flags_.size() + 1].name = NULL;
  get().options_[get().flags_.size() + 1].has_arg = 0;
  get().options_[get().flags_.size() + 1].flag = NULL;
  get().options_[get().flags_.size() + 1].val = 0;
}

void FlagSingleton::Print() {
  FlagSingleton fs = get();
  cout << "Command line options: " << endl;
  cout << "------------------------------------------" <<
    "--------------------------------------" << endl;

  int flag_width = 15;
  int value_width = 25;
  int description_width = 40;
  cout.width(flag_width); cout << "Flag name"; //NOLINT
  cout.width(value_width); cout << "Value"; //NOLINT
  cout.width(description_width); cout << "Description" << endl; // NOLINT
  cout << "------------------------------------------" <<
    "--------------------------------------" << endl;

  for (int ii = 0; ii < (int)fs.flags_.size(); ++ii) {
    cout.width(flag_width);
    cout << fs.flags_[ii]->name();
    cout.width(value_width);
    fs.flags_[ii]->Print(cout);
    cout.width(description_width);
    string description = fs.flags_[ii]->description();

    int last_printed_char = 0;
    while (last_printed_char < (int)description.size()) {
      if (last_printed_char > 0) {
        cout << endl;
        cout.width(flag_width + value_width + description_width);
      } else {
        cout.width(description_width);
      }

      // Just print the rest if we can
      int next = -1;
      if (last_printed_char + description_width > (int)description.size()) {
        next = description.size();
      } else {
        // We break at a space if we can
        next = description.rfind(' ', last_printed_char + description_width);
        // Otherwise, let the chips fall where they may
        if (next <= last_printed_char)
          next = last_printed_char + description_width;
      }
      cout << description.substr(last_printed_char, next);
      last_printed_char = next;
    }
    cout << endl;
  }
  cout << "------------------------------------------" <<
    "--------------------------------------" << endl;
}

void FlagSingleton::Parse(int index, const char* arg) {
  get().flags_[index]->Parse(arg);
}


void PrintFlags() {
  FlagSingleton::Print();
}

const vector<string>& UnparsedFlagArguments() {
  return FlagSingleton::unparsed_arguments();
}

void InitFlags(int argc, char *argv[]) {
  FlagSingleton::Build();

  for (int ii = 1; ii < argc; ++ii) {
    if (argv[ii][0] != '-')
      FlagSingleton::AddUnparsedArgument(argv[ii]);
  }

  int long_index;
  while (1) {
    int c = getopt_long(argc,
                        argv,
                        "",
                        FlagSingleton::longopts(),
                        &long_index);
    if (c == -1) {
      break;
    }

    // There's a better way to handle help requests, but this will
    // have to do for now
    if (long_index == 0) FlagSingleton::Print();
    assert(long_index > 0);
    // We subtract 1 because "--help" is in the zeroth position in
    // options, but doesn't actually appear in the index
    FlagSingleton::Parse(long_index - 1, optarg);
  }
}
