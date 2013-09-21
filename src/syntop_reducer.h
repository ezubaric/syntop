/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef SYNTOP_REDCUER_INCLUDED
#define SYNTOP_REDCUER_INLCUDED

#include "hadoop/Pipes.hh"
#include "hadoop/TemplateFactory.hh"
#include "hadoop/StringUtils.hh"

#include "syntop_standalone_reducer.h"

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */
class SyntopReducer : public SyntopStandaloneReducer, public
  HadoopPipes::Reducer {
  // protected:
  // HadoopPipes::ReduceContext *reduceContext;

 public:
  SyntopReducer(HadoopPipes::TaskContext& context); // NOLINT
  void reduce(HadoopPipes::ReduceContext& context); // NOLINT
  // void close(); // NOLINT
  ~SyntopReducer() {}
};

#endif
