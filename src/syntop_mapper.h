/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef SYNTOP_MAPPER_INCLUDED
#define SYNTOP_MAPPER_INCLUDED

#include "syntop_standalone_mapper.h"

#include "hadoop/Pipes.hh"
#include "hadoop/TemplateFactory.hh"
#include "hadoop/StringUtils.hh"

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */

class SyntopMapper : public SyntopStandaloneMapper, public HadoopPipes::Mapper {
 public:
  SyntopMapper(HadoopPipes::TaskContext& context); // NOLINT
  void map(HadoopPipes::MapContext& context); // NOLINT
};

#endif
