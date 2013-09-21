/*
 * Copyright 2010 Jordan Boyd-Graber
 */

#include "syntop_mapper.h"
#include "syntop_reducer.h"
#include "syntop_partitioner.cpp"

#include  <stdint.h>  // <--- to prevent uint64_t errors!

#include <algorithm>
#include <limits>
#include <fstream>
#include <string>

#include "hadoop/Pipes.hh"
#include "hadoop/TemplateFactory.hh"
#include "hadoop/StringUtils.hh"

using std::string;

int main(int argc, char *argv[]) {
  return HadoopPipes::runTask(HadoopPipes::TemplateFactory<
                              SyntopMapper,
                              SyntopReducer,
                              SyntopPartitioner>());

  /*
  {
    cout << "update beta phase started..." << endl;
    
    shared_ptr<const SyntopParameters> params_ = param_container;
    scoped_ptr<VariationalParameters> vars_;
    vars_.reset(new VariationalParameters(*params_));

    if (!params_->finite()) {
      display_vector(vars_->beta_.get(), "updated beta is\n");
      optimize_beta(vars_->beta_.get(), vars_->gamma_.get(),
                    vars_->nu_.get(), params_->ignore_docs(),
                    params_->ignore_trans(), params_->alpha_top(),
                    params_->alpha_doc(), params_->alpha_trans(),
                    params_->model_name());
      display_vector(vars_->beta_.get(), "updated beta is\n");

      ofstream yourfile;
      yourfile.open((FLAGS_directory + "_params/current.beta").c_str());
      // yourfile.open((FLAGS_directory + "_params/tmp.beta").c_str());
      for (int i = 0; i < params_->num_topics(); i++) {
	yourfile << gsl_vector_get(vars_->beta_.get(), i) << "\t";
      }
      yourfile.close();
    }

    cout << "update beta finished..." << endl;
  }
  */
}
