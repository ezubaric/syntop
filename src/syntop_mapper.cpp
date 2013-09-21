/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#include "syntop_mapper.h"
// This is included here to make life easier
#include "gradient.h"

#include <iostream>
#include <fstream>

using HadoopUtils::toString;

/*
 * corpus is a tree structured corpus with labeled roles
 *
 * model_name is a string representing this model, which is used to
 * ouput checkpoints and other information
 *
 * num_topics is the number of the topics used in the finite model or
 * the truncation level in the infinite case
 *
 * alpha_doc is the scalar controlling the Dirichlet (process /
 * distribution) that \theta comes from
 *
 * alpha_trans is the same, but for the transition probabilities
 *
 * alpha_top is the parameter for the GEM distribution used to
 * generate the top level stick-breaking weights, but is ignored in
 * the finite model
 *
 * finite is true when a finite model is desired
 */
SyntopMapper::SyntopMapper(HadoopPipes::TaskContext& context) { // NOLINT
  // TODO(zhaike): add in initialization, i.e., load parameters from distributed
  // cache.

  SyntopParameters* params = new SyntopParameters();

  params->set_vocab_size(context.getJobConf()->getInt("syntop.vocab.number"));
  params->set_num_docs(context.getJobConf()->getInt("syntop.doc.number"));
  params->set_num_topics(context.getJobConf()->getInt("syntop.topic.number"));
  params->set_model_name(context.getJobConf()->get("syntop.model.name"));

  params->set_finite(context.getJobConf()->getBoolean("syntop.model.finite"));
  params->set_ignore_trans(context.getJobConf()->
                           getBoolean("syntop.model.ignore.trans"));
  params->set_ignore_docs(context.getJobConf()->
                          getBoolean("syntop.model.ignore.docs"));
  params->set_shortcut_gsl(context.getJobConf()->
                           getBoolean("syntop.model.shortcut.gsl"));

  params->set_alpha_doc(context.getJobConf()->getFloat("syntop.alpha.doc"));
  params->set_alpha_trans(context.getJobConf()->getFloat("syntop.alpha.trans"));
  params->set_alpha_top(context.getJobConf()->getFloat("syntop.alpha.top"));

  params_.reset(params);
}

void SyntopMapper::map(HadoopPipes::MapContext& context) { // NOLINT
  /*
  cout << "parameter is\t"
       << params_->vocab_size() << "\t"
       << params_->num_docs() << "\t"
       << params_->num_topics() << "\t"
       << params_->model_name() << "\t"
       << params_->finite() << "\t"
       << params_->ignore_trans() << "\t"
       << params_->ignore_docs() << "\t"
       << params_->shortcut_gsl() << "\t"
       << params_->alpha_doc() << "\t"
       << params_->alpha_trans() << "\t"
       << params_->alpha_top() << "\t" << endl;
  */

  StringMap count;

  vars_.reset(new VariationalParameters(*params_));
  doc_mapper = new DocumentMapper(vars_.get(), params_.get());

  line = context.getInputValue();
  doc = ParseDocument(line, -1);

  // shall we emit the doc_id together with the key?
  int doc_id = doc->get_id();

  doc->ResetPhi(params_->num_topics());
  doc_mapper->EStep(difference_tolerance, doc, doc_id);
  doc_mapper->Emit(doc, &count);

  for (int ii = 0; ii < params_->num_topics(); ++ii) {
    // Tau
    string key = "tau_" + SimpleItoa(ii) + "_-1";
    context.emit(key, boost::lexical_cast<string>(count[key]));

    for (int vv = 0; vv < params_->vocab_size(); ++vv) {
      key = "tau_" + SimpleItoa(ii) + "_" + SimpleItoa(vv);
      context.emit(key, boost::lexical_cast<string>(count[key]));
    }

    string flag = "tau_" + SimpleItoa(ii) + "_~";
    context.emit(flag, boost::lexical_cast<string>(0));

    // Gamma
    key = "gamma_" + SimpleItoa(doc_id) + "_" + SimpleItoa(ii);
    context.emit(key, boost::lexical_cast<string>(count[key]));

    // Nu ss
    for (int jj = 0; jj < params_->num_topics() + 1; ++jj) {
      // It's jj before ii because I wanted to put it inside the outer loop, but
      // the first index of nu ss needs to go from 0 to K instead of 0 to K-1
      string key_doc = "nu_" + SimpleItoa(jj) + "_" + SimpleItoa(ii) + "_doc";
      context.emit(key_doc, boost::lexical_cast<string>(count[key_doc]));

      string key_trans = "nu_" + SimpleItoa(jj) + "_" +
        SimpleItoa(ii) + "_trans";
      context.emit(key_trans, boost::lexical_cast<string>(count[key_trans]));

      string key_flag = "nu_" + SimpleItoa(jj) + "_~";
      context.emit(key_flag, boost::lexical_cast<string>(count[key_flag]));
    }

    // context.emit("nu_" + SimpleItoa(params_->num_topics() + 1) + "_"
    // + SimpleItoa(params_->num_topics()),
    // boost::lexical_cast<string>(0));
  }

  // context.emit("tau_" + SimpleItoa(params_->num_topics() + 1) + "_"
  // + SimpleItoa(params_->vocab_size() + 1),
  // boost::lexical_cast<string>(0));

  // likelihood
  string key = "lhood_1_1";
  context.emit(key, boost::lexical_cast<string>(count[key]));

  /*
  // Gamma
  key = "gamma_" + SimpleItoa(doc_id) + "_~";
  context.emit(key, HadoopUtils::toString(0));

  // start emitting indicator
  for (int ii = 0; ii < params_->num_topics(); ++ii) {
    // Tau
    string key = "tau_" + SimpleItoa(ii) + "_~";
    context.emit(key, HadoopUtils::toString(0));
  }

  // Nu ss
  for (int jj = 0; jj < params_->num_topics() + 1; ++jj) {
    // It's jj before ii because I wanted to put it inside the outer loop, but
    // the first index of nu ss needs to go from 0 to K instead of 0 to K-1
    string key_doc = "nu_" + SimpleItoa(jj) + "_~_doc";
    context.emit(key_doc, HadoopUtils::toString(0));

    string key_trans = "nu_" + SimpleItoa(jj) + "_~_trans";
    context.emit(key_trans, HadoopUtils::toString(0));
  }

  key = "lhood_1_~";
  context.emit(key, HadoopUtils::toString(0));
  // finish emitting indicators
  */
}
