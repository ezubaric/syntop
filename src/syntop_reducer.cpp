/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#include "syntop_reducer.h"
#include "gradient.h"

/*
 * End implemented functions for local mode.
 */
SyntopReducer::SyntopReducer(HadoopPipes::TaskContext& context) { // NOLINT
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

  vars_ = new VariationalParameters(*params_);

  check_order_ = false;
  lhood_ = false;

  index = -1;

  tau_coordinate_ = -1;
  nu_coordinate_ = -1;
}

void SyntopReducer::reduce(HadoopPipes::ReduceContext& context) { // NOLINT
  boost::split(temp_string_components_, context.getInputKey(),
               boost::is_any_of("_"));

  float sum = 0;
  while (context.nextValue()) {
    sum += HadoopUtils::toFloat(context.getInputValue());
  }



  if (boost::starts_with(temp_string_components_[0], "gamma")) {
    context.emit(context.getInputKey(), boost::lexical_cast<string>(sum));
  } else if (boost::starts_with(temp_string_components_[0], "lhood")) {
    // sum += GlobalLikelihoodTerm();
    sum += GlobalWeightTerm(vars_->nu_.get(),
                            vars_->beta_.get(),
                            params_->alpha_trans(),
                            params_->alpha_top(),
                            params_->finite());

    context.emit(context.getInputKey(), boost::lexical_cast<string>(sum));
  } else {
    if (boost::starts_with(temp_string_components_[2], "~")) {
      // cout << "optimizing" << endl;

      Optimize();
      Emit(&output);

      StringMap::const_iterator last = (output).end();
      for (StringMap::const_iterator itr = (output).begin();
           itr != last; itr++) {
        // cout << itr->first << "\t" << itr->second << endl;
        context.emit(itr->first, boost::lexical_cast<string>(itr->second));
      }

      output.clear();

      last = (output).end();
      for (StringMap::const_iterator itr = (output).begin();
           itr != last; itr++) {
        // cout << "output is\t" << itr->first << "\t" << itr->second << endl;
        // context.emit(itr->first, boost::lexical_cast<string>(itr->second));
      }

      index = boost::lexical_cast<int>(temp_string_components_[1]);
      vars_ = new VariationalParameters(*params_);

      display_matrix(vars_->tau_est_top_.get(), "tau_est_top is\n");
      display_vector(vars_->tau_est_bottom_.get(), "tau_est_bottom is\n");

      tau_coordinate_ = -1;
      nu_coordinate_ = -1;
    } else {
      ProcessKey(context.getInputKey(), sum);
      // cout << "processing\t" << context.getInputKey() << "\t" << sum << endl;
    }

    /*
    if (index == -1) {
      index = boost::lexical_cast<double>(temp_string_components_[1]);
      // reduceContext = context;
    } else {
      if (index != boost::lexical_cast<int>(temp_string_components_[1])) {
        Optimize();
        Emit(&output);

        StringMap::const_iterator last = (output).end();
        for (StringMap::const_iterator itr = (output).begin();
             itr != last; itr++) {
          // cout << itr->first << "\t" << itr->second << endl;
          context.emit(itr->first, boost::lexical_cast<string>(itr->second));
        }

        output.clear();

        index = boost::lexical_cast<int>(temp_string_components_[1]);
        vars_ = new VariationalParameters(*params_);

        tau_coordinate_ = -1;
        nu_coordinate_ = -1;
      }
    }
    */
  }
  // }
}

/*
void SyntopReducer::close() {
  if (tau_coordinate_ != -1 || nu_coordinate_ != -1) {
    StringMap::const_iterator last = (output).end();
    for (StringMap::const_iterator itr = (output).begin();
         itr != last; itr++) {
      // cout << itr->first << "\t" << itr->second << endl;
      reduceContext.emit(itr->first, boost::lexical_cast<string>(itr->second));
    }
  }
}
*/
