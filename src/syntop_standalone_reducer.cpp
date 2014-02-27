/*
 * Copyright 2009 Jordan Boyd-Graber
 */

// #include "syntop_reducer.h"
#include "syntop_standalone_reducer.h"
#include "gradient.h"

using std::cerr;

/*
 * Begin implemented functions for local mode.
 */
SyntopStandaloneReducer::SyntopStandaloneReducer
(shared_ptr<const SyntopParameters> params) {
  params_ = params;

  vars_ = new VariationalParameters(*params_);

  check_order_ = false;
  lhood_ = false;

  tau_coordinate_ = -1;
  nu_coordinate_ = -1;
}

StringMap SyntopStandaloneReducer::reduce(string type, int index,
                                          StringMap* buffer) {
  StringMap output;

  StringMap::const_iterator last = (*buffer).end();
  for (StringMap::const_iterator itr = (*buffer).begin(); itr != last; itr++) {
    boost::split(temp_string_components_, itr->first, boost::is_any_of("_"));
    if (boost::starts_with(temp_string_components_[0], type) &&
        boost::lexical_cast<double>(temp_string_components_[1]) == index) {
      ProcessKey(itr->first, itr->second);
    }
  }

  if (lhood_) {
    double adjustment = GlobalWeightTerm(vars_->nu_.get(),
                                         vars_->beta_.get(),
                                         params_->alpha_trans(),
                                         params_->alpha_top(),
                                         params_->finite());

    vct_accum(0, adjustment, vars_->lhood_.get());
  }

  double old_likelihood_test = 0;
  if (nu_coordinate_ >= 0) {
    {
      nu_params p;

      p.transition_counts = vars_->nu_exp_topic_counts_.get();
      p.nu_gamma_interaction = vars_->nu_exp_doc_counts_.get();

      gsl_vector* beta = vars_->beta_.get();
      shared_ptr<gsl_vector> alpha;
      alpha.reset(gsl_vector_alloc(beta->size), gsl_vector_free);
      p.alpha = alpha.get();
      gsl_vector_memcpy(p.alpha, beta);
      gsl_vector_scale(p.alpha, params_->alpha_trans());

      gsl_vector_view nu_view = gsl_matrix_row(vars_->nu_.get(),
                                               nu_coordinate_);
      gsl_vector* nu = &nu_view.vector;

      p.topic = nu_coordinate_;

      old_likelihood_test = nu_lhood(nu, &p);
    }
  }

  Optimize();

  double new_likelihood_test = 0;
  if (nu_coordinate_ >= 0) {
    {
      nu_params p;

      shared_ptr<gsl_vector> nu_exp(gsl_vector_alloc(params_->num_topics()),
                                    gsl_vector_free);
      p.nu_exp = nu_exp.get();

      p.transition_counts = vars_->nu_exp_topic_counts_.get();
      p.nu_gamma_interaction = vars_->nu_exp_doc_counts_.get();

      gsl_vector* beta =  vars_->beta_.get();
      shared_ptr<gsl_vector> alpha;
      alpha.reset(gsl_vector_alloc(beta->size), gsl_vector_free);
      p.alpha = alpha.get();
      gsl_vector_memcpy(p.alpha, beta);
      gsl_vector_scale(p.alpha, params_->alpha_trans());

      gsl_vector_view nu_view = gsl_matrix_row(vars_->nu_.get(),
                                               nu_coordinate_);
      gsl_vector* nu = &nu_view.vector;

      p.topic = nu_coordinate_;

      new_likelihood_test = nu_lhood(nu, &p);
    }

    if (-old_likelihood_test + new_likelihood_test > 0.000001) {
      cerr << "something wrong with nu update: " << old_likelihood_test << "\t"
           << new_likelihood_test << " " << new_likelihood_test -
        old_likelihood_test << endl;
    }

    cout << "likelihood difference is (negative is increase) "
         << new_likelihood_test-old_likelihood_test << endl;
  }

  Emit(&output);
  return output;
}

/*
 * End implemented functions for local mode.
 */

void SyntopStandaloneReducer::ProcessKey(string key, double val) {
  boost::split(temp_string_components_, key, boost::is_any_of("_"));

  if (temp_string_components_[0] == "gamma") {
    return;
  }

  int ii = boost::lexical_cast<double>(temp_string_components_[1].c_str());
  int jj = boost::lexical_cast<double>(temp_string_components_[2].c_str());

  if (temp_string_components_[0] == "tau") {
    if (tau_coordinate_ == -1) tau_coordinate_ = ii;
    if (check_order_) {
      assert(tau_coordinate_ == ii);
      assert(nu_coordinate_ == -1);
      assert(!lhood_);
    }

    if (jj != -1) {
      mtx_accum(ii, jj, val, vars_->tau_est_top_.get());
    } else {
      vct_accum(ii, val, vars_->tau_est_bottom_.get());
    }
  }

  if (boost::starts_with(temp_string_components_[0], "nu")) {
    if (nu_coordinate_ == -1) nu_coordinate_ = ii;
    if (check_order_) {
      assert(nu_coordinate_ == ii);
      assert(tau_coordinate_ == -1);
      assert(!lhood_);
    }

    if (temp_string_components_[3] == "trans") {
      // cout << "accumated nu trans is " << ii << "\t" << jj << endl;
      // display_matrix(vars_->nu_exp_topic_counts_.get(), "nu_exp_topic_counts_
      // is ");
      mtx_accum(ii, jj, val, vars_->nu_exp_topic_counts_.get());
      // cout << "accumated nu trans is " << ii << "\t" << jj << endl;
      // display_matrix(vars_->nu_exp_topic_counts_.get(), "nu_exp_topic_counts_
      // is ");
    } else {
      assert(temp_string_components_[3] == "doc");
      mtx_accum(ii, jj, val, vars_->nu_exp_doc_counts_.get());
    }
  }

  if (temp_string_components_[0] == "lhood") {
    // KE ZHAI: comment this out
    lhood_ = true;
    if (check_order_) {
      assert(nu_coordinate_ == -1);
      assert(tau_coordinate_ == -1);
    }
    vct_accum(0, val, vars_->lhood_.get());
    // cout << "likelihood is " << val << endl;
  }
}

void SyntopStandaloneReducer::Emit(StringMap* emission_counts) {
  if (nu_coordinate_ >= 0) {
    shared_ptr<gsl_vector> nu_vector;
    nu_vector.reset(gsl_vector_alloc(params_->num_topics()));
    gsl_vector_set_all(nu_vector.get(), 0);

    nu_params p;

    double nu_lhood_ = 0;
    shared_ptr<gsl_vector> nu_exp(gsl_vector_alloc(params_->num_topics()),
                                  gsl_vector_free);
    p.nu_exp = nu_exp.get();

    p.transition_counts = vars_->nu_exp_topic_counts_.get();
    p.nu_gamma_interaction = vars_->nu_exp_doc_counts_.get();

    gsl_vector* beta =  vars_->beta_.get();
    shared_ptr<gsl_vector> alpha;
    alpha.reset(gsl_vector_alloc(beta->size), gsl_vector_free);
    p.alpha = alpha.get();

    gsl_vector_memcpy(p.alpha, beta);

    gsl_vector_scale(p.alpha, params_->alpha_trans());

    // This is less than or equal because there is a special start state that
    // cannot be transitioned into for each sentence head
    for (unsigned int topic = 0; topic <= params_->num_topics(); topic++) {
      // This scales self-transitions downward if the self-trans
      // parameter has been set (default is 1.0, which does nothing).
      double initial_self_alpha;
      if (topic < params_->num_topics()) {
        initial_self_alpha = gsl_vector_get(alpha.get(), topic);
        gsl_vector_set(alpha.get(), topic, initial_self_alpha);
      }

      gsl_vector_view nu_view = gsl_matrix_row(vars_->nu_.get(), topic);
      gsl_vector* nu = &nu_view.vector;

      p.topic = topic;

      nu_lhood_ += nu_lhood(nu, &p);
      if (topic < params_->num_topics()) {
        gsl_vector_set(alpha.get(), topic, initial_self_alpha);
      }
    }

    for (int ii = 0; ii < params_->num_topics(); ++ii) {
      (*emission_counts)["nu_" + SimpleItoa(nu_coordinate_) +
                         "_" + SimpleItoa(ii)] =
        gsl_matrix_get(vars_->nu_.get(), nu_coordinate_, ii);
      gsl_vector_set(nu_vector.get(), ii, gsl_matrix_get(vars_->nu_.get(),
                                                         nu_coordinate_, ii));
    }

    // ZKCOMMENT: what is this about?
    // (*emission_counts)["lhood_1_1"]=nu_lhood_;

    // (*emission_counts)["lhood_1_1"]=GlobalWeightTermSingle(nu_coordinate_,
    // vars_->nu_.get(), vars_->beta_.get(), params_->alpha_trans());
    (*emission_counts)["lhood_1_1"]=gsl_vector_get(vars_->lhood_.get(), 0);
  }

  if (tau_coordinate_ >= 0) {
    for (int ii = 0; ii < params_->vocab_size(); ++ii) {
      (*emission_counts)["tau_" + SimpleItoa(tau_coordinate_) + "_" +
                         SimpleItoa(ii)] =
        gsl_matrix_get(vars_->tau_.get(), tau_coordinate_, ii);
    }
  }

  if (lhood_) {
    (*emission_counts)["lhood_final"] = gsl_vector_get(vars_->lhood_.get(), 0);
    // cout << "emission counts is "
    // << (*emission_counts)["lhood_final"] << endl;
  }
}

void SyntopStandaloneReducer::Optimize() {
  if (nu_coordinate_ >= 0) {
    if (check_order_) {
      assert(tau_coordinate_ == -1);
      assert(!lhood_);
    }

    if (!params_->ignore_trans()) UpdateNu();
  }

  if (tau_coordinate_ >= 0) {
    if (check_order_) {
      assert(nu_coordinate_ == -1);
      assert(!lhood_);
    }

    // display_matrix(vars_->tau_est_top_.get(), "tau_est_top is\n");
    // display_vector(vars_->tau_est_bottom_.get(), "tau_est_bottom is\n");

    for (int ii = 0; ii < params_->vocab_size(); ++ii) {
      double val = params_->vocab_sigma() / params_->vocab_size();
      val += gsl_matrix_get(vars_->tau_est_top_.get(), tau_coordinate_, ii);
      gsl_matrix_set(vars_->tau_.get(), tau_coordinate_, ii,
                     log(val) - log(params_->vocab_sigma() +
                                    gsl_vector_get(vars_->tau_est_bottom_.get(),
                                                   tau_coordinate_)));
    }
  }

  if (lhood_) {
    assert(tau_coordinate_ == -1);
    assert(nu_coordinate_ == -1);

    // Don't have to do anything
    return;
  }
}

double SyntopStandaloneReducer::GlobalWeightTermPrior(const gsl_vector* beta,
                                                      double alpha_top,
                                                      bool finite) {
  double prior = 0.0;
  if (!finite) {
    // The second argument is NULL because we don't want to get a full
    // vector out
    shared_ptr<gsl_vector> T(gsl_vector_alloc(beta->size),
                             gsl_vector_free);
    create_T_full_beta(beta, beta->size, NULL, T.get());
    prior += beta_prior_lhood(T.get(), alpha_top);
  }
  return prior;
}

double SyntopStandaloneReducer::GlobalWeightTermSingle(int i,
                                                       const gsl_matrix* nu,
                                                       const gsl_vector* beta,
                                                       double alpha_trans) {
  double shared = 0.0;
  double not_shared = 0.0;

  gsl_vector_const_view nu_i = gsl_matrix_const_row(nu, i);
  double a;

  // We don't do extra work if there is no self-loop penalty or
  // we're in the start state
  a = var_dir_lhood(alpha_trans, beta, &nu_i.vector);

  double b = var_dir_lhood(&nu_i.vector, &nu_i.vector);
  shared += a;
  not_shared -= b;

  return shared + not_shared;
}

double SyntopStandaloneReducer::GlobalWeightTerm(const gsl_matrix* nu,
                                                 const gsl_vector* beta,
                                                 double alpha_trans,
                                                 double alpha_top,
                                                 bool finite) {
  double shared = 0.0;
  double not_shared = 0.0;
  double prior = 0.0;
  for (unsigned int i = 0; i < nu->size1; i++) {
    gsl_vector_const_view nu_i = gsl_matrix_const_row(nu, i);
    double a;

    // We don't do extra work if there is no self-loop penalty or
    // we're in the start state
    a = var_dir_lhood(alpha_trans, beta, &nu_i.vector);

    double b = var_dir_lhood(&nu_i.vector, &nu_i.vector);
    shared += a;
    not_shared -= b;
  }

  if (!finite) {
    // The second argument is NULL because we don't want to get a full
    // vector out
    shared_ptr<gsl_vector> T(gsl_vector_alloc(beta->size),
                             gsl_vector_free);
    create_T_full_beta(beta, beta->size, NULL, T.get());
    prior += beta_prior_lhood(T.get(), alpha_top);

    // cout << "beta_lhood_prior is " << prior << endl;
  }
  // cout << "BETA NU:" << shared << " BETA PRIOR: " << prior << endl;
  return shared + not_shared + prior;
}

/*
 * Compute the contribution of vocabulary to lhood
 */
double SyntopStandaloneReducer::GlobalVocabularyTerm(const gsl_matrix* tau,
                                                     const double vocab_sigma,
                                                     const int number_terms,
                                                     const int number_topics) {
  double vocab_lhood = 0.0;
  for (int v = 0; v < number_terms; v++) {
    for (int i = 0; i < number_topics; i++) {
      vocab_lhood -= gsl_matrix_get(tau, i, v) * gsl_matrix_get(tau, i, v);
    }
  }

  vocab_lhood = vocab_lhood / (2.0 * vocab_sigma * vocab_sigma);

  return vocab_lhood;
}

/*
 * While most of our likelihood terms are document-specific, some
 * likelihood terms cannot be expressed as a sum over document
 * likelihood terms.  Such terms are computed here for the current
 * state of variable assignments.
 */
double SyntopStandaloneReducer::GlobalLikelihoodTerm() {
  double weights = GlobalWeightTerm(vars_->nu_.get(),
                                    vars_->beta_.get(),
                                    params_->alpha_trans(),
                                    params_->alpha_top(),
                                    params_->finite());

  double vocab = GlobalVocabularyTerm(vars_->tau_.get(),
                                      params_->vocab_sigma(),
                                      params_->vocab_size(),
                                      params_->num_topics());

  return weights + vocab;
}

/*
 * Should only be called if we don't have interaction terms to gum up
 * the works.  Allows you to do closed form updates.
 */
void SyntopStandaloneReducer::LazyUpdateNu(int topic) {
  if (topic % 16 == 0) {
    cout << "Closed Nu update " << topic << endl;
  }

  // display_matrix(vars_->nu_.get(), "nu before updating is:");
  for (int j = 0; j < params_->num_topics(); j++) {
    double alpha = gsl_vector_get(vars_->beta_.get(), j) *
      params_->alpha_trans();

    gsl_matrix_set(vars_->nu_.get(), topic, j,
                   alpha + gsl_matrix_get(vars_->nu_exp_topic_counts_.get(),
                                          topic, j));
  }
}

void SyntopStandaloneReducer::UpdateNu() {
  if (params_->shortcut_gsl() && params_->ignore_docs()) {
    LazyUpdateNu(nu_coordinate_);
  } else {
    optimize_nu(vars_, params_.get(), nu_coordinate_);
  }
}

/*
double GlobalOptimizer::VocParUpdate(double numerator,
                                     double denominator,
                                     double sigma) {
  double c = 1.0 / (sigma * sigma);
  double z = denominator * safe_exp(numerator / c) / c;
  double answer = 0.0;

  if(sigma == 0.0) {
    answer = safe_log(numerator) - safe_log(denominator);
    return answer;
  } else if(z == numeric_limits<double>::infinity()) {
    // z is asymptotic to log at this point, so we can just work in log space
    answer = c * (log(denominator) + (numerator / c) - log(c));
    cout << "IFFY";
  } else if(z >= 0.0) {
    assert(check_valid_double(z, "vocab lambert input"));
    answer = c * gsl_sf_lambert_W0(z);
  } else {
    cout << "INVALID Lambert argument: " << z << endl;
  }
  assert(z >= 0.0);

  answer = (numerator - answer) / c;
  return answer;
}

void LensModel::UpdateRho() {
  for(int t=0; t<corpus_->get_number_roles_(); t++) {
    for(int v=0; v<corpus_->get_number_terms_(); v++) {
      if(v % 500 == 0) {
        cout << "Updating rho,  role=" << t << ", word=" << v << endl;
      }

      double before;
      if(check_lhood_) {
        before = GlobalVocabularyTerm(tau_.get(),
                                      rho_.get(),
                                      vocab_sigma_,
                                      corpus_->get_number_terms_(),
                                      corpus_->get_number_roles_(),
                                      number_topics_,
                                      ignore_tags_);
        before += DocumentRhoLikelihood(rho_.get(),
                                        rho_est_top_.get(),
                                        rho_est_bottom_.get());
      }

      double rho = VocParUpdate(gsl_matrix_get(rho_est_top_.get(), t, v),
                                gsl_matrix_get(rho_est_bottom_.get(), t, v),
                                vocab_sigma_);

      gsl_matrix_set(rho_.get(), t, v, rho);
      assert(check_valid_double(rho, "rho_{" + SimpleItoa(t) + ","
                                      + SimpleItoa(v) + "}"));

      double after;
      if(check_lhood_) {
        after = GlobalVocabularyTerm(tau_.get(),
                                     rho_.get(),
                                     vocab_sigma_,
                                     corpus_->get_number_terms_(),
                                     corpus_->get_number_roles_(),
                                     number_topics_,
                                     ignore_tags_);
        after += DocumentRhoLikelihood(rho_.get(),
                                       rho_est_top_.get(),
                                       rho_est_bottom_.get());

        check_increase(("rho_{" + SimpleItoa(t) + "," +
                        SimpleItoa(v) + "}").c_str(),
                       lhood_file_,
                       after,
                       before);

      }
    }
  }
}
*/
