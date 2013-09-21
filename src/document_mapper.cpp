/*
 * Copyright Jordan Boyd-Graber 2008
 *
 */

#include "document_mapper.h"
// This is included here to make life easier
#include "gradient.h"

using std::cerr;

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
DocumentMapper::DocumentMapper(VariationalParameters* vars,
                               const SyntopParameters* params) {
  // Some things don't make sense with a limited number of topics
  BOOST_ASSERT(params->num_topics() > 0);
  BOOST_ASSERT(params->finite() || params->num_topics() > 2);

  vars_ = vars;
  params_ = params;
}

/*
double DocumentMapper::WordProbabilityUnnormalized(int v,
                                                   int role,
                                                   int num_topics,
                                                   gsl_vector* phi,
                                                   gsl_matrix* tau) {
  double val = 0.0;
  for(int i=0; i<num_topics; i++) {
    val += gsl_vector_get(phi, i) * exp(gsl_matrix_get(tau, i, v));
  }
  return val;
}
*/

void DocumentMapper::UpdatePhi(Document* d) {
  // Updating these variables are part and parcel with updating phi,
  // and we want to do it uniformly, so we do it here
  UpdateOmega(d);
  // display_vector(d->get_omega_(), "omega a is\n");

  shared_ptr<gsl_vector> temp_phi_ptr(gsl_vector_alloc(params_->num_topics()),
                                      gsl_vector_free);
  gsl_vector_set_all(temp_phi_ptr.get(), 0.0);
  gsl_vector* temp_phi = temp_phi_ptr.get();

  int doc = d->get_id();
  int num_words = d->get_number_words_();
  gsl_vector_view gamma_view = gsl_matrix_row(vars_->gamma_.get(), doc);
  gsl_vector* gamma = &gamma_view.vector;

  for (int n = 0; n < num_words; n++) {
    LensNode* w = d->get_word_(n);
    for (int i = 0; i < params_->num_topics(); i++) {
      double new_phi = ComputeNewPhi(d,
                                     n,
                                     i,
                                     vars_->nu_.get(),
                                     vars_->tau_.get(),
                                     gamma,
                                     params_->vocab_size(),
                                     params_->num_topics());
      gsl_vector_set(temp_phi, i, new_phi);
    }

    // display_vector(w->phi_.get(), "phi before updating\n");
    w->ReplacePhiByUnnormalizedLogSpace(temp_phi);
    // display_vector(w->phi_.get(), "phi after updating\n");

    /*
    if(check_lhood_) {
      double new_lhood = Likelihood(doc);
      check_increase((SimpleItoa(n) + "'s phi").c_str(),
                     lhood_file_,
                     new_lhood,
                     lhood);
    }
    */

    /*
    display_vector(w->phi_.get(), ("NEW phi" + SimpleItoa(doc) + "," +
                                   SimpleItoa(n)).c_str());
    */
  }

  // Not sure this is needed, but we want it to be accurate because it
  // depends heavily on phi
  /*
  display_vector(d->omega_.get(), ("NEW omega" + SimpleItoa(doc)).c_str());
  display_vector(d->xi_.get(), ("NEW xi" + SimpleItoa(doc)).c_str());
  */
  UpdateOmega(d);
  // display_vector(d->get_omega_(), "omega a is\n");
}

void DocumentMapper::Emit(Document* d, StringMap* count) {
  int doc = d->get_id();

  for (int ii = 0; ii < params_->num_topics(); ++ii) {
    // Tau
    string key = "tau_" + SimpleItoa(ii) + "_-1";
    (*count)[key] += gsl_vector_get(vars_->tau_est_bottom_.get(), ii);

    for (int vv = 0; vv < params_->vocab_size(); ++vv) {
      double val = gsl_matrix_get(vars_->tau_est_top_.get(), ii, vv);
      if (val > 0.0) {
        key = "tau_" + SimpleItoa(ii) + "_" + SimpleItoa(vv);
        (*count)[key] += val;
      }
    }

    // Gamma
    key = "gamma_" + SimpleItoa(doc) + "_" + SimpleItoa(ii);
    (*count)[key] += gsl_matrix_get(vars_->gamma_.get(), doc, ii);

    // Nu ss
    for (int jj = 0; jj < params_->num_topics() + 1; ++jj) {
      // It's jj before ii because I wanted to put it inside the outer loop, but
      // the first index of nu ss needs to go from 0 to K instead of 0 to K-1
      string key_doc = "nu_" + SimpleItoa(jj) + "_" + SimpleItoa(ii) + "_doc";
      string key_trans = "nu_" + SimpleItoa(jj) + "_" + SimpleItoa(ii) +
        "_trans";
      (*count)[key_doc] += gsl_matrix_get(vars_->nu_exp_doc_counts_.get(),
                                         jj, ii);
      (*count)[key_trans] += gsl_matrix_get(vars_->nu_exp_topic_counts_.get(),
                                           jj, ii);

      string key_flag = "nu_" + SimpleItoa(jj) + "_-1";
      (*count)[key_flag] += 0;
    }
  }
  // cout << "emit all parameters..." << endl;
  (*count)["lhood_1_1"] = Likelihood(d);
}

void DocumentMapper::UpdateOmega(Document* d) {
  int num_words = d->get_number_words_();
  gsl_vector_view gamma_view = gsl_matrix_row(vars_->gamma_.get(), d->get_id());
  gsl_vector* gamma = &gamma_view.vector;
  double gamma_normalizer = gsl_blas_dsum(gamma);

  /*
   * If gamma is always one, then we don't have to worry about the
   * topic normalizer, and thus the slack variable for the log sum.
   * Thus, omega can always be zero and we can ignore the update.
   */
  if (params_->ignore_docs() || params_->ignore_trans()) {
    gsl_vector_set_all(d->get_omega_(), 0.0);
    return;
  }

  for (int n = 0; n < num_words; n++) {
    LensNode* w = d->get_word_(n);
    int parent = w->get_parent_index_();
    gsl_vector* parent_phi = NULL;
    if (parent != -1) {
      parent_phi = d->get_word_(parent)->phi_.get();
    }

    // The likelihood computation igores omega if it's zero, which
    // it starts out as.  Thus, we have to not compute it initially.

    double new_omega = 0;
    for (int i = 0; i < params_->num_topics(); i++) {
      double contribution = 1;
      if (parent != -1) {
        contribution = 0;
        for (int j = 0; j < params_->num_topics(); j++) {
          gsl_vector_view nu_j = gsl_matrix_row(vars_->nu_.get(), j);
          double nu_j_sum = gsl_blas_dsum(&nu_j.vector);

          contribution += gsl_vector_get(parent_phi, j) *
            (gsl_matrix_get(vars_->nu_.get(), j, i) /
             nu_j_sum);
        }
      } else {
        contribution = 0;
        gsl_vector_view nu_j = gsl_matrix_row(vars_->nu_.get(),
                                              params_->num_topics());
        double nu_j_sum = gsl_blas_dsum(&nu_j.vector);

        contribution += (gsl_matrix_get(vars_->nu_.get(),
                                        params_->num_topics(), i) / nu_j_sum);
      }
      contribution *= gsl_vector_get(gamma, i) / gamma_normalizer;
      new_omega += contribution;
    }

    BOOST_ASSERT(check_valid_double(new_omega, "omega_{" +
                                    SimpleItoa(d->get_id()) + "," +
                                    SimpleItoa(n) + "}"));
    BOOST_ASSERT(new_omega > 0);

    gsl_vector_set(d->get_omega_(), n, 1.0/new_omega);
  }
}

/*
 * We only save the information we need to update the other parameters
 * on the last step, when update_params is set to true.
 */
double DocumentMapper::EStep(double tolerance, Document* d, int doc) {
  double likelihood = -1e256;
  double new_likelihood = -1e255;

  gsl_vector_view gamma_view = gsl_matrix_row(vars_->gamma_.get(), doc);
  gsl_vector* gamma = &gamma_view.vector;

  shared_ptr<gsl_vector>
    phi_total(gsl_vector_alloc(params_->num_topics()), gsl_vector_free);
  shared_ptr<gsl_vector>
    nu_interaction(gsl_vector_alloc(params_->num_topics()), gsl_vector_free);

  int doc_iter = 0;
  double diff = fabs((likelihood - new_likelihood));
  assert(diff > 0);

  while (diff > tolerance &&
         doc_iter++ < params_->max_doc_iterations()) {
    double old_likelihood_test = phi_lhood(d);
    // display_matrix();
    UpdatePhi(d);
    // display_matrix();
    double new_likelihood_test = phi_lhood(d);
    if (old_likelihood_test > new_likelihood_test) {
      cerr << "Something wrong with phi updating: " <<
        old_likelihood_test << "\t" << new_likelihood_test << endl;
    }

    // Now we use conjugate gradient to compute gamma
    // We don't do it on every iteration because it's costly
    if (!params_->ignore_docs()) {
      CumulateGammaStatistics(d,
                              vars_->nu_.get(),
                              phi_total.get(),
                              nu_interaction.get());

      old_likelihood_test = gamma_likelihood(d);
      if (params_->shortcut_gsl() && params_->ignore_trans()) {
        for (int i = 0; i < params_->num_topics(); i++) {
          gsl_vector_set(gamma, i, gsl_vector_get(vars_->beta_.get(), i) *
                         params_->alpha_doc() +
                         gsl_vector_get(phi_total.get(), i));
        }
      } else {
        optimize_gamma(gamma, phi_total.get(),
                       nu_interaction.get(), params_->model_name(),
                       vars_->beta_.get(), params_->alpha_doc(), doc);
      }
      new_likelihood_test = gamma_likelihood(d);

      if (new_likelihood_test-old_likelihood_test < -0.000001) {
        cerr << "Something wrong with gamma updating: " <<
          old_likelihood_test << "\t" << new_likelihood_test << endl;
      }

      /*
      if(check_lhood_) {
        double sum = dot_product(nu_interaction.get(), gamma) /
          gsl_blas_dsum(gamma);
        double actual = gamma_nu_interaction(d, gamma, nu_.get(), false);
        double diff = abs(sum - actual);

        cout << "INTERACTION DIFF:" << sum << " - " << actual << " = " <<
          diff << endl;
        assert(diff < ERROR_TOLERANCE);

        sum = 0.0;
        for(unsigned i=0; i < phi_total->size; i++) {
          double contribution = gsl_vector_get(phi_total.get(), i) *
            (gsl_sf_psi(gsl_vector_get(gamma, i)) -
             gsl_sf_psi(gsl_blas_dsum(gamma)));
          sum += contribution;
        }
        actual = gamma_lhood_exp_counts(d, gamma);
        diff = abs(sum - actual);

        cout << "EX COUNTS DIFF:" << sum << " - " << actual << " = " <<
          diff << endl;
        assert(diff < ERROR_TOLERANCE);

        likelihood = check_increase("Gamma_" + SimpleItoa(doc),
                                    lhood_file_,
                                    Likelihood(doc),
                                    likelihood);
      }
      */
    } else {
      cout << "Ignoring gamma update." << endl;
    }

    likelihood = new_likelihood;
    new_likelihood = Likelihood(d);
    diff = fabs((likelihood - new_likelihood));

    // shared_ptr<gsl_vector> lhood_vector(gsl_vector_alloc(7),
    // gsl_vector_free); gsl_vector_set_all(lhood_vector.get(), 0);
    // new_likelihood = Likelihood(d, lhood_vector.get());
    // display_vector(lhood_vector.get(), "updated lhood is ");

    /*
    check_increase("doc_" + SimpleItoa(doc) + " E-step iteration " +
                   SimpleItoa(doc_iter),
                   lhood_file_,
                   Likelihood(doc),
                   likelihood);
    */
  }

  CumulateNuStatistics(d,
                       doc,
                       gamma,
                       vars_->nu_exp_topic_counts_.get(),
                       vars_->nu_exp_doc_counts_.get());

  CumulateTauUpdate(d,
                    params_->num_topics(),
                    params_->vocab_size(),
                    vars_->tau_est_top_.get(),
                    vars_->tau_est_bottom_.get());

  return Likelihood(d, vars_->lhood_.get());
}

/*
 * It is possible to pass tau_est_top as NULL, which doesn't update
 * the numerator.  This is hoped to speed up the compute time when we
 * do multiple iterations during the M step.
 */
void DocumentMapper::CumulateTauUpdate(Document* d,
                                       const int number_topics,
                                       const int number_terms,
                                       gsl_matrix* tau_est_top,
                                       gsl_vector* tau_est_bottom) {
  int num_words = d->get_number_words_();

  if (tau_est_top != NULL) {
    for (int i = 0; i < number_topics; i++) {
      for (int n = 0; n < num_words; n++) {
        LensNode* w = d->get_word_(n);
        gsl_matrix_set(tau_est_top, i, w->term_,
                       gsl_matrix_get(tau_est_top, i, w->term_) +
                       gsl_vector_get(w->phi_.get(), i));
      }
    }
  }

  for (int i = 0; i < number_topics; i++) {
    double new_tau_bottom = 0.0;
    for (int n = 0; n < num_words; n++) {
      LensNode* w = d->get_word_(n);
      double contribution = gsl_vector_get(w->phi_.get(), i);
      new_tau_bottom += contribution;
      /*
         cout << "TB:" << i << "," << v << " " << contribution << "<- phi="
         << gsl_vector_get(w->phi_.get(), i) << " rho="
         << gsl_matrix_get(rho, w->relation_, v) << " xi="
         << 1/gsl_vector_get(d->xi_.get(), n) << endl;
      */
    }

    // display_matrix(tau_est_top, "tau_est_top is\n");

    double check = 0.0;
    for (int v = 0; v < number_terms; v++) {
      check += gsl_matrix_get(tau_est_top, i, v);

      // KE ZHAI: move this block outside the for loop
      // gsl_vector_set(tau_est_bottom,
      //           i,
      //           gsl_vector_get(tau_est_bottom, i) +
      //           new_tau_bottom);
    }

    gsl_vector_set(tau_est_bottom, i,
                   gsl_vector_get(tau_est_bottom, i) +
                   new_tau_bottom);

    // cout << "i value is " << i << " and check value is " << check << endl;
    // display_vector(tau_est_bottom, "tau_est_bottom is\n");

    assert(abs(check - gsl_vector_get(tau_est_bottom, i)) < ERROR_TOLERANCE);
    // assert(abs(check - gsl_matrix_get(tau_est_bottom, i, 0)) <
    // ERROR_TOLERANCE);
  }
}

double DocumentMapper::LikelihoodTopic(Document* d,
                                  const int num_topics) {
  double val = 0.0;
  int num_words = d->get_number_words_();
  for (int n = 0; n < num_words; n++) {
    LensNode* w = d->get_word_(n);
    for (int i = 0; i < num_topics; i++) {
      // E_q [ log q(z) ]
      val += gsl_vector_get(w->phi_.get(), i) *
        safe_log(gsl_vector_get(w->phi_.get(), i));
    }
  }
  return val;
}

double DocumentMapper::LikelihoodVocab(Document* d,
                                       const gsl_matrix* tau,
                                       const int number_topics,
                                       const int number_terms) {
  double word_expectation = 0.0;
  double normalizer = 0.0;
  int number_words = d->get_number_words_();
  for (int n = 0; n < number_words; n++) {
    LensNode* w = d->get_word_(n);


    for (int i = 0; i < number_topics; i++) {
      // E_q log p(w | z)

      double phi_i = gsl_vector_get(w->phi_.get(), i);
      double tau_iw = gsl_matrix_get(tau, i, w->term_);

      BOOST_ASSERT(check_valid_double(phi_i, "phi"));
      BOOST_ASSERT(check_valid_double(tau_iw, "tau"));

      // This means that we have no roles
      word_expectation += phi_i * (tau_iw);
    }
  }

  return word_expectation - normalizer;
}

double DocumentMapper::Likelihood(Document* d,
                                  const gsl_vector* gamma,
                                  const double alpha_doc,
                                  const gsl_vector* beta,
                                  const gsl_matrix* nu,
                                  const int num_topics,
                                  const int num_terms,
                                  const gsl_matrix* tau,
                                  gsl_vector* lhood_total) {
  double val = 0.0;

  // E_q [ p (\theta | \alpha) ]
  double term1 = var_dir_lhood(alpha_doc, beta, gamma);
  val += term1;
  if (lhood_total) {
    vct_accum(0, term1, lhood_total);
  }

  // - E_q [ log q(\theta) ]
  double term2 = var_dir_lhood(gamma, gamma);
  val -= term2;
  if (lhood_total) {
    // vct_accum(1, -term2, lhood_total);
    vct_accum(0, -term2, lhood_total);
  }

  // E_q [ log p(z | \pi, \theta) ]
  double term3 = gamma_lhood_exp_counts(d, gamma);
  val += term3;
  if (lhood_total) {
    // vct_accum(2, term3, lhood_total);
    vct_accum(0, term3, lhood_total);
  }

  double term4 = 0;
  for (int i = 0; i <= num_topics; i++) {
    double row_contrib = nu_doc_transition_likelihood(d,
      &gsl_matrix_const_row(nu, i).vector, i);

    term4 += row_contrib;
  }
  val += term4;
  if (lhood_total) {
    // vct_accum(3, term4, lhood_total);
    vct_accum(0, term4, lhood_total);
  }

  // double term5 = 0;
  double term5 = gamma_nu_interaction(d, gamma, nu);
  // Ignore this term if omega is zero (which should only happen if
  // we've durned of document or transitions)
  if (gsl_vector_max(d->get_omega_()) > 0) {
    val -= term5;
    if (lhood_total) {
      // vct_accum(4, -term5, lhood_total);
      vct_accum(0, -term5, lhood_total);
    }
  }

  double term6 = LikelihoodVocab(d,
                                 tau,
                                 num_topics,
                                 num_terms);

  val += term6;
  if (lhood_total) {
    // vct_accum(5, term6, lhood_total);
    vct_accum(0, term6, lhood_total);
  }

  double term7 = LikelihoodTopic(d,
                                 num_topics);
  val -= term7;
  if (lhood_total) {
    // vct_accum(6, -term7, lhood_total);
    vct_accum(0, -term7, lhood_total);
  }

  // cout << "likelihood DOC LHOOD " << term1 << " + " << -term2 << " + "
  // << term3 << " + " << term4 << " + " << -term5 << " + " << term6
  // << " + " << -term7 << " = " << val << endl;

  BOOST_ASSERT(check_valid_double(term1, "Lhood term1"));
  BOOST_ASSERT(check_valid_double(term2, "Lhood term2"));
  BOOST_ASSERT(check_valid_double(term3, "Lhood term3"));
  BOOST_ASSERT(check_valid_double(term4, "Lhood term4"));
  BOOST_ASSERT(check_valid_double(term5, "Lhood term5"));
  BOOST_ASSERT(check_valid_double(term6, "Lhood term6"));
  BOOST_ASSERT(check_valid_double(term7, "Lhood term7"));

  return val;
}

double DocumentMapper::phi_lhood(Document* d,
                                 const gsl_vector* gamma,
                                 const double alpha_doc,
                                 const gsl_vector* beta,
                                 const gsl_matrix* nu,
                                 const int num_topics,
                                 const int num_terms,
                                 const gsl_matrix* tau,
                                 gsl_vector* lhood_total) {
  double val = 0.0;

  // E_q [ p (\theta | \alpha) ]
  double term1 = var_dir_lhood(alpha_doc, beta, gamma);
  term1 = 0;
  val += term1;
  if (lhood_total) {
    vct_accum(0, term1, lhood_total);
  }

  // - E_q [ log q(\theta) ]
  double term2 = var_dir_lhood(gamma, gamma);
  term2 = 0;
  val -= term2;
  if (lhood_total) {
    // vct_accum(1, -term2, lhood_total);
    vct_accum(0, -term2, lhood_total);
  }


  // E_q [ log p(z | \pi, \theta) ]
  double term3 = gamma_lhood_exp_counts(d, gamma);
  val += term3;
  if (lhood_total) {
    // vct_accum(2, term3, lhood_total);
    vct_accum(0, term3, lhood_total);
  }

  double term4 = 0;
  for (int i = 0; i <= num_topics; i++) {
    double row_contrib = nu_doc_transition_likelihood(d,
                                                      &gsl_matrix_const_row
                                                      (nu, i).vector,
                                                      i);

    term4 += row_contrib;
  }
  val += term4;
  if (lhood_total) {
    // vct_accum(3, term4, lhood_total);
    vct_accum(0, term4, lhood_total);
  }

  double term5 = gamma_nu_interaction(d, gamma, nu);
  // Ignore this term if omega is zero (which should only happen if
  // we've durned of document or transitions)
  if (gsl_vector_max(d->get_omega_()) > 0) {
    val -= term5;
    if (lhood_total) {
      // vct_accum(4, -term5, lhood_total);
      vct_accum(0, -term5, lhood_total);
    }
  }

  double term6 = LikelihoodVocab(d,
                                 tau,
                                 num_topics,
                                 num_terms);

  val += term6;
  if (lhood_total) {
    // vct_accum(5, term6, lhood_total);
    vct_accum(0, term6, lhood_total);
  }

  double term7 = LikelihoodTopic(d,
                                 num_topics);
  val -= term7;
  if (lhood_total) {
    // vct_accum(6, -term7, lhood_total);
    vct_accum(0, -term7, lhood_total);
  }

  // cout << "phi DOC LHOOD " << term1 << " + " << -term2 << " + "
  // << term3 << " + " << term4 << " + " << -term5 << " + "
  // << term6 << " + " << -term7 << " = " << val << endl;

  BOOST_ASSERT(check_valid_double(term1, "Lhood term1"));
  BOOST_ASSERT(check_valid_double(term2, "Lhood term2"));
  BOOST_ASSERT(check_valid_double(term3, "Lhood term3"));
  BOOST_ASSERT(check_valid_double(term4, "Lhood term4"));
  BOOST_ASSERT(check_valid_double(term5, "Lhood term5"));
  BOOST_ASSERT(check_valid_double(term6, "Lhood term6"));
  BOOST_ASSERT(check_valid_double(term7, "Lhood term7"));

  return val;
}

double DocumentMapper::gamma_likelihood(Document* d,
                                        const gsl_vector* gamma,
                                        const double alpha_doc,
                                        const gsl_vector* beta,
                                        const gsl_matrix* nu,
                                        const int num_topics,
                                        const int num_terms,
                                        const gsl_matrix* tau,
                                        gsl_vector* lhood_total) {
  double val = 0.0;

  // E_q [ p (\theta | \alpha) ]
  double term1 = var_dir_lhood(alpha_doc, beta, gamma);
  val += term1;
  if (lhood_total) {
    vct_accum(0, term1, lhood_total);
  }

  // - E_q [ log q(\theta) ]
  double term2 = var_dir_lhood(gamma, gamma);
  val -= term2;
  if (lhood_total) {
    // vct_accum(1, -term2, lhood_total);
    vct_accum(0, -term2, lhood_total);
  }

  // E_q [ log p(z | \pi, \theta) ]
  double term3 = gamma_lhood_exp_counts(d, gamma);
  val += term3;
  if (lhood_total) {
    // vct_accum(2, term3, lhood_total);
    vct_accum(0, term3, lhood_total);
  }

  double term4 = 0;
  /*
  for(int i=0; i<=num_topics; i++) {
    double row_contrib = nu_doc_transition_likelihood(d,
                                                      &gsl_matrix_const_row
                                                      (nu, i).vector,
                                                      i);

    term4 += row_contrib;
  }
  val += term4;
  if(lhood_total) {
    //vct_accum(3, term4, lhood_total);
    vct_accum(0, term4, lhood_total);
  }
  */

  double term5 = gamma_nu_interaction(d, gamma, nu);
  // Ignore this term if omega is zero (which should only happen if
  // we've durned of document or transitions)
  if (gsl_vector_max(d->get_omega_()) > 0) {
    val -= term5;
    if (lhood_total) {
      // vct_accum(4, -term5, lhood_total);
      vct_accum(0, -term5, lhood_total);
    }
  }

  double term6 = LikelihoodVocab(d,
                                 tau,
                                 num_topics,
                                 num_terms);
  term6 = 0;
  val += term6;
  if (lhood_total) {
    // vct_accum(5, term6, lhood_total);
    vct_accum(0, term6, lhood_total);
  }

  double term7 = LikelihoodTopic(d,
                                 num_topics);
  term7 = 0;
  val -= term7;
  if (lhood_total) {
    // vct_accum(6, -term7, lhood_total);
    vct_accum(0, -term7, lhood_total);
  }

  // cout << "gamma DOC LHOOD " << term1 << " + " << -term2 << " + "
  // << term3 << " + " << term4 << " + " << -term5 << " + "
  // << term6 << " + " << -term7 << " = " << val << endl;

  BOOST_ASSERT(check_valid_double(term1, "Lhood term1"));
  BOOST_ASSERT(check_valid_double(term2, "Lhood term2"));
  BOOST_ASSERT(check_valid_double(term3, "Lhood term3"));
  BOOST_ASSERT(check_valid_double(term4, "Lhood term4"));
  BOOST_ASSERT(check_valid_double(term5, "Lhood term5"));
  BOOST_ASSERT(check_valid_double(term6, "Lhood term6"));
  BOOST_ASSERT(check_valid_double(term7, "Lhood term7"));

  return val;
}

// Does not compute terms that depend only on \nu, \beta, and vocab
// priors, which are global and not related to individual documents
double DocumentMapper::gamma_likelihood(Document* d, gsl_vector* lhood) {
  int doc = d->get_id();
  gsl_vector_view gamma_view = gsl_matrix_row(vars_->gamma_.get(), doc);
  gsl_vector* gamma = &gamma_view.vector;

  return gamma_likelihood(d,
                    gamma,
                    params_->alpha_doc(),
                    vars_->beta_.get(),
                    vars_->nu_.get(),
                    params_->num_topics(),
                    params_->vocab_size(),
                    vars_->tau_.get(),
                    lhood);
}

// Does not compute terms that depend only on \nu, \beta, and vocab
// priors, which are global and not related to individual documents
double DocumentMapper::phi_lhood(Document* d, gsl_vector* lhood) {
  int doc = d->get_id();
  gsl_vector_view gamma_view = gsl_matrix_row(vars_->gamma_.get(), doc);
  gsl_vector* gamma = &gamma_view.vector;

  return phi_lhood(d,
                    gamma,
                    params_->alpha_doc(),
                    vars_->beta_.get(),
                    vars_->nu_.get(),
                    params_->num_topics(),
                    params_->vocab_size(),
                    vars_->tau_.get(),
                    lhood);
}


// Does not compute terms that depend only on \nu, \beta, and vocab
// priors, which are global and not related to individual documents
double DocumentMapper::Likelihood(Document* d, gsl_vector* lhood) {
  int doc = d->get_id();
  gsl_vector_view gamma_view = gsl_matrix_row(vars_->gamma_.get(), doc);
  gsl_vector* gamma = &gamma_view.vector;

  return Likelihood(d,
                    gamma,
                    params_->alpha_doc(),
                    vars_->beta_.get(),
                    vars_->nu_.get(),
                    params_->num_topics(),
                    params_->vocab_size(),
                    vars_->tau_.get(),
                    lhood);
}

void DocumentMapper::CumulateGammaStatistics(Document* d,
                                        const gsl_matrix* nu,
                                        gsl_vector* phi_total_vector,
                                        gsl_vector* nu_interaction_vector) {
  gsl_vector_set_all(phi_total_vector, 0.0);
  gsl_vector_set_all(nu_interaction_vector, 0.0);

  int number_words = d->get_number_words_();
  unsigned int number_topics = phi_total_vector->size;
  for (unsigned int i = 0; i < number_topics; i++) {
    double phi_sum = 0.0;
    double nu_interaction = 0.0;

    for (int n = 0; n < number_words; n++) {
      LensNode * w = d->get_word_(n);
      phi_sum += gsl_vector_get(w->phi_.get(), i);

      int parent = w->get_parent_index_();
      double weighted_nu = 0.0;
      if (parent >= 0) {
        for (unsigned int j = 0; j < number_topics; j++) {
          weighted_nu += gsl_matrix_get(nu, j, i) / gsl_matrix_rowsum(nu, j) *
            gsl_vector_get(d->get_word_(parent)->phi_.get(), j);
        }
      } else {
        weighted_nu += gsl_matrix_get(nu, number_topics, i) /
          gsl_matrix_rowsum(nu, number_topics);
      }
      nu_interaction += weighted_nu * d->get_omega_(n);
    }

    gsl_vector_set(phi_total_vector, i, phi_sum);
    gsl_vector_set(nu_interaction_vector, i, nu_interaction);
  }
}


void DocumentMapper::CumulateNuStatistics(Document* d,
                                          const int doc,
                                          const gsl_vector* gamma,
                                          gsl_matrix* nu_exp_topic_counts,
                                          gsl_matrix* nu_exp_doc_counts) {
  // In order to make computing updates for nu easier, we store extra
  // information so we don't have to do expensive sums for each
  // iteraction
  int num_words = d->get_number_words_();
  int number_topics = gamma->size;

  double gamma_sum = gsl_blas_dsum(gamma);

  // display_vector(gamma, "gamma before cumulate nu statis is\n");
  // display_matrix(nu_exp_topic_counts, "nu_exp_trans_counts before cumulate nu
  // statis is\n"); display_matrix(nu_exp_doc_counts, "nu_exp_doc_counts before
  // cumulate nu statis is\n");

  for (int n = 0; n < num_words; n++) {
    LensNode* w = d->get_word_(n);
    int parent = w->get_parent_index_();
    gsl_vector* parent_phi = NULL;
    if (parent != -1) {
      parent_phi = d->get_word_(parent)->phi_.get();
    }

    // Update expected transition counts, which will be used by nu
    // update later
    if (parent != -1) {
      for (int i = 0; i < number_topics; i++) {
        for (int j = 0; j < number_topics; j++) {
          double contribution = gsl_vector_get(parent_phi, i) *
            (gsl_vector_get(gamma, j) / gamma_sum) * d->get_omega_(n);
          mtx_accum(i, j, contribution, nu_exp_doc_counts);

          double trans = gsl_vector_get(parent_phi, i) *
            gsl_vector_get(w->phi_.get(), j);
          mtx_accum(i, j, trans, nu_exp_topic_counts);
        }
      }
    } else {
      // If there is no parent, it's the head of a tree, and we don't
      // have to do a second expectation over the parent states

      for (int i = 0; i < number_topics; i++) {
        double contribution = (gsl_vector_get(gamma, i) / gamma_sum) *
          d->get_omega_(n);
        mtx_accum(number_topics, i,
                  contribution,
                  nu_exp_doc_counts);
        mtx_accum(number_topics, i,
                  gsl_vector_get(w->phi_.get(), i),
                  nu_exp_topic_counts);
      }
    }
  }
}

// These functions are static because it's easier to test that way;
// there's probably a better way
double DocumentMapper::ComputeNewPhiIntrinsic(const int topic,
                                              const gsl_vector* gamma,
                                              double gamma_normalizer) {
  double new_phi = gsl_sf_psi(gsl_vector_get(gamma, topic)) -
    gsl_sf_psi(gamma_normalizer);
  BOOST_ASSERT(!isnan(new_phi));
  return new_phi;
}

double DocumentMapper::ComputeNewPhiParent(const gsl_vector* parent_phi,
                                           const gsl_matrix* nu,
                                           const int topic,
                                           const int number_topics) {
  BOOST_ASSERT(parent_phi);
  double new_phi = 0.0;
  for (int j = 0; j < number_topics; j++) {
    gsl_vector_const_view nu_j = gsl_matrix_const_row(nu, j);
    double phi = gsl_vector_get(parent_phi, j);
    double nu_val = gsl_matrix_get(nu, j, topic);
    double nu_j_sum = gsl_blas_dsum(&nu_j.vector);
    new_phi +=  phi * (gsl_sf_psi(nu_val) - gsl_sf_psi(nu_j_sum));
    BOOST_ASSERT(!isnan(new_phi));
  }
  return new_phi;
}

double DocumentMapper::ComputeNewPhiRoot(const gsl_matrix* nu,
                                         const int topic,
                                         const int number_topics) {
  gsl_vector_const_view nu_j = gsl_matrix_const_row(nu, number_topics);
  double nu_val = gsl_matrix_get(nu, number_topics, topic);
  double nu_j_sum = gsl_blas_dsum(&nu_j.vector);
  double new_phi =  (gsl_sf_psi(nu_val) - gsl_sf_psi(nu_j_sum));
  BOOST_ASSERT(!isnan(new_phi));
  return new_phi;
}


double DocumentMapper::ComputeNewPhiChildren(Document* d,
                                             LensNode* w,
                                             const gsl_matrix* nu,
                                             const gsl_vector* gamma,
                                             const double gamma_normalizer,
                                             const int topic,
                                             const int number_topics) {
  double new_phi = 0.0;

  double nu_i_sum = gsl_blas_dsum
    (&gsl_matrix_const_row(nu, topic).vector);

  for (int m = 0; m < w->number_children_; m++) {
    int child = w->children_positions_[m];
    double child_contribution = 0;
    for (int j = 0; j < number_topics; j++) {
      double nu_val = gsl_matrix_get(nu, topic, j);
      child_contribution += (gsl_vector_get(gamma, j)
                             / gamma_normalizer) * (nu_val / nu_i_sum);
      new_phi += gsl_vector_get(d->get_word_(child)->phi_.get(), j) *
        (gsl_sf_psi(nu_val) - gsl_sf_psi(nu_i_sum));
    }
    /*
    cout << "    3 new_phi-=" << child_contribution << "/"
         << gsl_vector_get(d->omega_.get(), child) << endl;
    */
    new_phi -= child_contribution * gsl_vector_get(d->get_omega_(), child);
    BOOST_ASSERT(!isnan(new_phi));
  }

  return new_phi;
}

/*
 * Returns log of what phi should be
 */
double DocumentMapper::ComputeNewPhi(Document* d,
                                     const int word_index,
                                     const int topic,
                                     const gsl_matrix* nu,
                                     const gsl_matrix* tau,
                                     const gsl_vector* gamma,
                                     int number_terms,
                                     int number_topics) {
  LensNode* w = d->get_word_(word_index);
  int parent = w->get_parent_index_();
  gsl_vector* parent_phi = NULL;
  if (parent != -1) {
    parent_phi = d->get_word_(parent)->phi_.get();
  }

  // if (parent_phi != NULL) {
  // display_vector(parent_phi, "parent phi is\n");
  // }

  double gamma_normalizer = gsl_blas_dsum(gamma);

  double new_phi = ComputeNewPhiIntrinsic(topic,
                                          gamma,
                                          gamma_normalizer);

  // cout << "Setting phi " << new_phi;
  if (parent != -1) {
    // cout << "(P) ";
    new_phi += ComputeNewPhiParent(parent_phi,
                                   nu,
                                   topic,
                                   number_topics);
    BOOST_ASSERT(!isnan(new_phi));
  } else {
    new_phi += ComputeNewPhiRoot(nu, topic, number_topics);

    BOOST_ASSERT(check_valid_double(new_phi, "word " + SimpleItoa(word_index) +
                                    " topic " + SimpleItoa(topic) + "'s phi"));
  }

  // cout << " " << new_phi;
  double children = ComputeNewPhiChildren(d,
                                          w,
                                          nu,
                                          gamma,
                                          gamma_normalizer,
                                          topic,
                                          number_topics);
  // cout << " CHILDREN=" << children;
  new_phi += children;

  BOOST_ASSERT(check_valid_double(new_phi, "word " + SimpleItoa(word_index)
                            + " topic " + SimpleItoa(topic) + "'s phi"));

  double voc = gsl_matrix_get(tau, topic, w->term_);
  // cout << " VOC=" << voc;
  new_phi += voc;

  BOOST_ASSERT(check_valid_double(new_phi, "word " + SimpleItoa(word_index)
                                  + " topic " + SimpleItoa(topic) + "'s phi"));

  // cout << " NP " << new_phi << endl;

  return new_phi;
}

DocumentMapper::~DocumentMapper() {
}
