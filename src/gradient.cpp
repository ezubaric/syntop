/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#include "gradient.h"

DEFINE_int(max_conj_grad_iter,
           150,
           "Number of conjugate gradient iterations.");

DEFINE_double(conj_grad_tol,
              1e-2,
              "Sensitivity to convergence of conjugate gradient.");

DEFINE_double(conj_grad_step,
              0.01,
              "Conjugate gradient step size");

DEFINE_int(optimization_iteration_output,
           10,
           "How often we print status.");

DEFINE_double(exp_max,
              1e10,
              "How large exp can be (when there are positivity constraints)");

DEFINE_bool(write_grid,
            false,
            "Do we write out the values during optimization?");

DEFINE_bool(uniform_start,
            // true,
            false,
            "Do we start from uniform (1/N)?");

DEFINE_bool(always_take_answer,
            false,
            "If we end up in a spot worse than our start, do we take it?");

using std::numeric_limits;

// Can this be optimized to handle scale = 1.0?
double var_dir_lhood(double scale,
                     const gsl_vector* prior,
                     const gsl_vector* dirichlet) {
  double sum_val = gsl_blas_dsum(dirichlet);
  double sum_prior;

  sum_prior = scale * gsl_blas_dsum(prior);
  double result = gsl_sf_lngamma(sum_prior);
  for (int i = 0; i < (int)dirichlet->size; i++) {
    double prior_val = scale * gsl_vector_get(prior, i);
    double val = gsl_vector_get(dirichlet, i);

    result -= gsl_sf_lngamma(prior_val);
    result += (prior_val - 1.0) * (gsl_sf_psi(val) -
                                   gsl_sf_psi(sum_val));
  }

  return result;
}

double var_dir_lhood(const gsl_vector* prior,
                     const gsl_vector* dirichlet) {
  return var_dir_lhood(1.0, prior, dirichlet);
}

/*
third term in log-likelihood
 */
double gamma_nu_total_interaction(const gsl_matrix* exp_doc_counts,
                                  const int topic,
                                  const gsl_vector* nu) {
  double result = 0.0;
  double nu_sum = gsl_blas_dsum(nu);

  for (unsigned int j = 0; j < exp_doc_counts->size2; j++) {
    result += gsl_matrix_get(exp_doc_counts, topic, j) *
      gsl_vector_get(nu, j) / nu_sum;
  }
  return result;
}

/*
 * We can ignore the terms that are constant in terms of gamma an nu
 * if we want (useful for checking to see if an optimization his
 * improved the likelihood).
 */
double gamma_nu_interaction(Document* doc,
                            const gsl_vector* gamma,
                            const gsl_matrix* nu,
                            bool include_constants) {
  double result = 0;

  int num_words = doc->get_number_words_();
  double gamma_sum = gsl_blas_dsum(gamma);

  for (unsigned int i = 0; i < gamma->size; i++) {
    double topic_contribution = 0.0;
    for (int n = 0; n < num_words; n++) {
      double word_contribution = 0.0;
      LensNode * w = doc->get_word_(n);
      int parent = w->get_parent_index_();

      if (parent >= 0) {
        for (unsigned int j = 0; j < gamma->size; j++) {
          gsl_vector_const_view nu_j = gsl_matrix_const_row(nu, j);
          double nu_sum = gsl_blas_dsum(&nu_j.vector);

          double val =
            gsl_vector_get(doc->get_word_(parent)->phi_.get(), j) *
            (gsl_vector_get(&nu_j.vector, i) / nu_sum);
          word_contribution += val;
        }
      } else {
        // display_matrix(nu, "nu matrix is\n");
        // gamma->size is the number of topics, and this is the topic
        // specially for root words
        gsl_vector_const_view nu_j = gsl_matrix_const_row(nu, gamma->size);
        // display_vector(&nu_j.vector, "nu vector is\n");
        double nu_sum = gsl_blas_dsum(&nu_j.vector);
        word_contribution = (gsl_vector_get(&nu_j.vector, i) / nu_sum);
      }

      // We ignore omega when it's set to zero, as that means we don't
      // have to worry about normalization
      double omega = doc->get_omega_(n);
      // WARNING: the way we store omega is actually its inverse
      word_contribution *= omega;

      topic_contribution += word_contribution;
    }

    result += topic_contribution * gsl_vector_get(gamma, i) / gamma_sum;

    // cout << "result in gamma_nu_interaction is " << result << endl;
  }

  if (include_constants) {
    for (int n = 0; n < num_words; n++) {
      double omega = doc->get_omega_(n);
      result += -safe_log(omega) - 1.0;
    }
  }

  return result;
}

// Computes document likelihood for \nu_i for just a single document.
double nu_doc_transition_likelihood(Document* d,
                                    gsl_vector* nu,
                                    unsigned int i) {
  int num_words = d->get_number_words_();
  double val_sum = gsl_blas_dsum(nu);
  double result = 0;
  for (unsigned int j = 0; j < nu->size; j++) {
    double val = gsl_vector_get(nu, j);
    double counts = 0.0;
    // cout << "(";
    for (int n = 0; n < num_words; n++) {
      LensNode* w = d->get_word_(n);
      int parent = w->get_parent_index_();
      if (parent != -1 && i != nu->size) {
        counts += gsl_vector_get(w->phi_.get(), j) *
          gsl_vector_get(d->get_word_(parent)->phi_.get(), i);
        // cout << gsl_vector_get(w->phi_.get(), j) << "*" <<
        // gsl_vector_get(d->get_word_(parent)->phi_.get(), i) << "+";
      } else if (parent == -1 && i == nu->size) {
        counts += gsl_vector_get(w->phi_.get(), j);
      }
    }
    // cout << ")*" << (gsl_sf_psi(val) - gsl_sf_psi(val_sum)) << endl;
    // cout << "COUNTS:" << counts << endl;

    // cout << val << "\t" << val_sum << "\t" << gsl_sf_psi(val) << "\t"
    // << gsl_sf_psi(val_sum) << endl;

    result += counts * (gsl_sf_psi(val) - gsl_sf_psi(val_sum));

    // cout << "result inside nu_doc_trans_lhood is " << result << endl;
  }
  return result;
}

// Gives the likelihood term for transitions using \nu_i, thus all
// transitions that start with their assigned topic = i.  Note that
// this covers all transitions, and duplicates what's computed in
// nu_doc_transition_likelihood.
double nu_transition_likelihood(const gsl_matrix* counts,
                                const gsl_vector* nu,
                                int i) {
  double val_sum = gsl_blas_dsum(nu);
  double result = 0.0;
  for (unsigned int j = 0; j < nu->size; j++) {
    double val = gsl_vector_get(nu, j);
    /*
    cout << "trans_term += " << gsl_matrix_get(counts, i, j)
         << "*(gsl_sf_psi(" << val << ") - gsl_sf_psi(" << val_sum
         << ")" << endl;
    */
    result += gsl_matrix_get(counts, i, j) * (gsl_sf_psi(val)
                                              - gsl_sf_psi(val_sum));
  }
  return result;
}

double gamma_lhood_exp_counts(Document* d, const gsl_vector* gamma) {
  double result = 0.0;
  double val_sum = gsl_blas_dsum(gamma);
  int num_words = d->get_number_words_();
  for (int n = 0; n < num_words; n++) {
    LensNode * w = d->get_word_(n);
    for (unsigned int i = 0; i < gamma->size; i++) {
      double val = gsl_vector_get(gamma, i);
      result += gsl_vector_get(w->phi_.get(), i) *
        (gsl_sf_psi(val) - gsl_sf_psi(val_sum));
    }
  }

  return result;
}

double gamma_gradient_interaction(double gamma_val,
                                  double gamma_sum,
                                  const gsl_vector* gamma,
                                  const gsl_vector* nu_interaction,
                                  unsigned int component) {
  double off_terms = 0.0;
  for (unsigned int k = 0; k < nu_interaction->size; k++) {
    if (k != component) {
      double contribution = gsl_vector_get(gamma, k) *
        gsl_vector_get(nu_interaction, k);
      off_terms += contribution;
    }
  }

  double primary_terms = gsl_vector_get(nu_interaction, component) *
    (gamma_sum - gamma_val);

  return (primary_terms - off_terms) / (gamma_sum * gamma_sum);
}


/*
 * topic - the topic corresponding to the nu vector we're considering
 * gradient_component - the component we're taking the gradient of
 */
double nu_gradient_interaction(double nu_val,
                               double nu_sum,
                               const gsl_vector* nu,
                               gsl_matrix* nu_gamma_interaction,
                               unsigned int topic,
                               unsigned int gradient_component) {
  // These are the components we're not taking the gradient of and how
  // they affect the gadient
  double off_terms = 0.0;
  for (unsigned int k = 0; k < nu_gamma_interaction->size2; k++) {
    if (k != gradient_component) {
      off_terms += gsl_vector_get(nu, k) *
        gsl_matrix_get(nu_gamma_interaction, topic, k);
    }
  }

  double primary_terms = (gsl_matrix_get(nu_gamma_interaction,
                                         topic,
                                         gradient_component) *
                          (nu_sum - nu_val));

  double result = (primary_terms - off_terms) / (nu_sum * nu_sum);
  return result;
}

double nu_double_gradient_second_term(nu_params* p,
                                      const gsl_vector* nu,
                                      double nu_sum) {
  double second_term = 0.0;
  gsl_vector_view counts_row = gsl_matrix_row(p->transition_counts, p->topic);
  second_term += gsl_blas_dsum(&counts_row.vector);
  second_term += gsl_blas_dsum(p->alpha);
  second_term -= nu_sum;
  second_term *= -gsl_sf_psi_n(2, nu_sum);
  second_term += gsl_sf_psi_n(1, nu_sum);
  return second_term;
}

double nu_double_gradient_j(nu_params* p,
                            const gsl_vector* nu,
                            unsigned int j,
                            double nu_sum) {
  double nu_val = gsl_vector_get(nu, j);
  double second_term = nu_double_gradient_second_term(p, nu, nu_sum);

  double first_term = gsl_sf_psi_n(2, nu_val)
    * (gsl_vector_get(p->alpha, j) +
       gsl_matrix_get(p->transition_counts, p->topic, j) -
       nu_val);
  first_term -= gsl_sf_psi_n(1, nu_val);

  // ke zhai: in reducer, the nu_gamma_interaction is missing...
  double third_term = 0.0;
  for (unsigned int d = 0; d < p->nu_gamma_interaction->size1; d++) {
    double document_term = 0.0;
    for (unsigned int k = 0; k < p->nu_gamma_interaction->size2; k++) {
      if (k != j) {
        document_term += gsl_matrix_get(p->nu_gamma_interaction, d, k) *
          gsl_vector_get(nu, k);
      }
    }

    document_term = (document_term - gsl_matrix_get(p->nu_gamma_interaction,
                                                    d,
                                                    j)) /
      pow(nu_sum, 3.0);
    third_term += document_term;
  }

  /*
  double primary_terms = (gsl_matrix_get(p->nu_gamma_interaction,
                                         p->topic,
                                         j) *
                          (nu_sum - nu_val));

  double third_term = -2.0 * primary_terms / (nu_sum * nu_sum * nu_sum);
  */

  double grad_j = first_term + second_term - third_term;
  return grad_j;
}

double nu_gradient_second_term(nu_params* p,
                               const gsl_vector* nu,
                               double nu_sum) {
  double second_term = 0.0;
  gsl_vector_view counts_row = gsl_matrix_row(p->transition_counts, p->topic);
  second_term += gsl_blas_dsum(&counts_row.vector);
  second_term += gsl_blas_dsum(p->alpha);
  second_term -= nu_sum;
  second_term *= gsl_sf_psi_n(1, nu_sum);
  return second_term;
}

double nu_gradient_j(nu_params* p,
                     const gsl_vector* nu,
                     int j,
                     double second_term,
                     double nu_sum) {
  double nu_val = gsl_vector_get(nu, j);
  BOOST_ASSERT(nu_val > 0);
  BOOST_ASSERT((nu->size + 1) == p->nu_gamma_interaction->size1);
  BOOST_ASSERT(p->nu_gamma_interaction->size1 == p->transition_counts->size1);
  BOOST_ASSERT(nu->size == p->nu_gamma_interaction->size2);
  BOOST_ASSERT(p->nu_gamma_interaction->size2 == p->transition_counts->size2);

  /*
  if(nu_val <= 0) {
    return numeric_limits<double>::infinity();
    }*/
  double first_term = gsl_vector_get(p->alpha, j) +
    gsl_matrix_get(p->transition_counts, p->topic, j) -
    nu_val;
  first_term *= gsl_sf_psi_n(1, nu_val);
  double third_term = nu_gradient_interaction(nu_val,
                                              nu_sum,
                                              nu,
                                              p->nu_gamma_interaction,
                                              p->topic,  // maximizing this nu
                                              j);        // this grad component

  double grad_j;
  grad_j = first_term - second_term - third_term;

  return grad_j;
}

double nu_gradient_j(nu_params* p, const gsl_vector* nu, int j, double nu_sum) {
  double second_term = nu_gradient_second_term(p, nu, nu_sum);
  return nu_gradient_j(p, nu, j, second_term, nu_sum);
}

void log_nu_gradient(const gsl_vector* x,
                     void* params,
                     gsl_vector* gradient) {
  nu_params* p = (nu_params*) params;
  vexp(x, p->nu_exp, FLAGS_exp_max);

  double nu_sum = gsl_blas_dsum(p->nu_exp);

  // The second term is shared by all components
  double second_term = nu_gradient_second_term(p, p->nu_exp, nu_sum);

  for (unsigned int j = 0; j < p->alpha->size; j++) {
    // double grad_j = first_term;
    // We multply by exp(x) because of the implicit differentiation
    double grad_j = nu_gradient_j(p, p->nu_exp, j, second_term, nu_sum) *
      exp(gsl_vector_get(x, j));
    gsl_vector_set(gradient, j, -1.0 * grad_j);
  }
}

void nu_gradient(const gsl_vector* nu,
                 void* params,
                 gsl_vector* gradient) {
  nu_params* p = (nu_params*) params;
  double nu_sum = gsl_blas_dsum(nu);

  // The second term is shared by all components
  double second_term = nu_gradient_second_term(p, nu, nu_sum);

  for (unsigned int j = 0; j < p->alpha->size; j++) {
    // double grad_j = first_term;
    double grad_j = nu_gradient_j(p, nu, j, second_term, nu_sum);
    gsl_vector_set(gradient, j, -1.0 * grad_j);
  }
}

double beta_prior_lhood(gsl_vector* T, double alpha) {
  int trunc = T->size - 1;
  double val = (alpha - 1.0) * safe_log(gsl_vector_get(T, trunc));
  for (int i = 0; i < trunc; i++) {
    val -= safe_log(gsl_vector_get(T, i));
  }

  return val;
}

double beta_trans_lhood(gsl_vector* beta,
                        gsl_vector* doc_terms,
                        gsl_vector* trans_terms,
                        int num_docs,
                        bool ignore_trans,
                        double alpha_doc,
                        double alpha_trans) {
  // We can't ignore both features
  BOOST_ASSERT(num_docs > 0 || !ignore_trans);

  double doc_val = 0.0;
  double trans_val = 0.0;
  int num_topics = beta->size;

  // Beta sums to 1.0
  doc_val += (double)num_docs * gsl_sf_lngamma(alpha_doc);
  trans_val += (double)(num_topics + 1) * gsl_sf_lngamma(alpha_trans);
  for (int i = 0; i < num_topics; i++) {
    double beta_val = gsl_vector_get(beta, i);
    trans_val -= (double)(num_topics + 1) *
      gsl_sf_lngamma(beta_val * alpha_trans);
    doc_val -= (double)num_docs *
      gsl_sf_lngamma(beta_val * alpha_doc);
    doc_val += (alpha_doc * beta_val - 1.0) *
      gsl_vector_get(doc_terms, i);
    trans_val += (alpha_trans * beta_val - 1.0) *
      gsl_vector_get(trans_terms, i);
  }

  if (num_docs == 0) {
    return trans_val;
  } else if (ignore_trans) {
    return doc_val;
  } else {
    return doc_val + trans_val;
  }
}

/*
 * This is the part of lhood bound that comes from both beta and gamma
 * as well as nu
 */
double beta_trans_lhood(gsl_vector* beta,
                        gsl_matrix* gamma,
                        gsl_matrix* nu,
                        double alpha_doc,
                        double alpha_trans) {
  double val = 0.0;
  for (unsigned int d = 0; d < gamma->size1; d++) {
    gsl_vector_view gamma_d = gsl_matrix_row(gamma, d);
    val += var_dir_lhood(alpha_doc, beta, &gamma_d.vector);
  }
  for (unsigned int k = 0; k < nu->size1; k++) {
    gsl_vector_view nu_d = gsl_matrix_row(nu, k);
    val += var_dir_lhood(alpha_trans, beta, &nu_d.vector);
  }
  return val;
}

/*
 * This is the part of the lhood bound that just comes from beta and gamma
 */
double beta_doc_lhood(gsl_vector* beta,
                      gsl_matrix* gamma,
                      double alpha_doc) {
  double val = 0.0;
  for (unsigned int d = 0; d < gamma->size1; d++) {
    gsl_vector_view gamma_d = gsl_matrix_row(gamma, d);
    val += var_dir_lhood(alpha_doc, beta, &gamma_d.vector);
  }
  // cout << "BETA GAMMA: " << val << endl;
  return val;
}

void beta_gradient(const gsl_vector* source,
                   const gsl_vector* T,
                   const gsl_vector* beta,
                   const gsl_vector* doc_terms,
                   const gsl_vector* trans_terms,
                   unsigned int num_docs,
                   unsigned int num_topics,
                   bool ignore_trans,
                   double alpha_top,
                   double alpha_doc,
                   double alpha_trans,
                   gsl_vector* grad) {
  unsigned int trunc = grad->size;
  BOOST_ASSERT(trunc == T->size - 1);
  double second_term = (alpha_top - 1.0) / gsl_vector_get(T, trunc);

  double first_term = 0.0;
  // display_vector(T, "T vector in beta_gradient is\n");

  double beta_K = gsl_vector_get(beta, trunc);
  // display_vector(beta, "beta vector in beta_gradient is\n");

  /*
   * The last position of the beta vector is implicitly defined in
   * terms of all of the other previous terms.  Thus, the final
   * position's setting effects all aother positions
   */
  double last_element_terms = 0;
  if (num_docs > 0) {
    last_element_terms += (double)num_docs *
      alpha_doc * gsl_sf_psi(alpha_doc * beta_K);
    last_element_terms -= alpha_doc * gsl_vector_get(doc_terms, trunc);
  }
  if (!ignore_trans) {
    last_element_terms += (double)(num_topics + 1) *
      alpha_trans * gsl_sf_psi(alpha_trans * beta_K);
    last_element_terms -= alpha_trans * gsl_vector_get(trans_terms, trunc);
  }

  // This is the normalizing constant; we add an additional zero
  // because of the implicit last term.

  display_vector(source, "source in beta_gradient is\n");
  double source_sum = log_sum(log_sum(source), 0.0);

  for (int k = trunc - 1; k >= 0; k--) {
    double beta_val = gsl_vector_get(beta, k);
    double third_term = 0;
    if (num_docs > 0) {
      third_term += (double)num_docs *
        -(alpha_doc * gsl_sf_psi(alpha_doc * beta_val));
      third_term += alpha_doc * gsl_vector_get(doc_terms, k);
    }

    double fourth_term = 0;
    if (!ignore_trans) {
      fourth_term += (double)(num_topics + 1) *
        -(alpha_trans * gsl_sf_psi(alpha_trans * beta_val));
      fourth_term += alpha_trans * gsl_vector_get(trans_terms, k);
    }

    double implicit_term = log_diff(source_sum + gsl_vector_get(source, k),
                                    2.0 * gsl_vector_get(source, k)) -
      2.0 * source_sum;

    /*
    cout << "BETA GRAD(" << k << "): " << first_term << " " << second_term
         << " " << third_term << " " << fourth_term << " "
         << last_element_terms << " " << implicit_term << endl;
    */

    gsl_vector_set(grad,
                   k,
                   - exp(implicit_term) * (first_term
                                           - second_term
                                           + third_term
                                           + fourth_term
                                           + last_element_terms));

    // We do this after setting the kth position because the sum is
    // from k+1 to K-1
    first_term += 1.0 / gsl_vector_get(T, k);
  }
}

void logspace_beta(const gsl_vector* full_beta,
                   gsl_vector* logspace_beta) {
  gsl_vector_const_view trunc_beta_view =
    gsl_vector_const_subvector(full_beta,
                               0,
                               full_beta->size - 1);
  const gsl_vector* trunc_beta = &trunc_beta_view.vector;

  gsl_vector_memcpy(logspace_beta, trunc_beta);

  // Make it so the last coordinate is 1.0 (if we weren't truncated)
  gsl_vector_scale(logspace_beta,
                   1.0 / gsl_vector_get(full_beta, full_beta->size - 1));

  // Now transform it into logspace.
  vlog(logspace_beta, logspace_beta);
}

/*
 * Creates T and the full beta from the truncated beta representation.
 * Because sometimes we only need T, if full_beta is NULL we don't try
 * to save the full_beta.
 */
void create_T_full_beta(const gsl_vector* trunc_beta_source,
                        const unsigned int num_topics,
                        gsl_vector* full_beta,
                        gsl_vector* T,
                        bool log_space) {
  double running_sum = 0.0;

  BOOST_ASSERT(trunc_beta_source->size == T->size ||
               trunc_beta_source->size == T->size - 1);
  BOOST_ASSERT(T->size == num_topics);

  BOOST_ASSERT((log_space && full_beta) || !log_space);

  gsl_vector_set_all(T, 0.0);

  // Because we don't know the length of the incomming beta, this normalizes it
  gsl_vector_const_view trunc_beta_view =
    gsl_vector_const_subvector(trunc_beta_source,
                               0,
                               num_topics - 1);
  const gsl_vector* trunc_beta = &trunc_beta_view.vector;
  gsl_vector_view truncated_T = gsl_vector_subvector(T,
                                                     0,
                                                     num_topics - 1);

  if (log_space) {
    double normalizer = log_sum(0.0, log_sum(trunc_beta_source));
    gsl_vector_memcpy(&truncated_T.vector, trunc_beta);
    // We subtract out from the whole vector because the last position
    // is 0.0, which defines the Kth component
    gsl_vector_add_constant(T, -normalizer);

    // Now we can exponentiate
    vexp(T, T, FLAGS_exp_max);
  } else {
    gsl_vector_memcpy(&truncated_T.vector, trunc_beta);
  }

  if (full_beta) {
    BOOST_ASSERT(full_beta->size == T->size);
  }

  for (unsigned int i = 0; i < num_topics-1; i++) {
    double beta_val = gsl_vector_get(T, i);

    gsl_vector_set(T, i, 1.0 - running_sum);
    if (full_beta) {
      gsl_vector_set(full_beta, i, beta_val);
    }
    running_sum += beta_val;
    // cout << "running sum is " << running_sum << endl;

    assert(running_sum <= (1.0));
  }

  gsl_vector_set(T, num_topics-1, 1.0 - running_sum);
  if (full_beta) {
    gsl_vector_set(full_beta, num_topics-1, 1.0 - running_sum);
    // display_vector(full_beta, "beta");
  }
}

double beta_lhood(const gsl_vector* beta, void* params) {
  beta_params* p = (beta_params*) params;
  // make sure beta strickly follows a stick breaking process
  create_T_full_beta(beta, beta->size + 1, p->full_beta, p->T, true);

  double prior = beta_prior_lhood(p->T, p->alpha_top);
  double trans = beta_trans_lhood(p->full_beta,
                                  p->doc_terms,
                                  p->trans_terms,
                                  p->num_docs,
                                  p->ignore_trans,
                                  p->alpha_doc,
                                  p->alpha_trans);

  // -prior? i think it should be +prior?
  double val =  - (prior + trans);

  return val;
}

void beta_lhood_gradient(const gsl_vector* beta,
                         void* params,
                         double* lhood,
                         gsl_vector* grad) {
  beta_params* p = (beta_params*) params;
  create_T_full_beta(beta, beta->size + 1, p->full_beta, p->T, true);

  *lhood = -(beta_prior_lhood(p->T, p->alpha_top) +
             beta_trans_lhood(p->full_beta,
                              p->doc_terms,
                              p->trans_terms,
                              p->num_docs,
                              p->ignore_trans,
                              p->alpha_doc,
                              p->alpha_trans));

  beta_gradient(beta,
                p->T,
                p->full_beta,
                p->doc_terms,
                p->trans_terms,
                p->num_docs,
                p->ignore_trans,
                p->full_beta->size,
                p->alpha_top,
                p->alpha_doc,
                p->alpha_trans,
                grad);
}

void beta_gradient_wrapper(const gsl_vector* beta,
                           void* params,
                           gsl_vector* grad) {
  beta_params* p = (beta_params*) params;
  create_T_full_beta(beta, beta->size + 1, p->full_beta, p->T, true);

  // display_vector(beta, "beta is\n");

  beta_gradient(beta,
                p->T,
                p->full_beta,
                p->doc_terms,
                p->trans_terms,
                p->num_docs,
                p->full_beta->size,
                p->ignore_trans,
                p->alpha_top,
                p->alpha_doc,
                p->alpha_trans,
                grad);
}

double nu_token_lhood(const gsl_vector* nu, void* params) {
  double result = 0;

  nu_params * p = (nu_params*) params;

  double first_term = var_dir_lhood(p->alpha, nu);

  // double second_term = nu_transition_likelihood(p->transition_counts,
  //                                               nu,
  //                                               p->topic);

  // We're ignoring the log(omega_n) - 1 term, but it's constant wrt nu
  // double third_term = gamma_nu_total_interaction(p->nu_gamma_interaction,
  //                                                p->topic,
  //                                                nu);

  double fourth_term = var_dir_lhood(nu, nu);

  // result = first_term + second_term - third_term - fourth_term;
  result = first_term - fourth_term;

  // Because we're maximizing, we negate
  // KE ZHAI: making sure the resulted nu follows a dirichlet distribution
  BOOST_ASSERT(gsl_vector_min(nu) > 0);
  if (gsl_vector_min(nu) <= 0) {
    cout << "BADNESS!" << endl;
    return numeric_limits<double>::infinity();
  }

  // KE ZHAI: why negate the likelihood?
  // return -1.0 * result;
  return result;
}

double nu_lhood(const gsl_vector* nu, void* params) {
  double result = 0;

  nu_params * p = (nu_params*) params;

  double first_term = var_dir_lhood(p->alpha, nu);

  double second_term = nu_transition_likelihood(p->transition_counts,
                                                nu,
                                                p->topic);

  // We're ignoring the log(omega_n) - 1 term, but it's constant wrt nu
  double third_term = gamma_nu_total_interaction(p->nu_gamma_interaction,
                                                 p->topic,
                                                 nu);

  double fourth_term = var_dir_lhood(nu, nu);

  result = first_term + second_term - third_term - fourth_term;

  // cout << "terms for nu likelihood are " << first_term << "\t"
  // << second_term << "\t" << third_term << "\t" << fourth_term << endl;

  // Because we're maximizing, we negate
  BOOST_ASSERT(gsl_vector_min(nu) > 0);
  if (gsl_vector_min(nu) <= 0) {
    cout << "BADNESS!" << endl;
    return numeric_limits<double>::infinity();
  }

  return -1.0 * result;
}

double log_nu_lhood(const gsl_vector* x, void* params) {
  nu_params * p = (nu_params*) params;
  vexp(x, p->nu_exp, FLAGS_exp_max);
  return nu_lhood(p->nu_exp, params);
}

void gamma_gradient_base(const gsl_vector* gamma,
                         const gsl_vector* alpha,
                         const gsl_vector* phi_sum,
                         const gsl_vector* nu_interaction,
                         gsl_vector* gradient) {
  double second_term = 0.0;
  double val_sum = gsl_blas_dsum(gamma);

  second_term = gsl_blas_dsum(alpha);
  second_term += gsl_blas_dsum(phi_sum);
  second_term -= gsl_blas_dsum(gamma);
  second_term *= gsl_sf_psi_1(val_sum);

  for (unsigned int i = 0; i < gamma->size; i++) {
    double val = gsl_vector_get(gamma, i);

    double first_term = gsl_sf_psi_1(val) * (gsl_vector_get(alpha, i) +
                                         gsl_vector_get(phi_sum, i) - val);

    double third_term = gamma_gradient_interaction(val,
                                                   val_sum,
                                                   gamma,
                                                   nu_interaction,
                                                   i);

    double gradient_i = first_term - second_term - third_term;

    // Because we're actually maximizing, we negate; because of the
    // implicit differentiation, we multiply by exp of the variable
    // we're maximizing, which happens to be "gamma"
    gsl_vector_set(gradient, i, -gsl_vector_get(gamma, i) * gradient_i);
  }
}

double gamma_lhood_base(const gsl_vector* gamma,
                        const gsl_vector* alpha,
                        const gsl_vector* phi_sum,
                        const gsl_vector* nu_interaction) {
  double result = 0;
  double gamma_sum = gsl_blas_dsum(gamma);

  if (!check_valid_double(gamma_sum, "gamma sum")) {
    display_vector(gamma, "gamma fault");
    display_vector(alpha, "alpha fault");
    display_vector(phi_sum, "phi sum fault");
    display_vector(nu_interaction, "nu int fault");
  }

  // double first_term = var_dir_lhood(alpha, gamma);
  // KE ZHAI: should this be beta, instead of gamma?
  double first_term = var_dir_lhood(alpha, gamma);

  double second_term = 0.0;
  for (unsigned i = 0; i < phi_sum->size; i++) {
    second_term += gsl_vector_get(phi_sum, i) *
      (gsl_sf_psi(gsl_vector_get(gamma, i)) - gsl_sf_psi(gamma_sum));
  }
  // We want to use the full nu matrix, so final two arguments are
  // ignored.
  double third_term = dot_product(nu_interaction, gamma) / gamma_sum;
  double fourth_term = var_dir_lhood(gamma, gamma);

  result = first_term + second_term - third_term - fourth_term;

  // Because we're actually maximizing, we negate
  return -1.0 * result;
}

double gamma_lhood(const gsl_vector* gamma, void* params) {
  gamma_params* p = (gamma_params*)params;
  vexp(gamma, p->gamma_exp, FLAGS_exp_max);
  return gamma_lhood_base(p->gamma_exp,
                          p->alpha,
                          p->phi_sum,
                          p->nu_interaction);
}

void gamma_lhood_gradient(const gsl_vector* gamma,
                          void* params,
                          double* lhood,
                          gsl_vector* grad) {
  gamma_params* p = (gamma_params*) params;
  vexp(gamma, p->gamma_exp, FLAGS_exp_max);
  *lhood = gamma_lhood_base(p->gamma_exp,
                            p->alpha,
                            p->phi_sum,
                            p->nu_interaction);
  gamma_gradient_base(p->gamma_exp,
                      p->alpha,
                      p->phi_sum,
                      p->nu_interaction,
                      grad);
}

void gamma_gradient(const gsl_vector* gamma,
                    void* params,
                    gsl_vector* gradient) {
  gamma_params* p = (gamma_params*) params;
  vexp(gamma, p->gamma_exp, FLAGS_exp_max);
  gamma_gradient_base(p->gamma_exp,
                      p->alpha,
                      p->phi_sum,
                      p->nu_interaction,
                      gradient);
}


void log_nu_lhood_gradient(const gsl_vector* x,
                       void* params,
                       double* lhood,
                       gsl_vector* grad) {
  nu_params* p = (nu_params*) params;
  vexp(x, p->nu_exp, FLAGS_exp_max);

  double nu_sum = gsl_blas_dsum(p->nu_exp);

  // The second term is shared by all components
  double second_term = nu_gradient_second_term(p, p->nu_exp, nu_sum);

  // cout << "inside log nu lhood gradient..." << endl;

  for (unsigned int j = 0; j < p->alpha->size; j++) {
    // double grad_j = first_term;
    // We multply by exp(x) because of the implicit differentiation
    double grad_j = nu_gradient_j(p, p->nu_exp, j, second_term, nu_sum) *
      exp(gsl_vector_get(x, j));
    gsl_vector_set(grad, j, -1.0 * grad_j);
  }

  // cout << "about to leave log nu lhood gradient..." << endl;

  *lhood = nu_lhood(p->nu_exp, params);
}


void nu_lhood_gradient(const gsl_vector* nu,
                       void* params,
                       double* lhood,
                       gsl_vector* grad) {
  *lhood = nu_lhood(nu, params);
  // cout << "Computed lhood = " << *lhood << " for nu update." << endl;
  nu_gradient(nu, params, grad);
}

void beta_const_terms(const gsl_matrix* gamma,
                      const gsl_matrix* nu,
                      const int num_topics,
                      gsl_vector* doc_terms,
                      gsl_vector* trans_terms) {
  for (int k = 0; k < num_topics; k++) {
    double const_term = 0.0;
    for (unsigned int d = 0; d < gamma->size1; d++) {
      gsl_vector_const_view gamma_row = gsl_matrix_const_row(gamma, d);
      const gsl_vector* gamma_d = &gamma_row.vector;
      const_term += (gsl_sf_psi(gsl_vector_get(gamma_d, k)) -
                     gsl_sf_psi(gsl_blas_dsum(gamma_d)));
    }
    gsl_vector_set(doc_terms, k, const_term);

    const_term = 0.0;
    for (unsigned int z = 0; z < nu->size1; z++) {
      gsl_vector_const_view nu_row = gsl_matrix_const_row(nu, z);
      const gsl_vector* nu_z = &nu_row.vector;
      const_term += (gsl_sf_psi(gsl_vector_get(nu_z, k)) -
                     gsl_sf_psi(gsl_blas_dsum(nu_z)));
    }
    gsl_vector_set(trans_terms, k, const_term);
  }

  // display_vector(trans_terms, "BETA  TRANS");
  // display_vector(doc_terms, "BETA DOC");
}


double maximize(gsl_vector* start_position,
                double (*lhood)(const gsl_vector*, void*),
                void (*grad)(const gsl_vector*, void*, gsl_vector*),
                void (*lhood_grad)(const gsl_vector*,
                                   void*,
                                   double*,
                                   gsl_vector*),
                void* params,
                string base_name,
                string name,
                const int size) {
  const gsl_multimin_fdfminimizer_type* T;
  gsl_multimin_fdfminimizer* s;

  gsl_multimin_function_fdf f;
  f.f = lhood;
  f.df = grad;
  f.fdf = lhood_grad;
  f.params = params;
  f.n = size;

  shared_ptr<gsl_vector> x(gsl_vector_alloc(start_position->size),
                           gsl_vector_free);
  gsl_vector_memcpy(x.get(), start_position);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc(T, size);

  gsl_multimin_fdfminimizer_set(s,
                                &f,
                                x.get(),
                                FLAGS_conj_grad_step,
                                FLAGS_conj_grad_tol);

  int status;
  int iter = 0;

  double start = lhood(start_position, params);

  do {
    /*
    if(iter % FLAGS_optimization_iteration_output == 0) {
      cout << iter << "(" << s->f << ")" << " " << endl;
    }
    */
    iter++;

    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;

    status = gsl_multimin_test_gradient(s->gradient, FLAGS_conj_grad_tol);
  } while (status == GSL_CONTINUE && iter < FLAGS_max_conj_grad_iter);

  if (status == GSL_SUCCESS) {
    cout << " SUC " << endl;
  } else if (status != GSL_CONTINUE) {
    switch (status) {
    case GSL_ENOPROG:
      cout << " NP " << endl;
    };
  }

  // display_vector(s->x, "End x");
  // display_vector(start_position, "final val");

  double final = lhood(s->x, params);
  if (final < start) {
    // cout << "Keeping new value for " << name <<
    //  " (start=" << start << ", final=" << final << ")" << endl;
    gsl_vector_memcpy(start_position, s->x);
  } else {
    cout << "Throwing out value, no increase for " << name <<
      " (start=" << start << ", final=" << final << ")" << endl;
  }

  // cout << "*** Overall " << name << " increase = " << start - final << " ***
  // " << endl;

  gsl_multimin_fdfminimizer_free(s);

  return start - final;
}

double optimize_beta(gsl_vector* beta,
                     gsl_matrix* gamma,
                     gsl_matrix* nu,
                     bool ignore_docs,
                     bool ignore_trans,
                     double alpha_top,
                     double alpha_doc,
                     double alpha_trans,
                     string model_name) {
  beta_params p;

  p.alpha_top = alpha_top;
  p.alpha_doc = alpha_doc;
  p.alpha_trans = alpha_trans;

  shared_ptr<gsl_vector> doc_terms(gsl_vector_alloc(beta->size),
                                   gsl_vector_free);
  display_vector(doc_terms.get(), "doc terms are\n");

  shared_ptr<gsl_vector> trans_terms(gsl_vector_alloc(beta->size),
                                     gsl_vector_free);
  display_vector(trans_terms.get(), "trans terms are\n");

  beta_const_terms(gamma,
                   nu,
                   beta->size,
                   doc_terms.get(),
                   trans_terms.get());
  p.doc_terms = doc_terms.get();
  p.trans_terms = trans_terms.get();
  if (ignore_docs) {
    p.num_docs = 0;
  } else {
    p.num_docs = gamma->size1;
  }

  p.ignore_trans = ignore_trans;

  shared_ptr<gsl_vector> full_beta(gsl_vector_alloc(beta->size),
                                   gsl_vector_free);
  shared_ptr<gsl_vector> tail(gsl_vector_alloc(beta->size),
                              gsl_vector_free);
  gsl_vector_set_all(full_beta.get(), 0);
  gsl_vector_set_all(tail.get(), 0);

  // display_vector(tail.get(), "tail in beta optimize\n");
  // display_vector(full_beta.get(), "full beta in beta optimize\n");

  p.T = tail.get();
  p.full_beta = full_beta.get();

  shared_ptr<gsl_vector> x(gsl_vector_alloc(beta->size - 1),
                           gsl_vector_free);
  shared_ptr<gsl_vector> T(gsl_vector_alloc(beta->size),
                           gsl_vector_free);
  gsl_vector_set_all(x.get(), 0);
  gsl_vector_set_all(T.get(), 0);

  logspace_beta(beta, x.get());

  // display_vector(T.get(), "T in beta optimize is\n");
  // display_vector(x.get(), "x in beta optimize is\n");
  // display_vector(beta, "beta in beta optimize is\n");

  double val = maximize(x.get(),
                        &beta_lhood,
                        &beta_gradient_wrapper,
                        &beta_lhood_gradient,
                        &p,
                        model_name,
                        "beta",
                        beta->size - 1);

  // Need to convert back into normal scale
  create_T_full_beta(x.get(),
                     beta->size,
                     beta,
                     T.get(),
                     true);

  return val;
}

/*
 * The second optional argument allows you to ignore all but one topic
 * for optimization
 */
double optimize_nu(VariationalParameters* vars,
                   const SyntopParameters* params,
                   int ignore_other_topics) {
  unsigned int num_topics = params->num_topics();
  double difference = 0;
  nu_params p;

  shared_ptr<gsl_vector> x(gsl_vector_alloc(num_topics), gsl_vector_free);
  shared_ptr<gsl_vector> nu_exp(gsl_vector_alloc(num_topics), gsl_vector_free);
  // display_vector(nu_exp.get(), "NU_EXP is\n");
  p.nu_exp = nu_exp.get();

  p.transition_counts = vars->nu_exp_topic_counts_.get();
  p.nu_gamma_interaction = vars->nu_exp_doc_counts_.get();

  gsl_vector* beta =  vars->beta_.get();
  shared_ptr<gsl_vector> alpha;
  alpha.reset(gsl_vector_alloc(beta->size), gsl_vector_free);
  p.alpha = alpha.get();
  // display_vector(alpha.get(), "alpha at checkpoint 1");

  gsl_vector_memcpy(p.alpha, beta);
  // display_vector(p.alpha, "alpha at checkpoint 3");
  // cout << "params->alpha_trans()\t" << params->alpha_trans() << endl;

  gsl_vector_scale(p.alpha, params->alpha_trans());

  // This is less than or equal because there is a special start state that
  // cannot be transitioned into for each sentence head
  for (unsigned int topic = 0; topic <= num_topics; topic++) {
    if (ignore_other_topics >= 0 &&
        (unsigned int) ignore_other_topics != topic) {
      continue;
    }

    // This scales self-transitions downward if the self-trans
    // parameter has been set (default is 1.0, which does nothing).
    double initial_self_alpha;
    if (topic < num_topics) {
      initial_self_alpha = gsl_vector_get(alpha.get(), topic);
      gsl_vector_set(alpha.get(), topic, initial_self_alpha);
    }

    gsl_vector_view nu_view = gsl_matrix_row(vars->nu_.get(), topic);
    gsl_vector* nu = &nu_view.vector;

    p.topic = topic;

    double start = nu_lhood(nu, &p);
    if (FLAGS_uniform_start) {
      gsl_vector_set_all(x.get(), params->alpha_trans() / (float)num_topics);
      vlog(x.get(), x.get());
    } else {
      vlog(nu, x.get());
    }

    maximize(x.get(),
             log_nu_lhood,
             log_nu_gradient,
             log_nu_lhood_gradient,
             &p,
             params->model_name(),
             "nu." + SimpleItoa(topic),
             num_topics);

    vexp(x.get(), nu, FLAGS_exp_max);

    difference += nu_lhood(nu, &p) - start;
    // cout << "likelihood difference inside optimize nu is " <<
    // "(negative is increasing) " << difference << endl;

    BOOST_ASSERT(gsl_vector_min(nu) > 0);

    if (topic < num_topics) {
      gsl_vector_set(alpha.get(), topic, initial_self_alpha);
    }
  }
  return -difference;
}

// Returns true if the value improved
double optimize_gamma(gsl_vector* gamma,
                      gsl_vector* phi_total,
                      gsl_vector* nu_interaction,
                      string model_name,
                      gsl_vector* beta,
                      double alpha_doc,
                      int doc_num) {
  gamma_params p;
  int num_topics = gamma->size;

  shared_ptr<gsl_vector> alpha(gsl_vector_alloc(beta->size), gsl_vector_free);
  shared_ptr<gsl_vector> x(gsl_vector_alloc(num_topics), gsl_vector_free);
  shared_ptr<gsl_vector> gamma_exp(gsl_vector_alloc(num_topics),
                                   gsl_vector_free);
  p.gamma_exp = gamma_exp.get();

  p.phi_sum = phi_total;
  p.nu_interaction = nu_interaction;
  p.alpha = alpha.get();
  gsl_vector_memcpy(p.alpha, beta);
  gsl_vector_scale(p.alpha, alpha_doc);

  vlog(gamma, x.get());
  double initial_value = gamma_lhood(x.get(), &p);
  if (FLAGS_uniform_start) {
    gsl_vector_set_all(x.get(), alpha_doc / (double)gamma->size);
    vlog(x.get(), x.get());
  }

  maximize(x.get(),
           gamma_lhood,
           gamma_gradient,
           gamma_lhood_gradient,
           &p,
           model_name,
           "gamma" + SimpleItoa(doc_num),
           num_topics);

  vexp(x.get(), gamma, FLAGS_exp_max);

  BOOST_ASSERT(gsl_vector_min(gamma) > 0);
  // cout << "gamma likelihood after maximization is " << -initial_value << "\t"
  // << -gamma_lhood(x.get(), &p) << "\t"
  // << -gamma_lhood(x.get(), &p)+initial_value << endl;
  return -(gamma_lhood(x.get(), &p) - initial_value);
}
