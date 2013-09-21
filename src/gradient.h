/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef GRADIENT_INCLUDED
#define GRADIENT_INCLUDED

#include <string>

#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"

#include "topicmod/lib/util/flags.h"
#include "topicmod/lib/util/strings.h"

#include "util.h"
#include "lens_doc.h"
#include "variational_parameters.h"

using lib_util::SimpleItoa;
using lib_prob::log_diff;
using lib_prob::safe_log;

#define psiexp(x) if (x > 20.0) x else gsl_sf_psi(exp(x)); // NOLINT

struct beta_params {
  gsl_vector* full_beta;
  gsl_vector* T;
  gsl_vector* doc_terms;
  gsl_vector* trans_terms;
  int num_docs;
  bool ignore_trans;
  double alpha_top;
  double alpha_doc;
  double alpha_trans;
};

struct nu_params {
  gsl_matrix* transition_counts;
  gsl_matrix* nu_gamma_interaction;
  gsl_vector* alpha;
  gsl_vector* nu_exp;
  unsigned int topic;
};

struct gamma_params {
  gsl_vector* gamma_exp;
  gsl_vector* alpha;
  gsl_vector* phi_sum;
  gsl_vector* nu_interaction;
};


// Can this be optimized to handle scale = 1.0?
double var_dir_lhood(double scale,
                     const gsl_vector* prior,
                     const gsl_vector* dirichlet);

double var_dir_lhood(const gsl_vector* prior,
                     const gsl_vector* dirichlet);

double gamma_nu_total_interaction(const gsl_matrix* exp_doc_counts,
                                  const int topic,
                                  const gsl_vector* nu);

/*
 * We can ignore the terms that are constant in terms of gamma an nu
 * if we want (useful for checking to see if an optimization his
 * improved the likelihood).
 */
double gamma_nu_interaction(Document* doc,
                            const gsl_vector* gamma,
                            const gsl_matrix* nu,
                            bool include_constants = true);

// Computes document likelihood for \nu_i for just a single document.
double nu_doc_transition_likelihood(Document* d,
                                    gsl_vector* nu,
                                    unsigned int i);

// Gives the likelihood term for transitions using \nu_i, thus all
// transitions that start with their assigned topic = i.  Note that
// this covers all transitions, and duplicates what's computed in
// nu_doc_transition_likelihood.
double nu_transition_likelihood(const gsl_matrix* counts,
                                const gsl_vector* nu,
                                int i);

double gamma_lhood_exp_counts(Document* d, const gsl_vector* gamma);

double gamma_gradient_interaction(double gamma_val,
                                  double gamma_sum,
                                  const gsl_vector* gamma,
                                  const gsl_vector* nu_interaction,
                                  unsigned int component);


/*
 * topic - the topic corresponding to the nu vector we're considering
 * gradient_component - the component we're taking the gradient of
 */
double nu_gradient_interaction(double nu_val,
                               double nu_sum,
                               const gsl_vector* nu,
                               gsl_matrix* nu_gamma_interaction,
                               unsigned int topic,
                               unsigned int gradient_component);

double nu_double_gradient_second_term(nu_params* p,
                                      const gsl_vector* nu,
                                      double nu_sum);

double nu_double_gradient_j(nu_params* p,
                            const gsl_vector* nu,
                            unsigned int j,
                            double nu_sum);

double nu_gradient_second_term(nu_params* p,
                               const gsl_vector* nu,
                               double nu_sum);

double nu_gradient_j(nu_params* p,
                     const gsl_vector* nu,
                     int j,
                     double second_term,
                     double nu_sum);

double nu_gradient_j(nu_params* p, const gsl_vector* nu, int j, double nu_sum);

void log_nu_gradient(const gsl_vector* x,
                     void* params,
                     gsl_vector* gradient);

void nu_gradient(const gsl_vector* nu,
                 void* params,
                 gsl_vector* gradient);

double beta_prior_lhood(gsl_vector* T, double alpha);

double beta_trans_lhood(gsl_vector* beta,
                        gsl_vector* doc_terms,
                        gsl_vector* trans_terms,
                        int num_docs,
                        bool ignore_trans,
                        double alpha_doc,
                        double alpha_trans);

/*
 * This is the part of lhood bound that comes from both beta and gamma
 * as well as nu
 */
double beta_trans_lhood(gsl_vector* beta,
                        gsl_matrix* gamma,
                        gsl_matrix* nu,
                        double alpha_doc,
                        double alpha_trans);

/*
 * This is the part of the lhood bound that just comes from beta and gamma
 */
double beta_doc_lhood(gsl_vector* beta,
                      gsl_matrix* gamma,
                      double alpha_doc);

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
                   gsl_vector* grad);

void logspace_beta(const gsl_vector* full_beta,
                   gsl_vector* logspace_beta);

/*
 * Creates T and the full beta from the truncated beta representation.
 * Because sometimes we only need T, if full_beta is NULL we don't try
 * to save the full_beta.
 */
void create_T_full_beta(const gsl_vector* trunc_beta_source,
                        const unsigned int num_topics,
                        gsl_vector* full_beta,
                        gsl_vector* T,
                        bool log_space = false);

double beta_lhood(const gsl_vector* beta, void* params);

void beta_lhood_gradient(const gsl_vector* beta,
                         void* params,
                         double* lhood,
                         gsl_vector* grad);

void beta_gradient_wrapper(const gsl_vector* beta,
                           void* params,
                           gsl_vector* grad);

double nu_token_lhood(const gsl_vector* nu, void* params);
double nu_lhood(const gsl_vector* nu, void* params);


double log_nu_lhood(const gsl_vector* x, void* params);

void gamma_gradient_base(const gsl_vector* gamma,
                         const gsl_vector* alpha,
                         const gsl_vector* phi_sum,
                         const gsl_vector* nu_interaction,
                         gsl_vector* gradient);

double gamma_lhood_base(const gsl_vector* gamma,
                        const gsl_vector* alpha,
                        const gsl_vector* phi_sum,
                        const gsl_vector* nu_interaction);

double gamma_lhood(const gsl_vector* gamma, void* params);

void gamma_lhood_gradient(const gsl_vector* gamma,
                          void* params,
                          double* lhood,
                          gsl_vector* grad);

void gamma_gradient(const gsl_vector* gamma,
                    void* params,
                    gsl_vector* gradient);

void log_nu_lhood_gradient(const gsl_vector* x,
                       void* params,
                       double* lhood,
                           gsl_vector* grad);

void nu_lhood_gradient(const gsl_vector* nu,
                       void* params,
                       double* lhood,
                       gsl_vector* grad);

void beta_const_terms(const gsl_matrix* gamma,
                      const gsl_matrix* nu,
                      const int num_topics,
                      gsl_vector* doc_terms,
                      gsl_vector* trans_terms);

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
                const int size);

double optimize_beta(gsl_vector* beta,
                     gsl_matrix* gamma,
                     gsl_matrix* nu,
                     bool ignore_docs,
                     bool ignore_trans,
                     double alpha_top,
                     double alpha_doc,
                     double alpha_trans,
                     string model_name);

/*
 * The second optional argument allows you to ignore all but one topic
 * for optimization
 */
double optimize_nu(VariationalParameters* vars,
                   const SyntopParameters* params,
                   int ignore_other_topics=-1);

// Returns true if the value improved
double optimize_gamma(gsl_vector* gamma,
                      gsl_vector* phi_total,
                      gsl_vector* nu_interaction,
                      string model_name,
                      gsl_vector* beta,
                      double alpha_doc,
                      int doc_num);

#endif
