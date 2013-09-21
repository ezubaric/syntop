/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef DOCUMENTMAPPER_INCLUDED
#define DOCUMENTMAPPER_INCLUDED

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

#include "util.h"
#include "vectorops.h"
#include "variational_parameters.h"
#include "gradient.h"

#include "syntop_parameters.pb.h"

#define MAX_SENT 100
#define MAX_DOC_ITER 2
#define ALWAYS_PLOT_DOC 0

using lib_prob::safe_log;

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */
class DocumentMapper {
 public:
  DocumentMapper(VariationalParameters* vars,
                 const SyntopParameters* params);
  ~DocumentMapper();
  double EStep(double tolerance, Document* d, int doc);
  void Emit(Document* d, StringMap* count);

  // These are substeps:
  void UpdatePhi(Document* d);
  void UpdateOmega(Document* d);
  void UpdateBeta();
  void LazyUpdateNu(int topic);
  void UpdateNu();
  void UpdateRho();
  void UpdateTau();

  double MStep();
  void RunIterations(int nrounds, double tolerance);
  // double Likelihood();
  double Likelihood(Document* d, gsl_vector* lhood = NULL);
  double phi_lhood(Document* d, gsl_vector* lhood = NULL);
  double gamma_likelihood(Document* d, gsl_vector* lhood = NULL);
  void InitializeParameters();
  void ResetLhoodFile();
  double GlobalLikelihoodTerm(bool print_lhood = false);
  void WriteModel();
  void ReadModel(bool resume);

  double get_alpha_doc_();
  double get_alpha_trans_();
  double get_alpha_trans_self_();
  int get_total_words_();
  int get_number_topics_();

  void ClearProgressFiles();

  gsl_vector* get_beta_();

  gsl_matrix* get_trans_counts_();
  gsl_matrix* get_nu_interaction_();
  gsl_matrix* get_nu_();
  gsl_matrix* get_gamma_();

  string get_model_name_();

  double ReduceAndMaximize();

  static void InitializeBeta(gsl_vector* beta, double alpha_top);

  static double NuLhood(const gsl_matrix* nu,
                        const gsl_matrix* nu_exp_topic_counts,
                        const gsl_matrix* nu_exp_doc_counts,
                        const gsl_vector* beta,
                        bool finite,
                        double alpha_top,
                        double alpha_trans,
                        double alpha_self_weight);

  static double VocParUpdate(double numerator, double denominator,
                             double sigma);

  static void AccumulateMatrixFromFiles(gsl_matrix* m,
                                        const string base,
                                        const string extension,
                                        const int num_files);

  static double GlobalWeightTerm(const gsl_matrix* nu,
                                 const gsl_vector* beta,
                                 double alpha_trans,
                                 double alpha_top,
                                 double self_weight,
                                 bool finite);

  static double GlobalVocabularyTerm(const gsl_matrix* tau,
                                     const gsl_matrix* rho,
                                     const double vocab_sigma,
                                     const int number_terms,
                                     const int number_roles,
                                     const int number_topics,
                                     const bool ignore_tags);

  // Static functions for easier testing
  static double ComputeNewPhi(Document* d,
                              const int word,
                              const int topic,
                              const gsl_matrix* nu,
                              const gsl_matrix* tau,
                              const gsl_vector* gamma,
                              int number_terms,
                              int number_topics);

  static double ComputeNewPhiIntrinsic(const int topic,
                                       const gsl_vector* gamma,
                                       double gamma_normalizer);

  static double ComputeNewPhiRoot(const gsl_matrix* nu,
                                    const int topic,
                                  const int number_topics);

  static double ComputeNewPhiParent(const gsl_vector* parent_phi,
                                    const gsl_matrix* nu,
                                    const int topic,
                                    const int number_topics);

  static double ComputeNewPhiChildren(Document* d,
                                      LensNode* w,
                                      const gsl_matrix* nu,
                                      const gsl_vector* gamma,
                                      const double gamma_normalizer,
                                      const int topic,
                                      const int number_topics);

  static double ComputeNewPhiVocab(const gsl_matrix* tau,
                                   const gsl_matrix* rho,
                                   const Document* d,
                                   const LensNode* w,
                                   const int word_index,
                                   const int topic,
                                   const int number_terms);

  static double LikelihoodTopic(Document* d,
                                const int num_topics);

  static double LikelihoodVocab(Document* d,
                                const gsl_matrix* tau,
                                const int number_topics,
                                const int num_words);

  static double Likelihood(Document* d,
                           const gsl_vector* gamma,
                           const double alpha_doc,
                           const gsl_vector* beta,
                           const gsl_matrix* nu,
                           const int num_topics,
                           const int num_terms,
                           const gsl_matrix* tau,
                           gsl_vector* lhood_total);

  static double phi_lhood(Document* d,
                           const gsl_vector* gamma,
                           const double alpha_doc,
                           const gsl_vector* beta,
                           const gsl_matrix* nu,
                           const int num_topics,
                           const int num_terms,
                           const gsl_matrix* tau,
                           gsl_vector* lhood_total);

  static double gamma_likelihood(Document* d,
                           const gsl_vector* gamma,
                           const double alpha_doc,
                           const gsl_vector* beta,
                           const gsl_matrix* nu,
                           const int num_topics,
                           const int num_terms,
                           const gsl_matrix* tau,
                           gsl_vector* lhood_total);

  static void CumulateGammaStatistics(Document* d,
                                      const gsl_matrix* nu,
                                      gsl_vector* phi_total_vector,
                                      gsl_vector* nu_interaction_vector);

  static void CumulateNuStatistics(Document* d,
                                   const int doc,
                                   const gsl_vector* gamma,
                                   gsl_matrix* nu_exp_topic_counts,
                                   gsl_matrix* nu_exp_doc_counts);

  static void CumulateTauUpdate(Document* d,
                                const int number_topics,
                                const int number_terms,
                                gsl_matrix* tau_est_top,
                                gsl_vector* tau_est_bottom);

 private:
  VariationalParameters* vars_;
  const SyntopParameters* params_;
  DocumentMapper();  // No evil constructors
};

#endif
