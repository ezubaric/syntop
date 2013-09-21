#ifndef VARIATIONAL_PARAMETERS_INCLUDED
#define VARIATIONAL_PARAMETERS_INCLUDED

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <boost/shared_ptr.hpp>

#include "syntop_parameters.pb.h"
#include "vectorops.h"
#include "util.h"

using std::endl;
using boost::shared_ptr;

double matrix_accumulate(gsl_matrix* m, int ii, int jj, double val);
double vector_accumulate(gsl_vector* m, int ii, double val);

struct VariationalParameters {
  // Likelihood estimate (in a matrix, even though it's a scalar, in
  // order to parallelize better)
  shared_ptr<gsl_vector> lhood_;

  // Transition probabilities
  shared_ptr<gsl_matrix> nu_;
  shared_ptr<gsl_matrix> nu_exp_topic_counts_;
  shared_ptr<gsl_matrix> nu_exp_doc_counts_;

  // Per document topic distribution
  shared_ptr<gsl_matrix> gamma_;

  // Topic output parameters
  shared_ptr<gsl_matrix> tau_;
  shared_ptr<gsl_matrix> tau_est_top_;
  shared_ptr<gsl_vector> tau_est_bottom_;

  // Top-level stick-breaking weights
  shared_ptr<gsl_vector> beta_;
  shared_ptr<gsl_matrix> topic_counts_;
  int iteration;

  void InitializeParameters(const SyntopParameters& par) {
    //void InitializeParameters(SyntopParameters& par) {
    int number_documents = par.num_docs();
    int trunc = par.num_topics();

    // This is because we have seven components of the lhood that we sum
    // together; this is a structural number and independent of the
    // model
    lhood_.reset(gsl_vector_alloc(1), gsl_vector_free);
    gsl_vector_set_all(lhood_.get(), 0);

    // We have one extra position because of the need to have a special
    // start state
    nu_.reset(gsl_matrix_alloc(trunc + 1, trunc),
              gsl_matrix_free);

    topic_counts_.reset(gsl_matrix_alloc(trunc, 1), gsl_matrix_free);

    // These matrices store information needed for updating nu
    nu_exp_topic_counts_.reset(gsl_matrix_alloc(trunc + 1,
                                                trunc),
                               gsl_matrix_free);
    nu_exp_doc_counts_.reset(gsl_matrix_alloc(trunc + 1,
                                              trunc),
                             gsl_matrix_free);

    gsl_matrix_set_all(nu_exp_doc_counts_.get(), 0);
    gsl_matrix_set_all(nu_exp_topic_counts_.get(), 0);

    gsl_matrix_set_all(nu_.get(), par.alpha_trans());

    gamma_.reset(gsl_matrix_alloc(number_documents, trunc),
               gsl_matrix_free);

    if(par.ignore_docs()) {
      gsl_matrix_set_all(gamma_.get(), 1.0);
    } else {
      gsl_matrix_set_all(gamma_.get(), par.alpha_doc());
    }

    tau_.reset(gsl_matrix_alloc(trunc,
                                par.vocab_size()));
    tau_est_bottom_.reset(gsl_vector_alloc(trunc), gsl_vector_free);
    tau_est_top_.
      reset(gsl_matrix_alloc(trunc,
                             par.vocab_size()));

    gsl_matrix_set_all(tau_.get(), 0);
    gsl_vector_set_all(tau_est_bottom_.get(), 0);
    gsl_matrix_set_all(tau_est_top_.get(), 0);

    beta_.reset(gsl_vector_alloc(trunc), gsl_vector_free);
    gsl_vector_set_all(beta_.get(), 1.0);

    // display matrix to the console
    //    display_matrix(nu_.get(), "initial nu");
    //    display_matrix(gamma_.get(), "initial gamma");
    //    display_matrix(beta_.get(), "initial beta");
  }

  VariationalParameters(const SyntopParameters& par) {
    //  VariationalParameters(SyntopParameters& par) {
    InitializeParameters(par);
    
    if(!par.ignore_docs()) {
      // gamma
      mtx_fscanf((par.model_name() + ".gamma").c_str(), gamma_.get());
    }
    // display_matrix(gamma_.get(), "initialize gamma in reducer\n");
    cout << "successfully load gamma matrix..." << endl;

    if(!par.ignore_trans()) {
      // nu
      mtx_fscanf((par.model_name() + ".nu").c_str(), nu_.get());
      // display_matrix(nu_.get(), "initialize nu in reducer\n");
    }
    cout << "successfully load nu matrix..." << endl;

    // tau
    mtx_fscanf((par.model_name() + ".tau").c_str(), tau_.get());
    // display_matrix(tau_.get(), "initialize tau in reducer\n");
    cout << "successfully load tau matrix..." << endl;

    // beta
    if(!par.finite()) {
      vct_fscanf((par.model_name() + ".beta").c_str(), beta_.get());
      // display_vector(beta_.get(), "initialize beta in reducer\n");
    }
    cout << "successfully load beta matrix..." << endl;
  }
};

#endif  // VARIATIONAL_PARAMETERS_INCLUDED
