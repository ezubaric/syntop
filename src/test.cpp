/*
 * Copyright 2007 Jordan Boyd-Graber
 *
 * Tests correctness of the lensing code.  Probably should be more
 * atomic, but that would have also required more stubbing (which also
 * would have been good).
 *
 * TODO: put setting of variables into a seperate function
 */

#include <boost/test/included/unit_test_framework.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "gradient.h"
#include "syntop_mapper.h"
#include "syntop_reducer.h"
#include "document_mapper.h"
#include <gsl/gsl_vector.h>

#define CLOSE_TOL 0.001
#define LOOSE_TOL 1

using boost::unit_test_framework::test_suite;

Document* load_test_doc(string line, int num_topics, int id) {
  Document* d = new Document(id);
  d->SetNumberSentences(1);
  d->AddSentence(0, line);
  d->FillWordArray();
  d->ResetPhi(num_topics);
  // d->InitializeVariationalParameters(num_topics);

  return d;
}

void test_log_functions() {
  double neg_infinity = -std::numeric_limits<double>::infinity();
  double zero = 0;
  BOOST_CHECK_CLOSE(log_sum(neg_infinity, zero),
                    zero,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(log_sum(zero, neg_infinity),
                    zero,
                    CLOSE_TOL);

  double double_nan = log_sum(neg_infinity, neg_infinity);
  cout << "Two sums are: " << double_nan << endl;
  BOOST_CHECK(!isnan(double_nan));
  BOOST_CHECK(double_nan < -1000);

  shared_ptr<gsl_vector> v(gsl_vector_alloc(6), gsl_vector_free);
  gsl_vector_set(v.get(), 0, neg_infinity);
  gsl_vector_set(v.get(), 1, neg_infinity);
  gsl_vector_set(v.get(), 2, -5.53);
  gsl_vector_set(v.get(), 3, neg_infinity);
  gsl_vector_set(v.get(), 4, -9.08);
  gsl_vector_set(v.get(), 5, neg_infinity);


  display_vector(v.get(), "Unnormalized ");
  double sum = log_normalize(v.get());
  cout << "Sum = " << sum << endl;
  BOOST_CHECK(!isnan(sum));
  display_vector(v.get(), "Normalized ");

  BOOST_CHECK(gsl_vector_get(v.get(), 1) <= safe_log(0.0));
  BOOST_CHECK_CLOSE(gsl_vector_get(v.get(), 2), -0.02832, CLOSE_TOL);
  BOOST_CHECK(gsl_vector_get(v.get(), 3) <= safe_log(0.0));

  vexp(v.get(), v.get());
  display_vector(v.get(), "Exponentiated ");

  BOOST_CHECK_CLOSE(gsl_vector_get(v.get(), 1), 0.0, CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_vector_get(v.get(), 2), 0.972077, CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_vector_get(v.get(), 3), 0.0, CLOSE_TOL);

  gsl_sf_lngamma(32e6);
}

void test_tau_rho_setting() {
  int num_topics = 2;
  int num_roles = 2;
  int num_terms = 2;

  scoped_ptr<Document> d(load_test_doc("1 0 1 0 1 0 -1 -1", num_topics, 0));

  //                            (0.5, 0.5)
  LensNode* n1 = d->get_head_(0);
  gsl_vector_set(n1->phi_.get(), 0, 0.5);
  gsl_vector_set(n1->phi_.get(), 1, 0.5);
  //       (1.0, 0.0)
  LensNode* n2 = n1->GetChild(0);
  gsl_vector_set(n2->phi_.get(), 0, 1.0);
  gsl_vector_set(n2->phi_.get(), 1, 0.0);

  shared_ptr<gsl_matrix> rho(gsl_matrix_alloc(num_roles, num_terms));
  gsl_matrix_set(rho.get(), 0, 0, 1.0/3.0);
  gsl_matrix_set(rho.get(), 0, 1, 2.0/3.0);
  gsl_matrix_set(rho.get(), 1, 0, 0.75);
  gsl_matrix_set(rho.get(), 1, 1, 0.25);
  shared_ptr<gsl_matrix> tau(gsl_matrix_alloc(num_topics, num_terms));
  gsl_matrix_set(tau.get(), 0, 0, 1.0);
  gsl_matrix_set(tau.get(), 0, 1, 0);
  gsl_matrix_set(tau.get(), 1, 0, 0.5);
  gsl_matrix_set(tau.get(), 1, 1, 0.5);

  shared_ptr<gsl_matrix> tau_est_top(gsl_matrix_alloc(num_topics, num_terms));
  shared_ptr<gsl_vector> tau_est_bottom(gsl_vector_alloc(num_topics));

  gsl_matrix_set_all(tau_est_top.get(), 0.0);
  gsl_vector_set_all(tau_est_bottom.get(), 0.0);

  DocumentMapper::CumulateTauUpdate(d.get(),
                               num_topics,
                               num_terms,
                               tau_est_top.get(),
                               tau_est_bottom.get());

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 0, 0),
                    1.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 0, 1),
                    0.5,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 1, 0),
                    0.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 1, 1),
                     0.5,
                     CLOSE_TOL);

  // 1/2*1/2*exp(1/3)+1/4*1*exp(3/4)
  BOOST_CHECK_CLOSE(gsl_vector_get(tau_est_bottom.get(), 0),
                    0.8781531,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_vector_get(tau_est_bottom.get(), 1),
                    0.8079399,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_vector_get(tau_est_bottom.get(), 2),
                    0.3489031,
                    CLOSE_TOL);

  // Now we'll see if the likelihood calculation makes sense
  // 1/2*(0)+1/2*(1/2)+1*(1)+0*(1/2) -
  //   1/2*(1/2*exp(1+1/3)+1/2*exp(1/2+1/3)+1/2*exp(0+2/3)+1/2*exp(1/2+2/3)) -
  //   1/4*(1*exp(1+3/4)+0+1*exp(0+1/4)+0) = -3.169425

  // These methods seem to have disappeared

  /*

  BOOST_CHECK_CLOSE(DocumentMapper::DocumentTauLikelihood(tau.get(),
                                                          tau_est_top.get(),
                                                          tau_est_bottom.get(),
                                                          false),
                    -3.323069,
                    CLOSE_TOL);

  double vocab_sigma = 2.0;


  // Test to make sure that LambertW is working the way we think it is
  BOOST_CHECK_CLOSE(gsl_sf_lambert_W0(1.0),
                    .5671432904,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_sf_lambert_W0(1.6151101794618921),
                    .7573426554,
                    CLOSE_TOL);

  double a = 2.0;
  double b = 4.0;
  double c = 2.0;

  // solve(2.0-4.0*exp(x)-1/4*x=0, x); [Maple command]
  BOOST_CHECK_CLOSE(DocumentMapper::VocParUpdate(a, b, c),
                    -0.6186591954,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(DocumentMapper::VocParUpdate(gsl_matrix_get(tau_est_top.get(), 0, 0),
                                            gsl_matrix_get(tau_est_bottom.get(), 0, 0),
                                            vocab_sigma),
                    0.1036740189,
                    CLOSE_TOL);
  */
}

void test_lhood() {
  int num_topics = 2;
  int num_roles = 2;
  int num_terms = 2;

  scoped_ptr<Document> d(load_test_doc("1 0 1 0 1 0 -1 -1", num_topics, 0));

  // Should we check the initialization somehow?
  BOOST_CHECK_EQUAL(d->get_number_words_(), 2);
  gsl_vector_set(d->get_omega_(), 0, 1.0/1.0);
  gsl_vector_set(d->get_omega_(), 1, 1.0/2.0);

  //                            (0.5, 0.5)
  LensNode* n1 = d->get_head_(0);
  gsl_vector_set(n1->phi_.get(), 0, 0.5);
  gsl_vector_set(n1->phi_.get(), 1, 0.5);
  //       (1.0, 0.0)
  LensNode* n2 = n1->GetChild(0);
  gsl_vector_set(n2->phi_.get(), 0, 1.0);
  gsl_vector_set(n2->phi_.get(), 1, 0.0);

  shared_ptr<gsl_matrix> rho(gsl_matrix_alloc(num_roles, num_terms));
  gsl_matrix_set(rho.get(), 0, 0, 1.0/3.0);
  gsl_matrix_set(rho.get(), 0, 1, 2.0/3.0);
  gsl_matrix_set(rho.get(), 1, 0, 0.75);
  gsl_matrix_set(rho.get(), 1, 1, 0.25);
  shared_ptr<gsl_matrix> tau(gsl_matrix_alloc(num_topics, num_terms));
  gsl_matrix_set(tau.get(), 0, 0, 1.0);
  gsl_matrix_set(tau.get(), 0, 1, 0);
  gsl_matrix_set(tau.get(), 1, 0, 0.5);
  gsl_matrix_set(tau.get(), 1, 1, 0.5);
  shared_ptr<gsl_matrix> nu(gsl_matrix_alloc(3, 2), gsl_matrix_free);
  gsl_matrix_set(nu.get(), 0, 0, 0.5);
  gsl_matrix_set(nu.get(), 0, 1, 0.5);
  gsl_matrix_set(nu.get(), 1, 0, 0.9);
  gsl_matrix_set(nu.get(), 1, 1, 0.1);
  gsl_matrix_set(nu.get(), 2, 0, 0.25);
  gsl_matrix_set(nu.get(), 2, 1, 0.75);
  shared_ptr<gsl_vector> gamma(gsl_vector_alloc(2), gsl_vector_free);
  gsl_vector_set(gamma.get(), 0, 0.2);
  gsl_vector_set(gamma.get(), 1, 0.8);
  shared_ptr<gsl_vector> beta(gsl_vector_alloc(2), gsl_vector_free);
  gsl_vector_set(beta.get(), 0, 2.0);
  gsl_vector_set(beta.get(), 1, 1.0);

  double alpha_doc = 0.5;

  double dirichlet_lhood = var_dir_lhood(alpha_doc,
                                         beta.get(),
                                         gamma.get());

  // lgamma(0.5*2 + 0.5*1)
  //   - lgamma(0.5*2) - lgamma(0.5*1)
  //   + (0.5*2-1)*(digamma(2/10.) - digamma(1.))
  //   + (0.5*1-1)*(digamma(8/10.) - digamma(1.))
  BOOST_CHECK_CLOSE(dirichlet_lhood,
                    -0.4992507,
                    CLOSE_TOL);

  // gamma_lhood_exp_counts
  // .5*(digamma(.2)-digamma(1.0)) +
  //    .5*(digamma(.8)-digamma(1.0)) +
  //    1.0*(digamma(.2)-digamma(1.0))
  double gamma_lhood = gamma_lhood_exp_counts(d.get(),
                                              gamma.get());

  BOOST_CHECK_CLOSE(gamma_lhood,
                    -7.261633,
                    CLOSE_TOL);

  // nu_doc_transition_likelihood

  //  1*0.5*(digamma(0.5)-digamma(1.0))
  double transition_1 = nu_doc_transition_likelihood(d.get(),
                                                     &gsl_matrix_row
                                                     (nu.get(), 0).vector,
                                                     0);

  BOOST_CHECK_CLOSE(transition_1,
                    -0.6931472,
                    CLOSE_TOL);

  // 1*0.5*(digamma(.9)-digamma(1.0))
  double transition_2 = nu_doc_transition_likelihood(d.get(),
                                                     &gsl_matrix_row
                                                     (nu.get(), 1).vector,
                                                     1);

  BOOST_CHECK_CLOSE(transition_2,
                    -0.08885564,
                    CLOSE_TOL);

  // 0.5*(digamma(.25) - digamma(.75)) + 0.5*(digamma(.75)-digamma(.75))
  double transition_3 = nu_doc_transition_likelihood(d.get(),
                                                     &gsl_matrix_row
                                                     (nu.get(), 2).vector,
                                                     2);


  BOOST_CHECK_CLOSE(transition_3,
                    -2.079442,
                    CLOSE_TOL);

  double gamma_nu_interaction_lhood = gamma_nu_interaction(d.get(),
                                                           gamma.get(),
                                                           nu.get());

  // 1/1*(2/10*1/4 + 8/10*3/4) + safe_log(1) - 1 +
  // 1/2*(1/2*2/10*1/2+1/2*2/10*9/10+1/2*8/10*1/2+1/2*8/10*1/10) + safe_log(2.)
  // - 1
  BOOST_CHECK_CLOSE(gamma_nu_interaction_lhood,
                    -0.4668528,
                    CLOSE_TOL);

  double vocab_lhood = DocumentMapper::LikelihoodVocab(d.get(),
                                                  tau.get(),
                                                  num_topics,
                                                  num_terms);

  // 1/2.*(0+2/3.)+1/2.*(1/2.+2/3.)+1*(1+3/4.)
  // - (1/2.*(1/2.*(exp(1/3.+1)+exp(2/3.+0))+1/2.*(exp(1/3.+1/2.) +
  //                exp(2/3+1/2)))+safe_log(2.0)-1)
  // - (1/4.*(1*(exp(3/4.+1)+exp(1/4.+0))) + safe_log(4.0) - 1)
  BOOST_CHECK_CLOSE(vocab_lhood,
                    -1.985844,
                    CLOSE_TOL);

  double variational_dirichlet = var_dir_lhood(gamma.get(),
                                               gamma.get());

  // lgamma(2/10. + 8/10.)
  //    - lgamma(.2) - lgamma(.8)
  //    + (.2-1)*(digamma(2/10.) - digamma(1.))
  //    + (.8-1)*(digamma(8/10.) - digamma(1.))
  BOOST_CHECK_CLOSE(variational_dirichlet,
                    2.170894,
                    CLOSE_TOL);

  // 0.5 * safe_log(0.5) + 0.5*safe_log(0.5) + 1 * safe_log(1.0)
  double vocab_topic = DocumentMapper::LikelihoodTopic(d.get(),
                                                  num_topics);
  BOOST_CHECK_CLOSE(vocab_topic,
                    -0.6931472,
                    CLOSE_TOL);

  double final = dirichlet_lhood
    + gamma_lhood + transition_1 + transition_2 +transition_3
    - gamma_nu_interaction_lhood
    + vocab_lhood
    - variational_dirichlet
    - vocab_topic;
  BOOST_CHECK_CLOSE(final,
                    DocumentMapper::Likelihood(d.get(),
                                          gamma.get(),
                                          alpha_doc,
                                          beta.get(),
                                          nu.get(),
                                          num_topics,
                                          num_terms,
                                          tau.get(),
                                          NULL),
                    CLOSE_TOL);


  double vocab_sigma = 2.0;

  // This is now in a different class
  /*

  double vocab_param_lhood = DocumentMapper::GlobalVocabularyTerm(tau.get(),
                                                             rho.get(),
                                                             vocab_sigma,
                                                             num_terms,
                                                             num_roles,
                                                             num_topics,
                                                             false);

  // -((1/3)^2+(2/3)^2+(3/4)^2+(1/4)^2+(1)^2+(0)^2+(1/2)^2+(1/2)^2) /
  //             (2 * 2 * 2)
  BOOST_CHECK_CLOSE(vocab_param_lhood,
                    -0.3350695,
                    CLOSE_TOL);
  */


  // TODO(zhaike): Add tests for when we're ignoring the tags
}

void test_document_creation() {
  scoped_ptr<Document> d(load_test_doc("2 0 2 1 1 0 -1 0 1 0 -1 -1", 1, 0));

  d->FillWordArray();
  BOOST_CHECK_EQUAL(d->get_number_words_(), 3);

  LensNode* n1 = d->get_head_(0);
  LensNode* n2 = n1->GetChild(0);
  LensNode* n3 = n1->GetChild(1);
  BOOST_CHECK_EQUAL(n1->number_children_, 2);
  BOOST_CHECK_EQUAL(n2->number_children_, 0);
  BOOST_CHECK_EQUAL(n3->number_children_, 0);
  BOOST_CHECK_EQUAL(n3->term_, 0);
  BOOST_CHECK_EQUAL(n2->term_, 1);
  BOOST_CHECK_EQUAL(n1->term_, 2);
  BOOST_CHECK_EQUAL(n1->relation_, 0);
  BOOST_CHECK_EQUAL(n2->relation_, 1);
  BOOST_CHECK_EQUAL(n3->relation_, 1);
}

/*
void test_beta_initialization() {
  int num_topics = 3;

  shared_ptr<gsl_vector> beta_full(gsl_vector_alloc(num_topics),
                                   gsl_vector_free);
  double alpha_top = 3.0;

  DocumentMapper::InitializeBeta(beta_full.get(), alpha_top);

  BOOST_CHECK_CLOSE(gsl_vector_get(beta_full.get(), 0),
                    16.0/21.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_vector_get(beta_full.get(), 1),
                    4.0/21.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_vector_get(beta_full.get(), 2),
                    1.0/21.0,
                    CLOSE_TOL);
}
*/

void test_stick_breaking() {
  int num_topics = 3;
  int num_docs = 1;

  double alpha_doc = 2.0;
  double alpha_trans = 1.0;

  // Creating T and the full beta
  shared_ptr<gsl_vector> beta_full(gsl_vector_alloc(num_topics),
                                   gsl_vector_free);
  shared_ptr<gsl_vector> beta_full_doc(gsl_vector_alloc(num_topics),
                                       gsl_vector_free);
  shared_ptr<gsl_vector> beta_full_trans(gsl_vector_alloc(num_topics),
                                         gsl_vector_free);

  shared_ptr<gsl_vector> beta_trunc(gsl_vector_alloc(num_topics - 1),
                                    gsl_vector_free);
  shared_ptr<gsl_vector> T(gsl_vector_alloc(num_topics),
                           gsl_vector_free);
  shared_ptr<gsl_vector> log_beta(gsl_vector_alloc(num_topics - 1),
                                  gsl_vector_free);

  gsl_vector_set(beta_trunc.get(), 0, 0.5);
  gsl_vector_set(beta_trunc.get(), 1, 1/3.0);

  create_T_full_beta(beta_trunc.get(),
                     num_topics,
                     beta_full.get(),
                     T.get());

  logspace_beta(beta_full.get(), log_beta.get());

  display_vector(beta_trunc.get(), "beta trunc");
  BOOST_CHECK_CLOSE(gsl_vector_get(beta_full.get(), 0), 1/2.0, CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_vector_get(beta_full.get(), 1), 1/3.0, CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_vector_get(beta_full.get(), 2), 1/6.0, CLOSE_TOL);

  gsl_vector_set(beta_full_doc.get(), 0, alpha_doc*1/2.0);
  gsl_vector_set(beta_full_doc.get(), 1, alpha_doc*1/3.0);
  gsl_vector_set(beta_full_doc.get(), 2, alpha_doc*1/6.0);
  gsl_vector_set(beta_full_trans.get(), 0, alpha_trans*1/2.0);
  gsl_vector_set(beta_full_trans.get(), 1, alpha_trans*1/3.0);
  gsl_vector_set(beta_full_trans.get(), 2, alpha_trans*1/6.0);

  BOOST_CHECK_CLOSE(gsl_vector_get(T.get(), 0), 1.0,     CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_vector_get(T.get(), 1), 1/2.0,   CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_vector_get(T.get(), 2), 1/6.0,   CLOSE_TOL);

  shared_ptr<gsl_matrix> gamma(gsl_matrix_alloc(num_docs, num_topics),
                               gsl_matrix_free);
  gsl_matrix_set(gamma.get(), 0, 0, 1/6.0);
  gsl_matrix_set(gamma.get(), 0, 1, 1/4.0);
  gsl_matrix_set(gamma.get(), 0, 2, 7/12.0);
  shared_ptr<gsl_matrix> nu(gsl_matrix_alloc(num_topics + 1, num_topics),
                            gsl_matrix_free);
  gsl_matrix_set_all(nu.get(), 1/3.0);
  gsl_matrix_set(nu.get(), 1, 0, 1/2.0);
  gsl_matrix_set(nu.get(), 1, 2, 1/6.0);
  gsl_matrix_set(nu.get(), 2, 0, 1/3.0);
  gsl_matrix_set(nu.get(), 2, 1, 1/4.0);
  gsl_matrix_set(nu.get(), 2, 2, 5/12.0);
  gsl_matrix_set(nu.get(), 3, 0, 1/9.0);
  gsl_matrix_set(nu.get(), 3, 1, 1/9.0);
  gsl_matrix_set(nu.get(), 3, 2, 7/9.0);

  // loading constant terms

  shared_ptr<gsl_vector> doc_terms(gsl_vector_alloc(num_topics),
                                   gsl_vector_free);
  shared_ptr<gsl_vector> trans_terms(gsl_vector_alloc(num_topics),
                                     gsl_vector_free);
  beta_const_terms(gamma.get(),
                   nu.get(),
                   num_topics,
                   doc_terms.get(),
                   trans_terms.get());

  // digamma(1/6)-digamma(1)
  BOOST_CHECK_CLOSE(gsl_vector_get(doc_terms.get(), 0),
                    -5.754912,
                    CLOSE_TOL);

  // digamma(1/4)-digamma(1)
  BOOST_CHECK_CLOSE(gsl_vector_get(doc_terms.get(), 1),
                    -3.650238,
                    CLOSE_TOL);

  // digamma(7/12)-digamma(1)
  BOOST_CHECK_CLOSE(gsl_vector_get(doc_terms.get(), 2),
                    -1.025428,
                    CLOSE_TOL);

  // digamma(1/3)+digamma(1/2)+digamma(1/3)+digamma(1/9)-4*digamma(1.0)
  BOOST_CHECK_CLOSE(gsl_vector_get(trans_terms.get(), 0),
                    -15.32666,
                    CLOSE_TOL);

  // digamma(1/3)+digamma(1/3)+digamma(1/4)+digamma(1/9)-4*digamma(1.0)
  BOOST_CHECK_CLOSE(gsl_vector_get(trans_terms.get(), 1),
                    -17.5906,
                    CLOSE_TOL);

  // digamma(1/3)+digamma(1/6)+digamma(5/12)+digamma(7/9)-4*digamma(1.0)
  BOOST_CHECK_CLOSE(gsl_vector_get(trans_terms.get(), 2),
                    -10.61696,
                    CLOSE_TOL);

  // beta_lhood

  // This is where we test var_dir_lhood with scaling
  double lhood = 0.0;
  for (int i = 0; i <= num_topics; i++) {
    gsl_vector_view nu_i = gsl_matrix_row(nu.get(), i);
    lhood += var_dir_lhood(beta_full_trans.get(), &nu_i.vector);
  }
  gsl_vector_view gamma_d = gsl_matrix_row(gamma.get(), 0);
  lhood += var_dir_lhood(beta_full_doc.get(), &gamma_d.vector);

  // Testing scaling
  BOOST_CHECK_CLOSE(var_dir_lhood(alpha_doc,
                                  beta_full.get(),
                                  &gamma_d.vector),
                    var_dir_lhood(beta_full_doc.get(),
                                  &gamma_d.vector),
                    CLOSE_TOL);


  BOOST_CHECK_CLOSE(var_dir_lhood(beta_full_doc.get(), &gamma_d.vector),
                    beta_doc_lhood(beta_full.get(),
                                   gamma.get(),
                                   alpha_doc),
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(beta_trans_lhood(beta_full.get(),
                                     gamma.get(),
                                     nu.get(),
                                     alpha_doc,
                                     alpha_trans),
                    lhood,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(beta_trans_lhood(beta_full.get(),
                                     doc_terms.get(),
                                     trans_terms.get(),
                                     num_docs,
                                     false,
                                     alpha_doc,
                                     alpha_trans),
                    lhood,
                    CLOSE_TOL);

  double alpha_root = 3.0;

  double lhood_prior = beta_prior_lhood(T.get(),
                                        alpha_root);

  beta_params p;
  p.doc_terms = doc_terms.get();
  p.trans_terms = trans_terms.get();
  p.num_docs = num_docs;
  p.ignore_trans = false;
  p.alpha_top = alpha_root;
  p.alpha_doc = alpha_doc;
  p.alpha_trans = alpha_trans;
  p.full_beta = beta_full.get();
  p.T = T.get();

  // vsafe_log(beta_trunc.get(), beta_trunc.get());

  // We negate because gsl does minimization
  BOOST_CHECK_CLOSE(-(lhood_prior + lhood),
                    beta_lhood(log_beta.get(),
                               &p),
                    CLOSE_TOL);


  // (3-1)*safe_log(1/6)-safe_log(1.)-safe_log(1/2.)
  BOOST_CHECK_CLOSE(lhood_prior,
                    -2.890372,
                    CLOSE_TOL);

  // beta_lhood_gradient

  shared_ptr<gsl_vector> gradient(gsl_vector_alloc(num_topics-1),
                                  gsl_vector_free);
  beta_gradient_wrapper(log_beta.get(),
                        &p,
                        gradient.get());

  display_vector(log_beta.get(), "log beta");
  // A <- 2-(3-1)/(1/6)-4*(1.0*digamma(1.0*1/2)-1*digamma(1*1/6)) -
  // 1*(2*digamma(2*1/2)-2*digamma(2*1/6)) +
  // 2*-5.755+-2*(-1.0254)+1*-15.32666-1*(-10.61696)

  // Implicit term:
  // B <- exp(1.09861)
  // C <- exp(1.09861) + exp(0.693147) + exp(0.0)
  // A * (C * B - B * B) / (C * C)
  // = -11.68825
  BOOST_CHECK_CLOSE(-(-11.68825),
                    gsl_vector_get(gradient.get(), 0),
                    CLOSE_TOL);

  // D <- ??? (TODO: write derivation for dL / dBeta_2)
  // E <- exp(0.693147)
  // D * (C * E - E * E) / (C * C)
  BOOST_CHECK_CLOSE(-(-9.033611),
                    gsl_vector_get(gradient.get(), 1),
                    CLOSE_TOL);
}


void test_phi_setting() {
  int num_topics = 2;
  scoped_ptr<Document> d(load_test_doc("2 0 2 1 1 0 -1 0 1 0 -1 -1",
                                       num_topics, 0));

  // Should we check the initialization somehow?
  BOOST_CHECK_EQUAL(d->get_number_words_(), 3);
  gsl_vector_set(d->get_omega_(), 0, 1.0/1.0);
  gsl_vector_set(d->get_omega_(), 1, 1.0/2.0);
  gsl_vector_set(d->get_omega_(), 2, 1.0/3.0);

  //                            (0.25, 0.75)
  LensNode* n1 = d->get_head_(0);
  gsl_vector_set(n1->phi_.get(), 0, 0.25);
  gsl_vector_set(n1->phi_.get(), 1, 0.75);
  //       (1.0, 0.0)
  LensNode* n2 = n1->GetChild(0);
  gsl_vector_set(n2->phi_.get(), 0, 1.0);
  gsl_vector_set(n2->phi_.get(), 1, 0.0);
  //                                                 (0.5, 0.5)

  shared_ptr<gsl_vector> gamma(gsl_vector_alloc(2), gsl_vector_free);
  gsl_vector_set(gamma.get(), 0, 0.2);
  gsl_vector_set(gamma.get(), 1, 0.8);
  shared_ptr<gsl_matrix> nu(gsl_matrix_alloc(3, 2), gsl_matrix_free);
  gsl_matrix_set(nu.get(), 0, 0, 0.5);
  gsl_matrix_set(nu.get(), 0, 1, 0.5);
  gsl_matrix_set(nu.get(), 1, 0, 0.9);
  gsl_matrix_set(nu.get(), 1, 1, 0.1);
  gsl_matrix_set(nu.get(), 2, 0, 0.25);
  gsl_matrix_set(nu.get(), 2, 1, 0.75);
  shared_ptr<gsl_matrix> rho(gsl_matrix_alloc(2, 3));
  gsl_matrix_set(rho.get(), 0, 0, 0.1);
  gsl_matrix_set(rho.get(), 0, 1, 0.1);
  gsl_matrix_set(rho.get(), 0, 2, 0.8);
  gsl_matrix_set(rho.get(), 1, 0, 0.4);
  gsl_matrix_set(rho.get(), 1, 1, 0.4);
  gsl_matrix_set(rho.get(), 1, 2, 0.2);
  shared_ptr<gsl_matrix> tau(gsl_matrix_alloc(2, 3));
  gsl_matrix_set(tau.get(), 0, 0, 2.0/3.0);
  gsl_matrix_set(tau.get(), 0, 1, 1.0/6.0);
  gsl_matrix_set(tau.get(), 0, 2, 1.0/6.0);
  gsl_matrix_set(tau.get(), 1, 0, 0.2);
  gsl_matrix_set(tau.get(), 1, 1, 0.4);
  gsl_matrix_set(tau.get(), 1, 2, 0.4);

  LensNode* w = d->get_word_(0);
  double intrinsic = DocumentMapper::
    ComputeNewPhiIntrinsic(0, gamma.get(), gsl_blas_dsum(gamma.get()));

  // digamma(2/10.) - digamma(1.0)
  BOOST_CHECK_CLOSE(intrinsic, -4.711824, CLOSE_TOL);

  double children = DocumentMapper::
    ComputeNewPhiChildren(d.get(), w, nu.get(), gamma.get(),
                          gsl_blas_dsum(gamma.get()), 0, 2);

  // -1/2.0 * (2/20.+8/20) - 1/3.0 * (8/20+2/20.) +
  // 2.0 * (digamma(1/2.) - digamma(1.0))
  BOOST_CHECK_CLOSE(children, -3.189255, CLOSE_TOL);

  double vocab = gsl_matrix_get(tau.get(), 0, 0);

  // 1/6. + 8/10. - 1/2.*(exp(1/10.+2/3.) + exp(1/10.+1/6.) + exp(8/10.+1/6.))
  BOOST_CHECK_CLOSE(vocab, -2.0770088, CLOSE_TOL);

  // digamma(.25) - digamma(1.0)
  double root = DocumentMapper::
    ComputeNewPhiRoot(nu.get(), 0, 2);

  BOOST_CHECK_CLOSE(root, -3.650238, CLOSE_TOL);

  double final = DocumentMapper::
    ComputeNewPhi(d.get(), 0, 0, nu.get(), tau.get(), gamma.get(), 3, 2);

  // -4.711824 + -3.189255 + -2.077008 + -3.650238
  BOOST_CHECK_CLOSE(final, -13.62833, CLOSE_TOL);

  double parent = DocumentMapper::ComputeNewPhiParent(n1->phi_.get(),
                                                 nu.get(),
                                                 1,
                                                 2);

  // 0.25 * (digamma(.5) - digamma(1.0)) + 0.75 * (digamma(0.1) - digamma(1.0))
  BOOST_CHECK_CLOSE(parent, -7.731478, CLOSE_TOL);

  // Now we'll test to see if values get cummulated correctly
  shared_ptr<gsl_matrix> tau_est_top(gsl_matrix_alloc(2, 3),
                                     gsl_matrix_free);
  shared_ptr<gsl_vector> tau_est_bottom(gsl_vector_alloc(2),
                                        gsl_vector_free);
  gsl_matrix_set_all(tau_est_top.get(), 0.0);
  gsl_vector_set_all(tau_est_bottom.get(), 0.0);

  DocumentMapper::CumulateTauUpdate(d.get(),
                               2,
                               3,
                               tau_est_top.get(),
                               tau_est_bottom.get());


  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 0, 0),
                    0.5,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 1, 0),
                    0.5,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 0, 1),
                    1.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 1, 1),
                    0.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 0, 2),
                    0.25,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(tau_est_top.get(), 1, 2),
                    0.75,
                    CLOSE_TOL);
}

void test_gamma_lhood_gradient() {
  scoped_ptr<Document> d(load_test_doc("2 0 2 1 1 0 -1 0 1 0 -1 -1", 2, 0));
  //                            (0.25, 0.75)
  LensNode* n1 = d->get_head_(0);
  gsl_vector_set(n1->phi_.get(), 0, 0.25);
  gsl_vector_set(n1->phi_.get(), 1, 0.75);
  //       (1.0, 0.0)
  LensNode* n2 = n1->GetChild(0);
  gsl_vector_set(n2->phi_.get(), 0, 1.0);
  gsl_vector_set(n2->phi_.get(), 1, 0.0);
  //                                                 (0.5, 0.5)
  LensNode* n3 = n1->GetChild(1);
  // ResetPhi sets n3 to be uniform by default, so this lets us test
  // that assumption
  BOOST_CHECK_EQUAL(gsl_vector_get(n3->phi_.get(), 0), 0.5);
  BOOST_CHECK_EQUAL(gsl_vector_get(n3->phi_.get(), 1), 0.5);

  // var_dir_lhood is tested in nu's section, so we're not doing that
  // here

  // gamma nu interaction
  shared_ptr<gsl_vector> gamma(gsl_vector_alloc(2), gsl_vector_free);
  gsl_vector_set(gamma.get(), 0, 0.2);
  gsl_vector_set(gamma.get(), 1, 0.8);
  shared_ptr<gsl_matrix> nu(gsl_matrix_alloc(3, 2), gsl_matrix_free);
  gsl_matrix_set(nu.get(), 0, 0, 0.5);
  gsl_matrix_set(nu.get(), 0, 1, 0.5);
  gsl_matrix_set(nu.get(), 1, 0, 0.9);
  gsl_matrix_set(nu.get(), 1, 1, 0.1);
  gsl_matrix_set(nu.get(), 2, 0, 0.25);
  gsl_matrix_set(nu.get(), 2, 1, 0.75);
  gsl_vector_view nu1_view = gsl_matrix_row(nu.get(), 0);

  // Is this no longer needed?
  // d->InitializeVariationalParameters(2);

  gsl_vector_set(d->get_omega_(), 0, 1.0/1.0);
  gsl_vector_set(d->get_omega_(), 1, 1.0/2.0);
  gsl_vector_set(d->get_omega_(), 2, 1.0/3.0);

  /* Should be
   * 1*(2/10*1/4 + 8/10*3/4) + (1/3 + 1/2)*(1/4 * 1/2 * 2/10 +
   *                                1/4 * 4/5 * 1/2  +
   *                                3/4*2/10*9/10    +
   *                                3/4*4/5*1/10) +
   * safe_log(1) + safe_log(2) + safe_log(3) - 1.0 - 1.0 - 1.0
   */
  BOOST_CHECK_CLOSE(gamma_nu_interaction(d.get(), gamma.get(), nu.get()),
                    0.9166667 + safe_log(6.0) - 3.0,
                    CLOSE_TOL);

  // test gamma gradient
  shared_ptr<gsl_vector> phi_sum(gsl_vector_alloc(2), gsl_vector_free);
  shared_ptr<gsl_vector> nu_interaction(gsl_vector_alloc(2), gsl_vector_free);

  DocumentMapper::CumulateGammaStatistics(d.get(),
                                     nu.get(),
                                     phi_sum.get(),
                                     nu_interaction.get());

  BOOST_CHECK_CLOSE(gsl_vector_get(phi_sum.get(), 0),
                    1.75,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_vector_get(phi_sum.get(), 1),
                    1.25,
                    CLOSE_TOL);


  // 1*1/4+1/2*(1/4*1/2+3/4*9/10)+1/3*(1/4*1/2+3/4*9/10)
  BOOST_CHECK_CLOSE(gsl_vector_get(nu_interaction.get(), 0),
                    0.9166667,
                    CLOSE_TOL);

  // 1*3/4+1/2*(1/4*1/2+3/4*1/10)+1/3*(1/4*1/2+3/4*1/10)
  BOOST_CHECK_CLOSE(gsl_vector_get(nu_interaction.get(), 1),
                    0.9166667,
                    CLOSE_TOL);

  shared_ptr<gsl_vector> alpha(gsl_vector_alloc(2), gsl_vector_free);
  gsl_vector_set(alpha.get(), 0, 0.1);
  gsl_vector_set(alpha.get(), 1, 0.01);

  // test gamma lhood
  double lhood = var_dir_lhood(alpha.get(), gamma.get()) -
    var_dir_lhood(gamma.get(), gamma.get());
  // lhood += -0.9166667 * (1.0) + 1.75*(digamma(.2) - digamma(1)) +
  //                             + 1.25*(digamma(.8) - digamma(1))
  BOOST_CHECK_CLOSE(lhood + -9.6471,
                    -gamma_lhood_base(gamma.get(),
                                      alpha.get(),
                                      phi_sum.get(),
                                      nu_interaction.get()),
                    CLOSE_TOL);

  // test gamma gradient

  // 1 * (1/4*8/10-3/4*8/10) + (1/2+1/3)*((1/4*1/2*8/10+3/4*9/10*8/10) -
  // (1/4*1/2*8/10+3/4*1/10*8/10))
  BOOST_CHECK_CLOSE(gamma_gradient_interaction(gsl_vector_get(gamma.get(), 0),
                                               gsl_blas_dsum(gamma.get()),
                                               gamma.get(),
                                               nu_interaction.get(),
                                               0),
                    -0.0,
                    CLOSE_TOL);

  // (3/4*2/10-1/4*2/10) + (1/2+1/3)*((1/4*1/2*2/10+3/4*1/10*2/10) -
  // (1/4*1/2*2/10+3/4*9/10*2/10))
  BOOST_CHECK_CLOSE(gamma_gradient_interaction(gsl_vector_get(gamma.get(), 1),
                                               gsl_blas_dsum(gamma.get()),
                                               gamma.get(),
                                               nu_interaction.get(),
                                               1),
                    -5.55e-17,
                    LOOSE_TOL);
}

/*
 * This is not a unit test and should be stubbed out, but I'm relying
 * on code that doesn't really respond to that
void test_matrix_accum() {
  shared_ptr<gsl_matrix> m1(gsl_matrix_alloc(2,3), gsl_matrix_free);
  shared_ptr<gsl_matrix> m2(gsl_matrix_alloc(2,3), gsl_matrix_free);
  gsl_matrix_set_all(m1.get(), 0);
  gsl_matrix_set_all(m2.get(), 0);

  gsl_matrix_set(m1.get(), 0, 1, 1.0);
  gsl_matrix_set(m2.get(), 0, 1, 1.0);
  gsl_matrix_set(m2.get(), 1, 2, 1.0);
  gsl_matrix_set(m1.get(), 0, 0, 3.0);

  mtx_fprintf("temp.0", m1.get());
  mtx_fprintf("temp.1", m2.get());

  shared_ptr<gsl_matrix> m(gsl_matrix_alloc(2,3), gsl_matrix_free);
  gsl_matrix_set_all(m.get(), 1.0);

  DocumentMapper::AccumulateMatrixFromFiles(m.get(), "temp", "", 2);
  DocumentMapper::AccumulateMatrixFromFiles(NULL, "temp", "", 2);
  display_matrix(m.get(), "M");

  BOOST_CHECK_CLOSE(gsl_matrix_get(m.get(), 1, 2), 1.0, CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_matrix_get(m.get(), 0, 1), 2.0, CLOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_matrix_get(m.get(), 0, 0), 3.0, CLOSE_TOL);
}
*/

void test_nu_lhood_gradient() {
  nu_params p;
  int num_topics = 2;
  scoped_ptr<Document> d(load_test_doc("2 0 2 1 1 0 -1 0 1 0 -1 -1",
                                       num_topics, 0));

  // First will fill these via an algorithm, then put in set values
  // that will make verifying them easier.
  shared_ptr<gsl_matrix> transition_counts(gsl_matrix_alloc(num_topics + 1,
                                                            num_topics),
                                           gsl_matrix_free);

  // There's only one doc
  shared_ptr<gsl_matrix> nu_gamma_interaction(gsl_matrix_alloc(num_topics + 1,
                                                               num_topics),
                                              gsl_matrix_free);

  shared_ptr<gsl_vector> gamma(gsl_vector_alloc(num_topics), gsl_vector_free);

  gsl_vector_set(d->get_omega_(), 0, 1.0/1.0);
  gsl_vector_set(d->get_omega_(), 1, 1.0/2.0);
  gsl_vector_set(d->get_omega_(), 2, 1.0/3.0);


  gsl_vector_set(gamma.get(), 0, 2.0/3.0);
  gsl_vector_set(gamma.get(), 1, 4.0/3.0);

  //                            (0.25, 0.75)
  LensNode* n1 = d->get_head_(0);
  gsl_vector_set(n1->phi_.get(), 0, 0.25);
  gsl_vector_set(n1->phi_.get(), 1, 0.75);
  //       (1.0, 0.0)
  LensNode* n2 = n1->GetChild(0);
  gsl_vector_set(n2->phi_.get(), 0, 1.0);
  gsl_vector_set(n2->phi_.get(), 1, 0.0);
  //                                                 (0.5, 0.5)

  gsl_matrix_set_all(transition_counts.get(), 0.0);
  gsl_matrix_set_all(nu_gamma_interaction.get(), 0.0);
  DocumentMapper::CumulateNuStatistics(d.get(),
                                   0,
                                   gamma.get(),
                                   transition_counts.get(),
                                   nu_gamma_interaction.get());


  // Test the gamma interactions

  // 1/2*1/4*1/3 + 1/3*1/4*1/3
  BOOST_CHECK_CLOSE(gsl_matrix_get(nu_gamma_interaction.get(), 0, 0),
                    0.06944444,
                    CLOSE_TOL);

  // 1/2*1/4*2/3 + 1/3*1/4*2/3
  BOOST_CHECK_CLOSE(gsl_matrix_get(nu_gamma_interaction.get(), 0, 1),
                    0.1388889,
                    CLOSE_TOL);

  // 1/2*3/4*1/3 + 1/3*3/4*1/3
  BOOST_CHECK_CLOSE(gsl_matrix_get(nu_gamma_interaction.get(), 1, 0),
                    0.2083333,
                    CLOSE_TOL);

  // 1/2*3/4*2/3 + 1/3*3/4*2/3
  BOOST_CHECK_CLOSE(gsl_matrix_get(nu_gamma_interaction.get(), 1, 1),
                    0.4166667,
                    CLOSE_TOL);

  // 1 * 1/3
  BOOST_CHECK_CLOSE(gsl_matrix_get(nu_gamma_interaction.get(), 2, 0),
                    1.0/3.0,
                    CLOSE_TOL);

  // 1 * 2/3
  BOOST_CHECK_CLOSE(gsl_matrix_get(nu_gamma_interaction.get(), 2, 1),
                    2.0/3.0,
                    CLOSE_TOL);

  // Test the transitions

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 0, 0),
                    1.0/4.0+1.0/8.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 0, 1),
                    1.0/8.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 1, 0),
                    3.0/4.0 + 3.0/8.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 1, 1),
                    3.0/8.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 2, 0),
                    1.0/4.0,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 2, 1),
                    3.0/4.0,
                    CLOSE_TOL);


  /*
   * 0.1 0.05
   */
  shared_ptr<gsl_vector> alpha(gsl_vector_alloc(2), gsl_vector_free);
  gsl_vector_set(alpha.get(), 0, 0.1);
  gsl_vector_set(alpha.get(), 1, 0.05);

  shared_ptr<gsl_matrix> nu(gsl_matrix_alloc(3, 2), gsl_matrix_free);
  gsl_matrix_set(nu.get(), 0, 0, 0.5);
  gsl_matrix_set(nu.get(), 0, 1, 0.5);
  gsl_matrix_set(nu.get(), 1, 0, 0.9);
  gsl_matrix_set(nu.get(), 1, 1, 0.1);
  gsl_matrix_set(nu.get(), 2, 0, 0.25);
  gsl_matrix_set(nu.get(), 2, 1, 0.75);
  gsl_vector_view nu1_view = gsl_matrix_row(nu.get(), 0);
  gsl_vector_view nu2_view = gsl_matrix_row(nu.get(), 1);
  gsl_vector_view nu3_view = gsl_matrix_row(nu.get(), 2);
  gsl_vector* nu1 = &nu1_view.vector;
  gsl_vector* nu2 = &nu2_view.vector;
  gsl_vector* nu3 = &nu3_view.vector;

  p.transition_counts = transition_counts.get();
  p.nu_gamma_interaction = nu_gamma_interaction.get();
  p.alpha = alpha.get();
  p.topic = 0;

  // Test the lhood components at an arbitrary location
  //
  // testing var_dir_lhood with a scaling factor is done in the
  // document lhood testing section
  BOOST_CHECK_CLOSE(var_dir_lhood(p.alpha, nu1),
                    -0.829133508,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(var_dir_lhood(nu1, nu1),
                    .2415644750,
                    CLOSE_TOL);

  double doc_lhood_term = nu_transition_likelihood(p.transition_counts,
                                                   nu1,
                                                   0);
  BOOST_CHECK_CLOSE(nu_transition_likelihood(p.transition_counts,
                                             nu1,
                                             0),
                    nu_doc_transition_likelihood(d.get(),
                                                 nu1,
                                                 0),
                    CLOSE_TOL);

  doc_lhood_term += nu_transition_likelihood(p.transition_counts,
                                             nu2,
                                             1);
  BOOST_CHECK_CLOSE(nu_transition_likelihood(p.transition_counts,
                                             nu2,
                                             1),
                    nu_doc_transition_likelihood(d.get(),
                                                 nu2,
                                                 1),
                    CLOSE_TOL);

  doc_lhood_term += nu_transition_likelihood(p.transition_counts,
                                             nu3,
                                             2);
  BOOST_CHECK_CLOSE(nu_transition_likelihood(p.transition_counts,
                                             nu3,
                                             2),
                    nu_doc_transition_likelihood(d.get(),
                                                 nu3,
                                                 2),
                    CLOSE_TOL);

  doc_lhood_term -= gamma_nu_total_interaction(nu_gamma_interaction.get(),
                                               0,
                                               nu1);

  // 1/4*(1/2+1/3)*(2/3*1/2+4/3*1/2)/2
  BOOST_CHECK_CLOSE(gamma_nu_total_interaction(nu_gamma_interaction.get(),
                                               0,
                                               nu1),
                    0.1041667,
                    CLOSE_TOL);

  doc_lhood_term -= gamma_nu_total_interaction(nu_gamma_interaction.get(),
                                               1,
                                               nu2);

  // 3/4*(1/2+1/3)*(2/3*9/10+4/3*1/10)/2
  BOOST_CHECK_CLOSE(gamma_nu_total_interaction(nu_gamma_interaction.get(),
                                               1,
                                               nu2),
                    0.2291667,
                    CLOSE_TOL);

  doc_lhood_term -= gamma_nu_total_interaction(nu_gamma_interaction.get(),
                                               2,
                                               nu3);

  // 1/1*(1/4*2/3+3/4*4/3)/2
  BOOST_CHECK_CLOSE(gamma_nu_total_interaction(nu_gamma_interaction.get(),
                                               2,
                                               nu3),
                    0.583333,
                    CLOSE_TOL);

  // The previous tests are only for consistency, not correctness.


  // (1/4*1 + 1/4*1/2 + 1/4*1/2) * (digamma(1/2)-digamma(1))
  BOOST_CHECK_CLOSE(nu_transition_likelihood(p.transition_counts,
                                             nu1,
                                             0),
                    -0.6931472,
                    CLOSE_TOL);

  // (3/4*1+3/4*1/2)*(digamma(.9)-digamma(1.0)) +
  // 3/4*1/2 *(digamma(.1)-digamma(1.0))
  BOOST_CHECK_CLOSE(nu_transition_likelihood(p.transition_counts,
                                             nu2,
                                             1),
                    -3.892377,
                    CLOSE_TOL);

  // .25*(digamma(.25)-digamma(1.)) + .75*(digamma(.75)-digamma(1.0))
  BOOST_CHECK_CLOSE(nu_transition_likelihood(p.transition_counts,
                                             nu3,
                                             2),
                    -1.294043,
                    CLOSE_TOL);

  p.topic = 0;
  double alt_lhood = nu_lhood(nu1, &p);
  p.topic = 1;
  alt_lhood += nu_lhood(nu2, &p);
  p.topic = 2;
  alt_lhood += nu_lhood(nu3, &p);

  // This method has been removed out of the document

  /*
  BOOST_CHECK_CLOSE(doc_lhood_term + DocumentMapper::GlobalWeightTerm(nu.get(),
                                                                 alpha.get(),
                                                                 1.0,
                                                                 1.0,
                                                                 1.0,
                                                                 false),
                    -alt_lhood,
                    CLOSE_TOL);
  */

  // TODO(jbg): Add a test of self-weighting

  // nu_gradient (0,0)

  p.topic = 0;
  // (psigamma(1.0,1))*(0.15 + (1.0/4.0+1.0/8.0+1.0/8.0) - 1.0)
  double second_term = -0.5757269;
  BOOST_CHECK_CLOSE(nu_gradient_second_term(&p,
                                            nu1,
                                            gsl_blas_dsum(nu1)),
                    second_term,
                    CLOSE_TOL);

  // first_term = psigamma(0.5, 1)*(0.1+(1.0/4.0+1.0/8.0)-0.5) = -0.1233701
  double first_term = -0.1233701;

  // third_term = 1/4*(1/2+1/3)*(2/3*1/2-4/3*1/2)/2 = -0.03472222
  double third_term = -0.03472222;

  BOOST_CHECK_CLOSE(nu_gradient_interaction(gsl_vector_get(nu1, 0),
                                            gsl_blas_dsum(nu1),
                                            nu1,
                                            nu_gamma_interaction.get(),
                                            0,
                                            0),
                    third_term,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(nu_gradient_j(&p,
                                  nu1,
                                  0,
                                  second_term,
                                  gsl_blas_dsum(nu1)),
                    first_term - second_term - third_term,
                    CLOSE_TOL);

  // nu_gradient (0, 1)

  // first_term = psigamma(0.5, 1)*(0.05+(.25*0+.25*.5)-0.5) = -1.603811
  first_term = -1.603811;

  // 0.25*(1/2+1/3)*(4/3*1/2-2/3*1/2)/2
  third_term = 0.03472222;

  BOOST_CHECK_CLOSE(nu_gradient_interaction(gsl_vector_get(nu1, 0),
                                            gsl_blas_dsum(nu1),
                                            nu1,
                                            nu_gamma_interaction.get(),
                                            0,
                                            1),
                    third_term,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(nu_gradient_j(&p,
                                  nu1,
                                  1,
                                  second_term,
                                  gsl_blas_dsum(nu1)),
                    first_term - second_term - third_term,
                    CLOSE_TOL);

  // nu_gradient (2, 0)

  // first_term = psigamma(0.25, 1)*(0.1+.25-0.25) = 1.719733
  first_term = 1.719733;

  // (1/1)*(2/3*3/4-4/3*3/4)/2
  third_term =-0.25;

  p.topic = 2;
  BOOST_CHECK_CLOSE(nu_gradient_interaction(gsl_vector_get(nu3, 0),
                                            gsl_blas_dsum(nu3),
                                            nu3,
                                            nu_gamma_interaction.get(),
                                            2,
                                            0),
                    third_term,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(nu_gradient_j(&p,
                                  nu3,
                                  0,
                                  second_term,
                                  gsl_blas_dsum(nu3)),
                    first_term - second_term - third_term,
                    CLOSE_TOL);
}

void test_nu_extreme_values() {
  shared_ptr<gsl_matrix> nu(gsl_matrix_alloc(2, 6), gsl_matrix_free);
  shared_ptr<gsl_vector> beta(gsl_vector_alloc(6), gsl_vector_free);

  int num_topics = 3;
  scoped_ptr<Document> d(load_test_doc("2 0 2 1 1 1 0 4 0 -1 -1 3 2 1 4 3" +
                                       " 1 0 4 0 -1 -1 -1 -1", num_topics, 0));

  gsl_matrix_set(nu.get(), 0, 0, 0.0166817);
  gsl_matrix_set(nu.get(), 0, 1, 0.0740372);
  gsl_matrix_set(nu.get(), 0, 2, 0.0166817);
  gsl_matrix_set(nu.get(), 0, 3, 0.196154);
  gsl_matrix_set(nu.get(), 0, 4, 57.8978);
  gsl_matrix_set(nu.get(), 0, 5, 0.0166817);

  gsl_matrix_set(nu.get(), 1, 0, 0.00495806);
  gsl_matrix_set(nu.get(), 1, 1, 0.00507185);
  gsl_matrix_set(nu.get(), 1, 2, 0.00495806);
  gsl_matrix_set(nu.get(), 1, 3, 0.00509654);
  gsl_matrix_set(nu.get(), 1, 4, 1.11406e+17);
  gsl_matrix_set(nu.get(), 1, 5, 0.00495806);

  gsl_vector_set_all(beta.get(), 0.1/6.0);

  gsl_vector_view nu1_view = gsl_matrix_row(nu.get(), 0);
  gsl_vector_view nu2_view = gsl_matrix_row(nu.get(), 1);
  gsl_vector* nu1 = &nu1_view.vector;
  gsl_vector* nu2 = &nu2_view.vector;

  BOOST_CHECK_CLOSE(gsl_sf_psi(gsl_blas_dsum(nu1)),
                    4.055582,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_sf_psi(gsl_vector_get(nu1, 4)),
                    4.055582,
                    LOOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_sf_psi(gsl_vector_get(nu1, 3)),
                    -5.39196,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_sf_lngamma(gsl_blas_dsum(beta.get())),
                    2.252713,
                    CLOSE_TOL);

  display_vector(beta.get(), "beta");
  BOOST_CHECK_CLOSE(var_dir_lhood(beta.get(), nu1),
                    195.1893,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(var_dir_lhood(beta.get(), nu2),
                    1155.342,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(var_dir_lhood(nu1, nu1),
                    199.3280,
                    CLOSE_TOL);

  BOOST_CHECK_CLOSE(var_dir_lhood(nu2, nu2),
                    235.0115,
                    CLOSE_TOL);

  // Now we'll test a vanilla nu update from a setting of phi that
  // caused problems
  shared_ptr<gsl_vector> gamma(gsl_vector_alloc(num_topics), gsl_vector_free);
  gsl_vector_set_all(gamma.get(), 1.0);
  nu.reset(gsl_matrix_alloc(num_topics + 1, num_topics));
  gsl_matrix_set_all(nu.get(), 1.0/(float)num_topics);

  int num_terms = 8;
  int num_roles = 5;
  double alpha = 0.1;
  shared_ptr<gsl_matrix> tau(gsl_matrix_alloc(num_topics, num_terms));
  shared_ptr<gsl_matrix> rho(gsl_matrix_alloc(num_roles, num_terms));
  gsl_matrix_set_all(tau.get(), 0.0);
  gsl_matrix_set_all(rho.get(), 0.0);

  /*
  gsl_vector_set(d->get_word_(0)->phi_.get(), 0, 0.583396);
  gsl_vector_set(d->get_word_(0)->phi_.get(), 1, 0.200404);
  gsl_vector_set(d->get_word_(0)->phi_.get(), 2, 0.2162);
  gsl_vector_set(d->get_word_(1)->phi_.get(), 0, 0.186152);
  gsl_vector_set(d->get_word_(1)->phi_.get(), 1, 0.697103);
  gsl_vector_set(d->get_word_(1)->phi_.get(), 2, 0.116745);
  gsl_vector_set(d->get_word_(2)->phi_.get(), 0, 0.143174);
  gsl_vector_set(d->get_word_(2)->phi_.get(), 1, 0.417413);
  gsl_vector_set(d->get_word_(2)->phi_.get(), 2, 0.439413);
  gsl_vector_set(d->get_word_(3)->phi_.get(), 0, 0.457082);
  gsl_vector_set(d->get_word_(3)->phi_.get(), 1, 0.37833);
  gsl_vector_set(d->get_word_(3)->phi_.get(), 2, 0.164589);
  gsl_vector_set(d->get_word_(4)->phi_.get(), 0, 0.336124);
  gsl_vector_set(d->get_word_(4)->phi_.get(), 1, 0.451827);
  gsl_vector_set(d->get_word_(4)->phi_.get(), 2, 0.212049);
  gsl_vector_set(d->get_word_(5)->phi_.get(), 0, 0.143174);
  gsl_vector_set(d->get_word_(5)->phi_.get(), 1, 0.417413);
  gsl_vector_set(d->get_word_(5)->phi_.get(), 2, 0.439413);
  */

  /*
  gsl_vector_set(d->get_xi_(), 0, 1.0/8.0);
  gsl_vector_set(d->get_xi_(), 1, 1.0/8.0);
  gsl_vector_set(d->get_xi_(), 2, 1.0/8.0);
  gsl_vector_set(d->get_xi_(), 3, 1.0/8.0);
  gsl_vector_set(d->get_xi_(), 4, 1.0/8.0);
  gsl_vector_set(d->get_xi_(), 5, 1.0/8.0);
  */

  gsl_vector_set_all(d->get_omega_(), 0.0);

  // First will fill these via an algorithm, then put in set values
  // that will make verifying them easier.
  shared_ptr<gsl_matrix> transition_counts(gsl_matrix_alloc(num_topics + 1,
                                                            num_topics),
                                           gsl_matrix_free);

  // There's only one doc
  shared_ptr<gsl_matrix> nu_gamma_interaction(gsl_matrix_alloc(num_topics + 1,
                                                               num_topics),
                                              gsl_matrix_free);

  gsl_matrix_set_all(transition_counts.get(), 0.0);
  gsl_matrix_set_all(nu_gamma_interaction.get(), 0.0);
  gsl_vector_set_all(beta.get(), alpha/(float)num_topics);
  DocumentMapper::CumulateNuStatistics(d.get(),
                                   0,
                                   gamma.get(),
                                   transition_counts.get(),
                                   nu_gamma_interaction.get());

  /*
  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 0, 0),
                    .603,
                    LOOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 0, 1),
                    1.051,
                    LOOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 2, 1),
                    .444,
                    LOOSE_TOL);

  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 3, 1),
                    0.200,
                    LOOSE_TOL);
  */

  // 1/3*1/3*5
  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 2, 0),
                    5.0/9.0,
                    LOOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 2, 1),
                    5.0/9.0,
                    LOOSE_TOL);
  BOOST_CHECK_CLOSE(gsl_matrix_get(transition_counts.get(), 2, 2),
                    5.0/9.0,
                    LOOSE_TOL);


  double init_doc = DocumentMapper::Likelihood(d.get(),
                                          gamma.get(),
                                          alpha,
                                          beta.get(),
                                          nu.get(),
                                          num_topics,
                                          num_terms,
                                          tau.get(),
                                          NULL);
  double init_global = SyntopReducer::GlobalWeightTerm(nu.get(),
                                                   beta.get(),
                                                   alpha,
                                                   alpha,
                                                   // 1.0,
                                                   true);

  cout << "BEFORE: doc=" << init_doc
       << " global=" << init_global << endl;

  for (int j = 0; j < num_topics; j++) {
    gsl_matrix_set(nu.get(), 0, j, gsl_vector_get(beta.get(), j) +
                   gsl_matrix_get(transition_counts.get(), 0, j));
  }

  double after1_doc = DocumentMapper::Likelihood(d.get(),
                                            gamma.get(),
                                            alpha,
                                            beta.get(),
                                            nu.get(),
                                            num_topics,
                                            num_terms,
                                            tau.get(),
                                            NULL);
  double after1_global = SyntopReducer::GlobalWeightTerm(nu.get(),
                                                     beta.get(),
                                                     alpha,
                                                     alpha,
                                                         // 1.0,
                                                     true);

  BOOST_CHECK(init_global + init_doc < after1_doc + after1_global);
  cout << "AFTER 1: doc=" << after1_doc
       << " global=" << after1_global << endl;

  for (int j = 0; j < num_topics; j++) {
    gsl_matrix_set(nu.get(), 1, j, gsl_vector_get(beta.get(), j) +
                   gsl_matrix_get(transition_counts.get(), 1, j));
  }

  double after2_doc = DocumentMapper::Likelihood(d.get(),
                                            gamma.get(),
                                            alpha,
                                            beta.get(),
                                            nu.get(),
                                            num_topics,
                                            num_terms,
                                            tau.get(),
                                            NULL);
  double after2_global = SyntopReducer::GlobalWeightTerm(nu.get(),
                                                     beta.get(),
                                                     alpha,
                                                     alpha,
                                                         // 1.0,
                                                     true);

  BOOST_CHECK(after1_doc + after1_global < after2_doc + after2_global);

  cout << "AFTER 2: doc=" << after2_doc
       << " global=" << after2_global << endl;


  for (int j = 0; j < num_topics; j++) {
    gsl_matrix_set(nu.get(), 2, j, gsl_vector_get(beta.get(), j) +
                   gsl_matrix_get(transition_counts.get(), 2, j));
  }

  double after3_doc = DocumentMapper::Likelihood(d.get(),
                                            gamma.get(),
                                            alpha,
                                            beta.get(),
                                            nu.get(),
                                            num_topics,
                                            num_terms,
                                            tau.get(),
                                            NULL);
  double after3_global = SyntopReducer::GlobalWeightTerm(nu.get(),
                                                     beta.get(),
                                                     alpha,
                                                     alpha,
                                                         // 1.0,
                                                     true);

  BOOST_CHECK(after2_doc + after2_global < after3_doc + after3_global);

  cout << "AFTER 3: doc=" << after3_doc
       << " global=" << after3_global << endl;
}

test_suite* init_unit_test_suite(int, char* []) {
  test_suite* test= BOOST_TEST_SUITE("Testing DocumentMapper");
  test->add(BOOST_TEST_CASE(&test_nu_lhood_gradient),    0);
  test->add(BOOST_TEST_CASE(&test_document_creation),    0);
  test->add(BOOST_TEST_CASE(&test_gamma_lhood_gradient), 0);
  test->add(BOOST_TEST_CASE(&test_phi_setting),          0);
  test->add(BOOST_TEST_CASE(&test_lhood),                0);
  test->add(BOOST_TEST_CASE(&test_log_functions),        0);
  test->add(BOOST_TEST_CASE(&test_tau_rho_setting),      0);
  test->add(BOOST_TEST_CASE(&test_stick_breaking),       0);
  // test->add( BOOST_TEST_CASE( &test_matrix_accum ),         0);
  test->add(BOOST_TEST_CASE(&test_nu_extreme_values),    0);
  // test->add( BOOST_TEST_CASE( &test_beta_initialization ),  0);
  return test;
}
