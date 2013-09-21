/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <assert.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::ofstream;

typedef std::map<string, double> StringMap;

#define ERROR_TOLERANCE 1e-2

string dtos(double d);

bool my_isnan(double x);

int rand_int(int max);

int sample_from_vector(vector<double> probabilities, bool normalized);

void matrix_rand_init(gsl_matrix* m, double min, double max, bool prob);

void vector_rand_init(gsl_vector* v, double min, double max, bool prob);

void display_vector(const gsl_vector* v, const char* name,
                    std::ostream& place = cout);

// This should be templated if I ever need anything other than an int
void display_vector(std::vector<int> v, const char* name);

double safe_exp(double val);

void display_matrix(const gsl_matrix* m, const char* name,
                    std::ostream& place = cout);

bool check_valid_double(double value, string name);

double check_increase(string variable,
                      ofstream& outfile,
                      double new_value,
                      double old_value,
                      string prefix = "lhood",
                      bool always_write = false);

#endif
