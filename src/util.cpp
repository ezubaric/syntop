/*
 * Copyright 2009 Jordan Boyd-Graber
 */
#include "util.h"

using std::numeric_limits;
using std::endl;

bool my_isnan(double x)
{ return x != x; }

/*
 * given log(a) and log(b), return log(a+b)
 *
 */

int rand_int(int max) {
  return (int)(((double)rand() / ((double)(RAND_MAX)+(double)(1.0)))*max);
}

double rand_double(double low, double high) {
  return (double)rand() /
    (((double)(RAND_MAX)+(double)(1.0))/(high-low)) + low;
}

/*
 * Warning: this is not normalized.
 */
void vector_rand_init(gsl_vector* v, double min, double max, bool prob) {
  double sum = 0.0;
  gsl_vector_set_all(v, 1.0 / ((double)v->size));

  for (unsigned int i = 0; i < v->size; i++) {
    double val = rand_double(min, max);
    sum += val;
    gsl_vector_set(v, i, val);
  }

  gsl_vector_scale(v, 1.0 / sum);
}

void matrix_rand_init(gsl_matrix* m, double min, double max, bool prob) {
  for (unsigned int i = 0; i < m->size1; i++) {
    double row_sum = 0.0;
    for (unsigned int j = 0; j < m->size2; j++) {
      double new_val = rand_double(min, max);
      row_sum += new_val;
      gsl_matrix_set(m, i, j, new_val);
    }
    gsl_vector_view row = gsl_matrix_row(m, i);
    gsl_vector_scale(&row.vector, 1.0 / row_sum);
  }
}

/*

This needs to be moved to common utils

int sample_from_vector(vector<double> probabilities, bool normalized) {
  double normalizer = safe_log(0.0);
  unsigned int i;

  if(!normalized) {
    for(i=0; i<probabilities.size(); i++) {
      normalizer = log_sum(normalizer, probabilities[i]);
    }
  }

  double cutoff = (double)rand() / ((double)RAND_MAX + 1.0);
  double sum = 0;
  double val = 0;

  for(i=0; i<probabilities.size(); i++) {
    if(normalized) {
      val = exp(probabilities[i]);
    } else {
      val = exp(probabilities[i] - normalizer);
    }
    sum += val;
    if(sum >= cutoff)
      break;
  }

  if(i>=probabilities.size()) {
    i--;
  }
  return i;
}
*/

bool check_valid_double(double value, string name) {
  if (isnan(value) ||
      value == numeric_limits<double>::infinity() ||
      value == -numeric_limits<double>::infinity()) {
    cout << name << " is not a valid double:" << value << endl;
    return false;
  } else {
    return true;
  }
}

double safe_exp(double val) {
  if (val > 250.0) {
    // cout << "@";
    return exp(250.0);
  } else {
    return exp(val);
  }
}

double check_increase(string variable,
                      ofstream& outfile,
                      double new_value,
                      double old_value,
                      string prefix,
                      bool always_write) {
  cout << "[" << prefix << "]\t";
  double diff =  old_value - new_value;
  if (new_value < old_value - ERROR_TOLERANCE) {
    cout << "Increase didn't happen for " << variable
         << ", actually fell by " << diff << " from "
         << old_value << endl;
  } else {
    cout << "Value increased by " << -diff << " to "
         << new_value << " for " << variable << endl;
  }

  if (always_write) {
    outfile << "     " << variable << " " << new_value << endl;
  } else if (new_value < old_value - ERROR_TOLERANCE) {
    outfile << " *** " << variable << " " << old_value << "->" << new_value
            << endl;
  }
  outfile.flush();
  assert(new_value >= old_value - ERROR_TOLERANCE);
  return new_value;
}

void display_vector(std::vector<int> v, const char* name) {
  cout << name << " = <";
  for (unsigned int i = 0; i < v.size(); i++) {
    cout << v[i] << ", ";
  }
  cout << ">    (" << v.size() << ")" << endl;
}

void display_vector(const gsl_vector* v, const char* name,
                    std::ostream& place) {
  place << name << " = <";
  double sum = 0.0;
  for (unsigned int i = 0; i < v->size; i++) {
    place << gsl_vector_get(v, i);
    sum += gsl_vector_get(v, i);
    if (i < v->size - 1) {
      place << ", ";
    }
  }
  place << ">    (" << v->size << ", sum=" << sum << ")" << endl;
}

void display_matrix(const gsl_matrix* m, const char* name,
                    std::ostream& place) {
  place << name << "\t = |";
  for (unsigned int i = 0; i < m->size1; i++) {
    if (i != 0) {
      place << "\t   |";
    }
    for (unsigned int j = 0; j < m->size2; j++) {
      place << gsl_matrix_get(m, i, j) << "\t";
    }
    place << "|" << endl;
  }
  place << "                                SIZE: " <<
    m->size1 << " x " << m->size2 << endl;
}
