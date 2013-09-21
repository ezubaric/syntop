/*
 *  node.h
 *
 *
 *  Created by Jordan Boyd-Graber on 27.09.07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SYNTOP_NODE_INCLUDED
#define SYNTOP_NODE_INCLUDED

#include <vector>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <gsl/gsl_vector.h>
#include "topicmod/lib/prob/logspace.h"
#include "vectorops.h"
#include "util.h"

using lib_prob::log_sum;
using lib_prob::log_normalize;

/*
 * This is templated because Nodes need to know the type of their
 * children.  I.e. LensNode needs to know its children is also
 * LensNode.
 *
 * TODO: This was a struct, so everything is public, that should
 * probably change.
 */
class LensNode {
 public:
  LensNode();
  boost::shared_ptr< gsl_vector > phi_;
  boost::shared_ptr< gsl_vector > temp_phi_;
  void ReplacePhiByUnnormalizedLogSpace(gsl_vector* new_phi);
  void SetNumberChildren(int number_children);
  void SetParentIndex(int parent_index);
  int get_parent_index_();
  LensNode* AddChild(int word, int relation, int number_children);
  LensNode* GetChild(int child);
  int FillWordArray(LensNode** words, int parent, int position,
                    LensNode* current);
  void PrintTree();
  int number_children_;
  int term_;
  int relation_;
  int parent_index_;
  // This is redundant, but sometimes we need the indices and sometimes
  // we need the child itself.
  vector<int> children_positions_;
  boost::scoped_array< boost::scoped_ptr<LensNode> > children_;
};


#endif
