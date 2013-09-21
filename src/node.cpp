/*
 * Copyright 2009 Jordan Boyd-Graber
 */

# include "node.h"
/*
 * Swap
 */
void LensNode::ReplacePhiByUnnormalizedLogSpace(gsl_vector* new_phi) {
  // display_vector(new_phi, "old phi");
  log_normalize(new_phi);
  // cout << "SUM: " << sum << endl;
  vexp(new_phi, phi_.get());
  // display_vector(phi_.get(), "phi is\n");
}

/*
struct SamplingNode: public Node {
        int synset_;
        double LogProbability(Model* model, int candidate_synset, int parent_synset);
        void Sample(Model* model, int parent_synset, bool init);
        std::string StateString();
};
*/

int LensNode::get_parent_index_() {
        return this->parent_index_;
}

void LensNode::PrintTree() {
  cout << this->term_ << " ";
  for (int i = 0; i < this->number_children_; i++) {
    cout << "(";
    this->children_[i]->PrintTree();
    cout << ")";
  }
}

int LensNode::FillWordArray(LensNode** words, int parent,
                            int position, LensNode* current) {
  // cout << this->word_ << " " << this->number_children_ << endl;
  this->children_positions_.resize(this->number_children_);
  if (words) {
    current->SetParentIndex(parent);
    words[position] = current;
    // cout << "Setting position " << position << " " << current << endl;
  }
  int child_position = position + 1;
  for (int i = 0; i < this->number_children_; i++) {
    this->children_positions_[i] = child_position;
    child_position =
      this->children_[i]->FillWordArray(words,
                                        position,
                                        child_position,
                                        this->children_[i].get());
  }
  return child_position;
}

void LensNode::SetParentIndex(int parent_index) {
        this->parent_index_ = parent_index;
}

void LensNode::SetNumberChildren(int number_children) {
  this->number_children_ = number_children;
  this->children_.reset(new boost::scoped_ptr<LensNode>[number_children]);
  for (int i = 0; i < number_children; i++) {
    this->children_[i].reset(new LensNode());
  }
}

/*
std::string SamplingNode::StateString() {
        std::ostringstream synset;
        synset << this->synset_ << " ";
        string state = synset.str();

        for(int i=0; i<this->number_children_; i++) {
                state += this->children_[i]->StateString();
        }
        return state;
}
*/

LensNode* LensNode::AddChild(int word, int relation, int number_children) {
  // Find the next free slot
  int child_id=-1;
  while (this->children_[++child_id].get()->term_ >= 0);
  assert(child_id < number_children_);

  LensNode* child = this->children_[child_id].get();
  child->relation_ = relation;

  child->term_ = word;
  child->SetNumberChildren(number_children);
  return child;
}

/*
double SamplingNode::LogProbability(Model* model, int candidate_synset, int parent_synset) {
        double p = 0.0;
        p += model->sucessor_log_probability_(parent_synset, this->relation_, candidate_synset);
        p += model->word_log_probability_(candidate_synset, this->word_);
        for(int i=0; i<this->number_children_; i++) {
                p += model->sucessor_log_probability_(candidate_synset,
                                this->children_[i]->relation_,
                                this->children_[i]->synset_);
        }

        return p;
}
*/

LensNode::LensNode() {
        this->term_ = -1;
        this->relation_ = -1;
        this->number_children_ = 0;
        this->parent_index_ = -1;
}

LensNode* LensNode::GetChild(int child) {
        return this->children_[child].get();
}

/*
void SamplingNode::Sample(Model* model, int parent_synset, bool init) {
        vector<double> probs;

        if(!init) {
                // On the first round, all counts are zero, so we can't decrement.
                model->UpdateSucessorCounts(parent_synset, this, -1);
                model->UpdateSynsetCounts(this->synset_, this->word_, -1);
        }

        // If we are initializing and we have a state, that means that we resumed from a checkpointed state
        if(init && this->synset_!=-1) {
                for(int i=0; i<model->get_num_synsets_(); i++) {
                        probs.push_back(this->LogProbability(model, i, parent_synset));
                }
                int new_synset = sample_from_vector(probs, false);
                if(this->synset_!=new_synset) {
                        cout << "Changing from synset " << this->synset_ << " to " << new_synset << endl;
                }
                this->synset_ = new_synset;
        }

        model->UpdateSucessorCounts(parent_synset, this, +1);
        model->UpdateSynsetCounts(this->synset_, this->word_, +1);

        for(int i=0; i<this->number_children_; i++) {
                this->children_[i]->Sample(model, this->synset_, init);
        }
}
*/
