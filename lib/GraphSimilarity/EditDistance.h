#pragma once
#include "Base/Logger.h"
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "DifferenceVector.h"
#include "GraphSimilarity/GSSEntry.h"
#include "GraphSimilaritySearch.h"
#include "State.h"
#include "WeisfeilerLehman.h"

namespace GraphLib::GraphSimilarity {
bool verbosity = true;
static int32_t LOG_EVERY = 50000;

class GraphEditDistanceSolver {
protected:
  GSSEntry *G1, *G2;
  int threshold = -1, current_best = 1e9;
  int64_t num_nodes = 0;
  std::vector<int> current_best_mapping;
  std::vector<std::vector<int>> branch_distances, branch_scores;
  std::priority_queue<State *, std::vector<State *>, StateComparator> queue;
  DifferenceVector *vlabel_diff = new DifferenceVector(100),
                   *elabel_diff = new DifferenceVector(20);
  std::vector<int> matching_order, inv_matching_order;
  ResultLogger log;

public:
  void InitializeSolver(GSSEntry *G1_, GSSEntry *G2_, int threshold_ = -1) {
    log.clear();
    current_best = 1e9;
    // ensure V(G1) <= V(G2)
    if (G1_->GetNumVertices() > G2_->GetNumVertices()) {
      this->G1 = G2_;
      this->G2 = G1_;
    } else {
      this->G1 = G1_;
      this->G2 = G2_;
    }
    this->threshold = threshold_;
  }
  /*
   * Function for filtering
   * used for GED verification tasks
   */
  int NaiveCountBound() {
    int bound = abs(G1->GetNumVertices() - G2->GetNumVertices()) +
                abs(G1->GetNumEdges() - G2->GetNumEdges());
    return bound;
  }

  int LabelSetDifferenceBound() {
    vlabel_diff->reset();
    elabel_diff->reset();
    for (auto &[l, t] : G2->GetVertexLabelFrequency())
      vlabel_diff->update(l, t);
    for (auto &[l, t] : G1->GetVertexLabelFrequency())
      vlabel_diff->update(l, -t);
    for (auto &[l, t] : G1->GetEdgeLabelFrequency())
      elabel_diff->update(l, t);
    for (auto &[l, t] : G2->GetEdgeLabelFrequency())
      elabel_diff->update(l, -t);
    int cost = 0;
    cost = vlabel_diff->GetDifference() + elabel_diff->GetDifference();
    return cost;
  }

  int DegreeSequenceBound() {
    int pos = 0, neg = 0;
    auto &g1_deg = G1->GetDegreeSequence();
    auto &g2_deg = G2->GetDegreeSequence();
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      int a = (i >= G2->GetNumVertices()) ? 0 : g1_deg[i];
      int b = g2_deg[i];
      if (a > b)
        pos += (a - b);
      if (a < b)
        neg += (b - a);
    }
    return (pos + 1) / 2 + (neg + 1) / 2;
  }

  int BranchBound() {
    branch_distances.clear();
    branch_distances.resize(G2->GetNumVertices(),
                            std::vector<int>(G2->GetNumVertices(), 0));
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      for (int j = 0; j < G2->GetNumVertices(); j++) {
        branch_distances[i][j] =
            BranchEditDistance(G1->GetBranch(i), G2->GetBranch(j));
      }
    }
    for (int j = 0; j < G2->GetNumVertices(); j++) {
      for (int i = G1->GetNumVertices(); i < G2->GetNumVertices(); i++) {
        branch_distances[i][j] = BranchEditDistanceFromNull(G2->GetBranch(j));
      }
    }

    Hungarian hungarian(branch_distances);
    hungarian.Solve();
    int cost = (hungarian.GetTotalCost() + 1) / 2;
    return cost;
  }

  bool GEDVerificiationFiltering() {
    if (NaiveCountBound() > threshold)
      return false;
    if (LabelSetDifferenceBound() > threshold)
      return false;
    if (DegreeSequenceBound() > threshold)
      return false;
    if (BranchBound() > threshold)
      return false;
    return true;
  }

  /*
   * Functions for GED Computation / Verification
   */

  int GetCurrentBestGED() const { return current_best; }
  ResultLogger GetLog() { return log; }

  void ComputeMatchingOrder() {
    int N = G1->GetNumVertices();
    std::vector<int> T(N, 0), w(N, 0);
    auto vlabel_freq = G1->GetVertexLabelFrequency();
    auto elabel_freq = G1->GetEdgeLabelFrequency();
    for (int i = 0; i < N; i++) {
      w[i] -= 2 * vlabel_freq[G1->GetVertexLabel(i)];
      for (int x : G1->GetNeighbors(i)) {
        w[i] -= elabel_freq[G1->GetEdgeLabel(i, x)];
      }
    }
    std::fill(T.begin(), T.end(), 0);
    int max_idx = std::max_element(w.begin(), w.end()) - w.begin();
    std::priority_queue<std::pair<int, int>> weighted_queue;
    weighted_queue.push({w[max_idx], max_idx});
    T[max_idx] = 1;
    while (!weighted_queue.empty()) {
      int u = weighted_queue.top().second;
      weighted_queue.pop();
      matching_order.push_back(u);
      for (int x : G1->GetNeighbors(u)) {
        if (T[x] == 1)
          continue;
        T[x] = 1;
        weighted_queue.push({w[x], x});
      }
    }
    for (int i = 0; i < N; i++) {
      if (T[i] != 1) {
        matching_order.push_back(i);
      }
    }
    inv_matching_order.resize(G1->GetNumVertices(), -1);
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      inv_matching_order[matching_order[i]] = i;
    }
  }

  virtual int GED() { return 0; };
};
} // namespace GraphLib::GraphSimilarity