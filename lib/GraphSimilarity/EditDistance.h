#pragma once
#include "Base/Logger.h"
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "DifferenceVector.h"
#include "GraphSimilarity/GSSEntry.h"
#include "GraphSimilaritySearch.h"
#include "State.h"
#include "WeisfeilerLehman.h"
#include "GraphSimilarity/GraphColoring/GWL.h"
#include "Base/DynamicHungarian.h"

namespace GraphLib::GraphSimilarity {
bool verbosity = true;
static int32_t LOG_EVERY = 50000;

class GraphEditDistanceSolver {
protected:
  GSSEntry *G1, *G2, *combined;
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

  int GetChildEditCost(State *parent_state, int u, int v) {
    int cost = parent_state->cost;
    int u_label = G1->GetVertexLabel(u), v_label = G2->GetVertexLabel(v);
    // update editorial cost
    if (u_label != v_label) {
      cost++;
    }
    int num_u_edges = 0;
    for (int u_nbr : G1->GetNeighbors(u)) {
      if (parent_state->mapping[u_nbr] != -1)
        num_u_edges++;
    }
    int ec = num_u_edges;
    for (int vprime : G2->GetNeighbors(v)) {
      int uprime = parent_state->inverse_mapping[vprime];
      if (uprime == -1)
        continue;
      ec++;
      int l1 = G1->GetEdgeLabel(u, uprime);
      int l2 = G2->GetEdgeLabel(v, vprime);
      if (l1 == -1)
        continue;
      if (l1 == l2)
        ec -= 2;
      else
        ec--;
    }
    cost += ec;
    return cost;
  }

  void PrepareGED(GSSEntry *combined_) {
    NumG1Vertices = G1->GetNumVertices();
    NumG2Vertices = G2->GetNumVertices();
    matching_order.clear();
    ComputeMatchingOrder();
    current_best_mapping.clear();
    current_best_mapping.resize(G1->GetNumVertices(), -1);
    num_nodes = 0;
    current_best = 1e9;
    combined = combined_;
  }

  virtual int GED() { return 0; };

  void ComputeBranchDistanceMatrix(
      State *state, std::vector<std::vector<int>> &branch_distance_matrix,
      std::vector<int> &rem_left, std::vector<int> &rem_right) {
    DifferenceVector diff;
    diff.init(20);
    for (int u = 0; u < G1->GetNumVertices(); u++) {
      if (state->mapping[u] == -1)
        rem_left.emplace_back(u);
    }
    for (int v = 0; v < G2->GetNumVertices(); v++) {
      if (state->inverse_mapping[v] == -1)
        rem_right.emplace_back(v);
    }

    for (int v_idx = 0; v_idx < rem_right.size(); v_idx++) {
      int v = rem_right[v_idx];
      auto &v_nbrs = G2->GetNeighbors(v);
      int u_idx = 0;
      for (u_idx = 0; u_idx < rem_left.size(); u_idx++) {
        int u = rem_left[u_idx];
        auto &u_nbrs = G1->GetNeighbors(u);
        diff.reset();
        if (G1->GetVertexLabel(u) != G2->GetVertexLabel(v)) {
          branch_distance_matrix[u_idx][v_idx] += 2;
        }
        for (int l = 0; l < u_nbrs.size(); l++) {
          int u_nbr = u_nbrs[l];
          int u_el = G1->GetEdgeLabel(u, u_nbr);
          if (state->mapping[u_nbr] == -1) {
            diff.update(u_el, 1);
          } else {
            int v_mapping_el = G2->GetEdgeLabel(v, state->mapping[u_nbr]);
            if (v_mapping_el != u_el) {
              branch_distance_matrix[u_idx][v_idx] += 2;
            }
          }
        }
        for (int r = 0; r < v_nbrs.size(); r++) {
          int v_nbr = v_nbrs[r];
          int v_el = G2->GetEdgeLabel(v, v_nbr);
          if (state->inverse_mapping[v_nbr] == -1) {
            diff.update(v_el, -1);
          } else {
            int u_mapping_el =
                G1->GetEdgeLabel(u, state->inverse_mapping[v_nbr]);
            if (u_mapping_el == -1) {
              branch_distance_matrix[u_idx][v_idx] += 2;
            }
          }
        }
        int inner_distance = diff.GetDifference();
        branch_distance_matrix[u_idx][v_idx] += inner_distance;
      }
      int from_null = BranchEditDistanceFromNull(G2->GetBranch(v));
      for (int v_nbr : G2->GetNeighbors(v)) {
        if (state->inverse_mapping[v_nbr] != -1) {
          from_null++;
        }
      }
      for (; u_idx < branch_distance_matrix.size(); u_idx++)
        branch_distance_matrix[u_idx][v_idx] = from_null;
    }
  }

  int ComputeDistance(std::vector<int> &mapping, bool verbose = false) {
    std::vector<int> inverse_mapping(G2->GetNumVertices(), -1);
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      inverse_mapping[mapping[i]] = i;
    }
    return ComputeDistance(mapping, inverse_mapping, verbose);
  }

  int ComputeDistance(std::vector<int> &mapping,
                      std::vector<int> &inverse_mapping, bool verbose = false) {
    int cost = 0;
    // vertex re-labeling cost
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      int l1 = G1->GetVertexLabel(i);
      int l2 = G2->GetVertexLabel(mapping[i]);
      if (l1 != l2) {
        if (verbose)
          printf("Vertex %d(%d)-%d(%d) re-labeling cost!\n", i, l1, mapping[i],
                 l2);
        cost++;
      }
    }
    if (verbose)
      printf("#Missing Vertices: %d\n",
             (G2->GetNumVertices() - G1->GetNumVertices()));
    cost += (G2->GetNumVertices() - G1->GetNumVertices());
    for (auto &[u, v] : G1->GetEdges()) {
      int fu = mapping[u], fv = mapping[v];
      int l1 = G1->GetEdgeLabel(u, v);
      int l2 = G2->GetEdgeLabel(fu, fv);
      if (l1 != l2) {
        if (verbose)
          printf("Edge (%d, %d)(%d)-(%d, %d)(%d) re-labeling cost!\n", u, v, l1,
                 fu, fv, l2);
        cost++;
      }
    }
    for (auto &[u, v] : G2->GetEdges()) {
      int inv_u = inverse_mapping[u], inv_v = inverse_mapping[v];
      if (inv_u == -1 || inv_v == -1) {
        if (verbose)
          printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n", u, v,
                 inv_u, inv_v);
        cost++;
      } else {
        int l = G1->GetEdgeLabel(inv_u, inv_v);
        if (l == -1) {
          if (verbose)
            printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n", u,
                   v, inv_u, inv_v);
          cost++;
        }
      }
    }
    if (verbose)
      printf("Total ED Cost: %d\n", cost);
    current_best = std::min(cost, current_best);
    return cost;
  }

  int ComputeDistance(State *state) {
    int cost = 0;
    // vertex re-labeling cost
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      int l1 = G1->GetVertexLabel(i);
      int l2 = G2->GetVertexLabel(state->mapping[i]);
      if (l1 != l2) {
        cost++;
      }
    }
    cost += (G2->GetNumVertices() - G1->GetNumVertices());
    for (auto &[u, v] : G1->GetEdges()) {
      int fu = state->mapping[u], fv = state->mapping[v];
      int l1 = G1->GetEdgeLabel(u, v);
      int l2 = G2->GetEdgeLabel(fu, fv);
      if (l1 != l2) {
        cost++;
      }
    }
    for (auto &[u, v] : G2->GetEdges()) {
      int inv_u = state->inverse_mapping[u], inv_v = state->inverse_mapping[v];
      if (inv_u == -1 || inv_v == -1) {
        cost++;
      } else {
        int l = G1->GetEdgeLabel(inv_u, inv_v);
        if (l == -1) {
          cost++;
        }
      }
    }
    current_best = std::min(cost, current_best);
    return cost;
  }
};
} // namespace GraphLib::GraphSimilarity