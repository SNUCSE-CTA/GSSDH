#pragma once
#include "Base/DynamicHungarian.h"
#include "Base/Logger.h"
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "DifferenceVector.h"
#include "GraphSimilarity/GSSEntry.h"
#include "GraphSimilarity/GraphColoring/GWL.h"
#include "GraphSimilaritySearch.h"
#include "State.h"
#include "WeisfeilerLehman.h"

namespace GraphLib::GraphSimilarity {
bool verbosity = true;
static int32_t LOG_EVERY = 50000;

class GraphEditDistanceSolver {
protected:
  GSSEntry *G1, *G2, *combined;
  // int threshold = -1, current_best = 1e9;
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
    // this->threshold = threshold_;
    this->threshold = threshold_;
  }
  /*
   * Function for filtering
   * used for GED verification tasks
   */

  using ui = unsigned int;
  ui size_based_bound() {

    ui r1 = G1->GetNumVertices() > G2->GetNumVertices()
                ? G1->GetNumVertices() - G2->GetNumVertices()
                : G2->GetNumVertices() - G1->GetNumVertices();
    ui r2 = G1->GetNumEdges() > G2->GetNumEdges()
                ? G1->GetNumEdges() - G2->GetNumEdges()
                : G2->GetNumEdges() - G1->GetNumEdges();
    return r1 + r2;
  }

  int *vlabel_cnt = new int[300]; // fix array size
  int *elabel_cnt = new int[300];
  int *degree_q = new int[300];
  int *degree_g = new int[300];
  int *tmp = new int[300];

  ui ged_lower_bound_filter() {
    ui verify_upper_bound = threshold;
    ui n = G1->GetNumVertices();
    ui m = 2 * G1->GetNumEdges();
    ui g_n = G2->GetNumVertices();
    ui g_m = 2 * G2->GetNumEdges();
    ui lb = size_based_bound();
    if (lb > verify_upper_bound)
      return lb;

    lb = (G1->GetNumVertices() > G2->GetNumVertices() ? G1->GetNumVertices()
                                                      : G2->GetNumVertices());
    for (ui i = 0; i < G1->GetNumVertices(); i++) {
      ++vlabel_cnt[G1->GetVertexLabel(i)];
    }
    for (ui i = 0; i < G2->GetNumVertices(); i++) {
      ui vl = G2->GetVertexLabel(i);
      if (vlabel_cnt[vl] > 0) {
        --vlabel_cnt[vl];
        --lb;
      }
    }
    for (ui i = 0; i < G1->GetNumVertices(); i++)
      vlabel_cnt[G1->GetVertexLabel(i)] = 0;
    if (lb > verify_upper_bound)
      return lb;

    // assert(m == pstarts[n] && g->m == g->pstarts[g->n]);
    for (ui i = 0; i < n; i++)
      degree_q[i] = G1->GetIncidentEdges(i).size();
    for (ui i = 0; i < G2->GetNumVertices(); i++)
      degree_g[i] = G2->GetIncidentEdges(i).size();
    int *degrees_cnt_q = tmp;
    int max_degree_q = 0, max_degree_g = 0;
    memset(degrees_cnt_q, 0, sizeof(int) * n);
    for (ui i = 0; i < n; i++) {
      int td = degree_q[i];
      ++degrees_cnt_q[td];
      if (td > max_degree_q)
        max_degree_q = td;
    }
    int *degrees_cnt_g = degree_q;
    memset(degrees_cnt_g, 0, sizeof(int) * g_n);
    for (ui i = 0; i < g_n; i++) {
      ui td = degree_g[i];
      ++degrees_cnt_g[td];
      if (td > max_degree_g)
        max_degree_g = td;
    }
    ui de = 0, ie = 0;
    while (max_degree_q > 0 && max_degree_g > 0) {
      if (degrees_cnt_q[max_degree_q] == 0) {
        --max_degree_q;
        continue;
      }
      if (degrees_cnt_g[max_degree_g] == 0) {
        --max_degree_g;
        continue;
      }
      ui td =
          std::min(degrees_cnt_q[max_degree_q], degrees_cnt_g[max_degree_g]);
      if (max_degree_q > max_degree_g)
        de += td * (max_degree_q - max_degree_g);
      else
        ie += td * (max_degree_g - max_degree_q);
      degrees_cnt_q[max_degree_q] -= td;
      degrees_cnt_g[max_degree_g] -= td;
    }
    while (max_degree_q > 0) {
      de += max_degree_q * degrees_cnt_q[max_degree_q];
      --max_degree_q;
    }
    while (max_degree_g > 0) {
      ie += max_degree_g * degrees_cnt_g[max_degree_g];
      --max_degree_g;
    }
    de = (de + 1) / 2;
    ie = (ie + 1) / 2;
    ui edge_lb = de + ie;
    if (de * 2 + g_m / 2 > m / 2 && de * 2 + g_m / 2 - m / 2 > edge_lb)
      edge_lb = de * 2 + g_m / 2 - m / 2;
    if (ie * 2 + m / 2 > g_m / 2 && ie * 2 + m / 2 - g_m / 2 > edge_lb)
      edge_lb = ie * 2 + m / 2 - g_m / 2;
    if (lb + edge_lb > verify_upper_bound)
      return lb + edge_lb;
    ui common_elabel_cnt = 0;
    for (ui i = 0; i < m / 2; i++) {
      elabel_cnt[G1->GetEdgeLabel(i)] += 2;
    }
    for (ui i = 0; i < g_m / 2; i++) {
      ui el = G2->GetEdgeLabel(i);
      if (elabel_cnt[el] > 0) {
        elabel_cnt[el] -= 2;
        common_elabel_cnt += 2;
      }
    }
    for (ui i = 0; i < G1->GetNumEdges(); i++) {
      elabel_cnt[G1->GetEdgeLabel(i)] = 0;
    }
    common_elabel_cnt /= 2;
    if (de + g_m / 2 - common_elabel_cnt > edge_lb) {
      edge_lb = de + g_m / 2 - common_elabel_cnt;
    }
    if (ie + m / 2 - common_elabel_cnt > edge_lb) {
      edge_lb = de + m / 2 - common_elabel_cnt;
    }
    ui e_cnt = m;
    if (g_m > e_cnt) {
      e_cnt = g_m;
    }
    e_cnt /= 2;
    if (e_cnt - common_elabel_cnt > edge_lb) {
      edge_lb = e_cnt - common_elabel_cnt;
    }

    // std::cout << lb << " " << edge_lb << "\n";
    return lb + edge_lb;
  }

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
    // if (NaiveCountBound() > threshold)
    //     return false;
    // if (LabelSetDifferenceBound() > threshold)
    //     return false;
    // if (DegreeSequenceBound() > threshold){
    //     return false;
    // }
    // std::cout << BranchBound();
    // if (BranchBound() > threshold)
    // {
    //     return false;
    // }
    // std::cout << __LINE__ << "\n";
    if (ged_lower_bound_filter() > threshold) {
      return false;
    }
    return true;
  }

  /*
   * Functions for GED Computation / Verification
   */

  int GetCurrentBestGED() const { return current_best; }
  ResultLogger GetLog() { return log; }


void heap_top_down(ui idx, ui heap_n, std::vector<std::pair<double, int> > &heap, std::vector<int> &pos) {
	std::pair<double, int> tmp = heap[idx];
	while(2*idx+1 < heap_n) {
		ui i = 2*idx+1;
		if(i+1 < heap_n&&heap[i+1].first > heap[i].first) ++ i;
		if(heap[i].first > tmp.first) {
			heap[idx] = heap[i];
			pos[heap[idx].second] = idx;
			idx = i;
		}
		else break;
	}
	heap[idx] = tmp;
	pos[tmp.second] = idx;
}
void heap_bottom_up(ui idx, std::vector<std::pair<double,int> > &heap, std::vector<int> & pos) {
	std::pair<double, int> tmp = heap[idx];
	while(idx > 0) {
		ui i = (idx-1)/2;
		if(heap[i].first < tmp.first) {
			heap[idx] = heap[i];
			pos[heap[idx].second] = idx;
			idx = i;
		}
		else break;
	}
	heap[idx] = tmp;
	pos[tmp.second] = idx;
}


  void ComputeMatchingOrder() {
    int N = G2->GetNumVertices();
    auto vlabel_freq = G2->GetVertexLabelFrequency();
    auto elabel_freq = G2->GetEdgeLabelFrequency();
    int root = 0;
    int cnt = 0;
    double root_weight = 0;
    for(int i = 0 ; i < G1->GetNumVertices() ; i++){
      double weight = 1 - vlabel_freq[G1->GetVertexLabel(i)]/double(G2->GetNumVertices());
      for(auto j : G1->GetNeighbors(i)){
        weight += 1 - elabel_freq[G1->GetEdgeLabel(i, j)]/double(G2->GetNumEdges() * 2);
      }
      if(weight > root_weight){
        root = i;
        root_weight = weight;
      }
    }
    int q_n = G1->GetNumVertices();
    std::vector<std::pair<double, int> > heap(q_n);
	for(ui i = 0;i < q_n;i ++) {
		if(i == root) heap[i].first = root_weight;
		else heap[i].first = 0;
		heap[i].second = i;
	}
  std::swap(heap[0], heap[root]);
  std::vector<int> pos(q_n);
	for(ui i = 0;i < q_n;i ++) pos[heap[i].second] = i;

  ui heap_n = q_n;

  for(ui i = 0;i < q_n;i ++) {
		ui u = heap[0].second;
		matching_order.push_back(u);
		pos[u] = heap_n-1;
		heap[0] = heap[-- heap_n];
		pos[heap[0].second] = 0;
		heap_top_down(0, heap_n, heap, pos);
		
  
		for(ui v : G1->GetNeighbors(u)){
      if(pos[v] < heap_n){
        ui idx = pos[v];
			  if(heap[idx].first < 10e-6) {
				  heap[idx].first += 1 - vlabel_freq[G1->GetVertexLabel(v)]/double(G2->GetNumVertices());
			  }
			  heap[idx].first += 1 - elabel_freq[G1->GetEdgeLabel(u, v)]/double(G2->GetNumEdges() * 2);
			  // cout << heap[idx].first << " " << heap[idx].second << "\n";
			  heap_bottom_up(idx, heap, pos);
		  }
    }
	}
    // std::vector<double> w(G1->GetNumVertices(), 0);
    // std::vector<int> T(G1->GetNumVertices(), 0);
    // w[root] = root_weight;

    // std::priority_queue<std::pair<double, int>> weighted_queue;
    // weighted_queue.push({w[root], root});

    // while(!weighted_queue.empty()){
    //   int u = weighted_queue.top().second;
    //   weighted_queue.pop();
    //   if(T[u] == 1){
    //     // std::cout << matching_order.size() << " " <<G1->GetNumVertices() <<"\n";
    //     continue;
    //   }

    //   matching_order.push_back(u);
    //   T[u] = 1;
    //   // std::cout <<"pop : " << u << " " <<w[u] << " " << matching_order.size() <<  "\n";
    //   for(auto v : G1->GetNeighbors(u)){
    //     // std::cout << T[v] << "\n";
    //     if(T[v] == 1){
    //       continue;
    //     }
    //     if(w[v] < 10e-6){
    //       w[v] += 1- vlabel_freq[G1->GetVertexLabel(v)]/double(G2->GetNumVertices()); 
    //     }
    //     w[v] += 1 - elabel_freq[G1->GetEdgeLabel(u, v)]/double(G2->GetNumEdges() * 2);
    //     weighted_queue.push({w[v], v});
    //     // std::cout << "push  : " << v << " " << w[v] << "\n";
    //   }
    // }
    // for(int i = 0 ; i < G1->GetNumVertices();i++){
    //   if(T[i] == 0){
    //     matching_order.push_back(i);
    //   }
    // }
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
      // if (inv_u == -1 || inv_v == -1) {
      if (inv_u == -1 || inv_v == -1 || inv_u >= G1->GetNumVertices() ||
          inv_v >= G1->GetNumVertices()) {
        if (verbose)
          printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n", u, v,
                 inv_u, inv_v);
        cost++;
      } else {
        // printf("%d %d\n", inv_u, inv_v);
        // printf("%d\n", G1->GetNumVertices());
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