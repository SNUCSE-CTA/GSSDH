#pragma once
#include "Base/Timer.h"
#include "StableColoring.h"
static int global_state_idx = 0;
namespace GraphLib {
class GraphIsomorphism {
  LabeledGraph *GL, *GR;
  PairwiseStableColoring *coloring;
  bool *color_mapped;
  int VL, VR;
  uint64_t num_backtrack_states = 0;

  struct State {
    int depth, u, v_candidate_index;
    std::vector<int> v_candidates, mapped_vertices;
    int state_idx;
    State(int depth, int v_candidate_idx) {
      this->depth = depth;
      this->v_candidate_index = v_candidate_idx;
      this->state_idx = global_state_idx++;
      this->u = -1;
    }
  };

 public:
  LabeledGraph *GCombined;
  GraphIsomorphism(LabeledGraph *g1, LabeledGraph *g2) {
    GL = g1, GR = g2;
    VL = GL->GetNumVertices(), VR = GR->GetNumVertices();
    color_mapped = new bool[GL->GetNumVertices()];
  }

  std::vector<int> matching_order;
  void GetMatchingOrder() {
    std::vector<std::pair<int, int>> color_order;
    for (int i = 0; i < coloring->GetNumColors(); i++) {
      if (coloring->GetLeftColorCount(i) != 1) {
        color_order.emplace_back(coloring->GetLeftColorCount(i), i);
      }
    }
    std::sort(color_order.begin(), color_order.end());
    std::reverse(color_order.begin(), color_order.end());
    for (auto &[_, c] : color_order) {
      for (auto &v : coloring->GetVerticesByColor(c)) {
        if (v < VL) matching_order.emplace_back(v);
      }
    }
  }

  std::stack<State> stk;
  bool Backtrack() {
    stk.emplace(0, 0);
    while (!stk.empty()) {
      State &current = stk.top();
      num_backtrack_states++;
      if (current.depth >= matching_order.size()) {
        return true;
      }
      if (current.u == -1) {
        current.u = matching_order[current.depth];
        int u_color = coloring->GetColor(current.u);
        for (int x : coloring->GetVerticesByColor(u_color)) {
          if (x >= VL) current.v_candidates.push_back(x);
        }
      }
      while (current.v_candidate_index < current.v_candidates.size()) {
        int v_cand = current.v_candidates[current.v_candidate_index++];
        bool ok = coloring->PairwiseIndividualize(current.u, v_cand);
        if (ok) {
          int next_depth = current.depth + 1;
          while (next_depth < matching_order.size() and
                 coloring->GetNumVertexClass(matching_order[next_depth]) == 1) {
            next_depth++;
          }
          stk.emplace(next_depth, 0);
          goto nxt_state;
        }
        coloring->RevertColoring();
      }
      if (current.v_candidate_index >= current.v_candidates.size()) {
        if (current.state_idx == 0) return false;
        coloring->RevertColoring();
        stk.pop();
      }
    nxt_state:
      continue;
    }
    return false;
  }

  bool isIsomorphic() {
    if (GL->GetNumVertices() != GR->GetNumVertices()) return false;
    if (GL->GetNumEdges() != GR->GetNumEdges()) return false;
    Timer rt;
    rt.Start();
    coloring = new PairwiseStableColoring(GL, GR, GCombined);
    coloring->Refine();
    rt.Stop();
    rt.GetTime();
    fprintf(stdout, "Stable #Colors: %d\n", coloring->GetNumColors());
    fprintf(stdout, "Refinement Time: %.02lf ms\n", rt.GetTime());
    fflush(stdout);
    if (!coloring->CheckColorCount()) return false;
    GetMatchingOrder();
    if (matching_order.empty()) return true;
    bool result = Backtrack();
    return result;
  }
};
}  // namespace GraphLib
