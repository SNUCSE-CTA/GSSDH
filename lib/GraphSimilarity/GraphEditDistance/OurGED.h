#pragma once
#include "Base/Hungarian.h"
#include "GraphSimilarity/EditDistance.h"

DEBUG debugger(2);
namespace GraphLib::GraphSimilarity {
int BranchSimilarity(const Branch& a, const Branch& b) {
  int ed = 0;
  if (a.vertex_label == b.vertex_label) ed += 2;
  ed += VectorIntersectionSize(a.edge_labels, b.edge_labels);
  return ed;
}

struct OurState : public State {
  int similarity = 0;

  OurState(OurState* p = NULL, int c = -1, int d = -1) : State(p, c, d) {
    similarity = 0;
    id = global_state_id++;
  };

  OurState(const OurState& other) : State(other) {
    id = global_state_id++;
    parent = other.parent;
    next_mapping_order = other.next_mapping_order;
    cost = other.cost;
    depth = other.depth;
    lower_bound = other.lower_bound;
    vertex_label_bound = other.vertex_label_bound;
    inner_edge_label_bound = other.inner_edge_label_bound;
    cross_edge_label_bound = other.cross_edge_label_bound;
    mapping = new int[NumG1Vertices];
    inverse_mapping = new int[NumG2Vertices];
    std::memcpy(mapping, other.mapping, sizeof(int) * NumG1Vertices);
    std::memcpy(inverse_mapping, other.inverse_mapping,
                sizeof(int) * NumG2Vertices);
  }

  std::string to_string() const {
    std::ostringstream oss;
    oss << "State " << this->id << ": ";
    for (int i = 0; i < NumG1Vertices; i++) {
      oss << this->mapping[i] << " ";
    }
    oss << "| Depth = " << depth << " | Cost = " << cost
        << " | Lb = " << lower_bound << " | Similarity = " << similarity;
    return oss.str();
  }
};

struct OurStateComparator {
  bool operator()(const OurState* a, const OurState* b) const {
    if (a->similarity == b->similarity) {
      if (a->lower_bound == b->lower_bound) {
        if (a->depth == b->depth) {
          return a->id < b->id;
        }
        return a->depth < b->depth;
      }
      return a->lower_bound > b->lower_bound;
    }
    return a->similarity < b->similarity;
  }
};

class OurGED : public GraphEditDistanceSolver {
  std::priority_queue<OurState*, std::vector<OurState*>, OurStateComparator>
      queue;
  std::vector<std::vector<int>> branch_similarities;

 public:
  void ExtendState(OurState* state) {
    if (state->cost >= current_best) return;
    int depth = state->depth + 1;
    int u = matching_order[depth];

    DifferenceVector unmapped_vertex_labels(vlabel_diff->size());
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      if (state->mapping[i] == -1)
        unmapped_vertex_labels.update(G1->GetVertexLabel(i), +1);
    }
    for (int i = 0; i < G2->GetNumVertices(); i++) {
      if (state->inverse_mapping[i] == -1)
        unmapped_vertex_labels.update(G2->GetVertexLabel(i), -1);
    }

    DifferenceVector unmapped_inner_edge_labels(elabel_diff->size());
    std::vector<DifferenceVector> unmapped_cross_edge_labels(
        G1->GetNumVertices(), DifferenceVector(elabel_diff->size()));
    for (int i = 0; i < G1->GetNumEdges(); i++) {
      auto [a, b] = G1->GetEdge(i);
      if (state->mapping[a] == -1) {
        if (state->mapping[b] == -1)
          unmapped_inner_edge_labels.update(G1->GetEdgeLabel(i), +1);
        else
          unmapped_cross_edge_labels[b].update(G1->GetEdgeLabel(i), +1);
      } else {
        if (state->mapping[b] == -1)
          unmapped_cross_edge_labels[a].update(G1->GetEdgeLabel(i), +1);
      }
    }
    for (int i = 0; i < G2->GetNumEdges(); i++) {
      auto [a, b] = G2->GetEdge(i);
      int fa = state->inverse_mapping[a], fb = state->inverse_mapping[b];
      if (state->inverse_mapping[a] == -1) {
        if (state->inverse_mapping[b] == -1)
          unmapped_inner_edge_labels.update(G2->GetEdgeLabel(i), -1);
        else
          unmapped_cross_edge_labels[fb].update(G2->GetEdgeLabel(i), -1);
      } else {
        if (state->inverse_mapping[b] == -1)
          unmapped_cross_edge_labels[fa].update(G2->GetEdgeLabel(i), -1);
      }
    }

    //            fprintf(stderr, "Start extending state %d by mapping %d:\n",
    //            state->id, u);
    for (int v = 0; v < G2->GetNumVertices(); v++) {
      if (state->inverse_mapping[v] != -1) continue;

      int child_cost = GetChildEditCost(state, u, v);
      int vertex_lower_bound = state->vertex_label_bound;
      int inner_edge_lower_bound = state->inner_edge_label_bound;
      int cross_edge_lower_bound = state->cross_edge_label_bound;

      int u_label = G1->GetVertexLabel(u);
      int v_label = G2->GetVertexLabel(v);

      // Update Vertex-label-based Lower Bound
      unmapped_vertex_labels.update(u_label, -1);
      unmapped_vertex_labels.update(v_label, +1);
      // Update Inner-edge-label-based Lower Bound
      for (int uprime : G1->GetNeighbors(u)) {
        int el = G1->GetEdgeLabel(u, uprime);
        //                if (verbosity) printf("Consider %d (label =
        //                %d)\n",uprime, el);
        if (state->mapping[uprime] == -1) {
          //                    inner_edge_lower_bound +=
          //                    flip(unmapped_inner_edge_labels, el, -1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[u], el, +1);
          unmapped_inner_edge_labels.update(el, -1);
          unmapped_cross_edge_labels[u].update(el, +1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[uprime].update(el, -1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[uprime], el,
          //                    -1);
        }
        //                if (verbosity) printf("Inner is %d, Cross is
        //                %d\n",inner_edge_lower_bound, cross_edge_lower_bound);
      }
      for (int vprime : G2->GetNeighbors(v)) {
        int el = G2->GetEdgeLabel(v, vprime);
        int fvprime = state->inverse_mapping[vprime];
        //                if (verbosity) printf("Consider (%d-%d) G2 (label =
        //                %d)\n",v,vprime,el);
        if (fvprime == -1) {
          unmapped_inner_edge_labels.update(el, +1);
          unmapped_cross_edge_labels[u].update(el, -1);
          //                    inner_edge_lower_bound +=
          //                    flip(unmapped_inner_edge_labels, el, +1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[u], el, -1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[fvprime].update(el, +1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[fvprime], el,
          //                    +1);
        }
        //                if (verbosity) printf("Inner is %d, Cross is
        //                %d\n",inner_edge_lower_bound, cross_edge_lower_bound);
      }
      vertex_lower_bound = unmapped_vertex_labels.GetDifference();
      inner_edge_lower_bound = unmapped_inner_edge_labels.GetDifference();
      cross_edge_lower_bound = 0;
      for (auto& it : unmapped_cross_edge_labels) {
        cross_edge_lower_bound += it.GetDifference();
      }
      int lb =
          vertex_lower_bound + inner_edge_lower_bound + cross_edge_lower_bound;

      // Revert everything
      unmapped_vertex_labels.update(u_label, +1);
      unmapped_vertex_labels.update(v_label, -1);
      for (int uprime : G1->GetNeighbors(u)) {
        int el = G1->GetEdgeLabel(u, uprime);
        if (state->mapping[uprime] == -1) {
          unmapped_inner_edge_labels.update(el, +1);
          unmapped_cross_edge_labels[u].update(el, -1);
          //                    flip(unmapped_inner_edge_labels, el, +1);
          //                    flip(unmapped_cross_edge_labels[u], el, -1);
        } else {
          unmapped_cross_edge_labels[uprime].update(el, +1);
          //                    flip(unmapped_cross_edge_labels[uprime], el,
          //                    +1);
        }
      }
      for (int vprime : G2->GetNeighbors(v)) {
        int el = G2->GetEdgeLabel(v, vprime);
        int fvprime = state->inverse_mapping[vprime];
        if (fvprime == -1) {
          unmapped_inner_edge_labels.update(el, -1);
          unmapped_cross_edge_labels[u].update(el, +1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[fvprime].update(el, -1);
        }
      }
      // If decided to proceed with this...
      //                fprintf(stderr, "  Consider %d-%d: cost %d + future lb
      //                %d = %d\n", u, v, child_cost, lb, child_cost+lb);
      if (child_cost + lb >= current_best) continue;
      if (threshold >= 0) {
        if (child_cost + lb > threshold) continue;
      }
      OurState* child_state = new OurState(*state);
      child_state->cost = child_cost;
      child_state->depth = state->depth + 1;
      child_state->vertex_label_bound = vertex_lower_bound;
      child_state->inner_edge_label_bound = inner_edge_lower_bound;
      child_state->cross_edge_label_bound = cross_edge_lower_bound;
      child_state->similarity = state->similarity + branch_similarities[u][v];
      child_state->ComputeLowerBound();
      child_state->mapping[u] = v;
      child_state->inverse_mapping[v] = u;
      debugger.log("Generate Child State", 1, 2);
      debugger.log(child_state->to_string(), 1, 2);
      if (depth == G1->GetNumVertices() - 1) {
        current_best = std::min(current_best, child_state->lower_bound);
        memcpy(&current_best_mapping[0], child_state->mapping,
               sizeof(int) * NumG1Vertices);
        continue;
      }
      //                fprintf(stderr,"  Push state %d with bound %d\n",
      //                child_state->id, child_state->lower_bound);
      queue.push(child_state);
    }
  }

  int GED() {
    PrepareGED();
    branch_similarities.clear();
    branch_similarities.resize(G1->GetNumVertices(),
                               std::vector<int>(G2->GetNumVertices(), 0));
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      for (int j = 0; j < G2->GetNumVertices(); j++) {
        branch_similarities[i][j] =
            BranchSimilarity(G1->GetBranch(i), G2->GetBranch(j));
      }
    }
    OurState* initial_state = new OurState(NULL);
    initial_state->cost = 0;
    initial_state->similarity = 0;
    initial_state->vertex_label_bound = 0;
    initial_state->vertex_label_bound = vlabel_diff->GetDifference();
    initial_state->inner_edge_label_bound = elabel_diff->GetDifference();
    initial_state->ComputeLowerBound();
    initial_state->depth = -1;
    queue.push(initial_state);
    int64_t max_qsize = 1;
    debugger.log("I am Here", 1);
    while (!queue.empty()) {
      OurState* current_state = queue.top();
      num_nodes++;
      queue.pop();
      if (current_state->lower_bound >= current_best) {
        continue;
      }
      debugger.log("Current QueueTop State", 1);
      debugger.log(current_state->to_string(), 1);
      ExtendState(current_state);
      max_qsize = std::max(max_qsize, (int64_t)queue.size());
    }
    if (threshold >= 0 and current_best > threshold) {
      current_best = -1;
    }
    log.AddResult("MaxQueueSize", max_qsize, RESULT_INT64);
    log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
    log.AddResult("EditDistance", current_best, RESULT_INT);
    return current_best;
  }
};

class OurGEDwithBM : public GraphEditDistanceSolver {
  std::priority_queue<OurState*, std::vector<OurState*>, OurStateComparator>
      queue;
  std::vector<std::vector<int>> branch_similarities;
  std::vector<int> inv_matching_order;
  long long num_pushed_nodes = 0;
  const bool DEBUG = false;
  std::map<int, int> num_hungarian;

 public:
  void ExtendState(OurState* state, int u = -1, int except = -1) {
    if (state->cost >= current_best) return;
    int depth = state->depth + 1;
    if (u == -1) u = matching_order[depth];

    DifferenceVector unmapped_vertex_labels(vlabel_diff->size());
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      if (state->mapping[i] == -1)
        unmapped_vertex_labels.update(G1->GetVertexLabel(i), +1);
    }
    for (int i = 0; i < G2->GetNumVertices(); i++) {
      if (state->inverse_mapping[i] == -1)
        unmapped_vertex_labels.update(G2->GetVertexLabel(i), -1);
    }

    DifferenceVector unmapped_inner_edge_labels(elabel_diff->size());
    std::vector<DifferenceVector> unmapped_cross_edge_labels(
        G1->GetNumVertices(), DifferenceVector(elabel_diff->size()));
    for (int i = 0; i < G1->GetNumEdges(); i++) {
      auto [a, b] = G1->GetEdge(i);
      if (state->mapping[a] == -1) {
        if (state->mapping[b] == -1)
          unmapped_inner_edge_labels.update(G1->GetEdgeLabel(i), +1);
        else
          unmapped_cross_edge_labels[b].update(G1->GetEdgeLabel(i), +1);
      } else {
        if (state->mapping[b] == -1)
          unmapped_cross_edge_labels[a].update(G1->GetEdgeLabel(i), +1);
      }
    }
    for (int i = 0; i < G2->GetNumEdges(); i++) {
      auto [a, b] = G2->GetEdge(i);
      int fa = state->inverse_mapping[a], fb = state->inverse_mapping[b];
      if (state->inverse_mapping[a] == -1) {
        if (state->inverse_mapping[b] == -1)
          unmapped_inner_edge_labels.update(G2->GetEdgeLabel(i), -1);
        else
          unmapped_cross_edge_labels[fb].update(G2->GetEdgeLabel(i), -1);
      } else {
        if (state->inverse_mapping[b] == -1)
          unmapped_cross_edge_labels[fa].update(G2->GetEdgeLabel(i), -1);
      }
    }

    //            fprintf(stderr, "Start extending state %d by mapping %d:\n",
    //            state->id, u);
    for (int v = 0; v < G2->GetNumVertices(); v++) {
      if (state->inverse_mapping[v] != -1) continue;

      int child_cost = GetChildEditCost(state, u, v);
      int vertex_lower_bound = state->vertex_label_bound;
      int inner_edge_lower_bound = state->inner_edge_label_bound;
      int cross_edge_lower_bound = state->cross_edge_label_bound;

      int u_label = G1->GetVertexLabel(u);
      int v_label = G2->GetVertexLabel(v);

      // Update Vertex-label-based Lower Bound
      unmapped_vertex_labels.update(u_label, -1);
      unmapped_vertex_labels.update(v_label, +1);
      // Update Inner-edge-label-based Lower Bound
      for (int uprime : G1->GetNeighbors(u)) {
        int el = G1->GetEdgeLabel(u, uprime);
        //                if (verbosity) printf("Consider %d (label =
        //                %d)\n",uprime, el);
        if (state->mapping[uprime] == -1) {
          //                    inner_edge_lower_bound +=
          //                    flip(unmapped_inner_edge_labels, el, -1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[u], el, +1);
          unmapped_inner_edge_labels.update(el, -1);
          unmapped_cross_edge_labels[u].update(el, +1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[uprime].update(el, -1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[uprime], el,
          //                    -1);
        }
        //                if (verbosity) printf("Inner is %d, Cross is
        //                %d\n",inner_edge_lower_bound, cross_edge_lower_bound);
      }
      for (int vprime : G2->GetNeighbors(v)) {
        int el = G2->GetEdgeLabel(v, vprime);
        int fvprime = state->inverse_mapping[vprime];
        //                if (verbosity) printf("Consider (%d-%d) G2 (label =
        //                %d)\n",v,vprime,el);
        if (fvprime == -1) {
          unmapped_inner_edge_labels.update(el, +1);
          unmapped_cross_edge_labels[u].update(el, -1);
          //                    inner_edge_lower_bound +=
          //                    flip(unmapped_inner_edge_labels, el, +1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[u], el, -1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[fvprime].update(el, +1);
          //                    cross_edge_lower_bound +=
          //                    flip(unmapped_cross_edge_labels[fvprime], el,
          //                    +1);
        }
        //                if (verbosity) printf("Inner is %d, Cross is
        //                %d\n",inner_edge_lower_bound, cross_edge_lower_bound);
      }
      vertex_lower_bound = unmapped_vertex_labels.GetDifference();
      inner_edge_lower_bound = unmapped_inner_edge_labels.GetDifference();
      cross_edge_lower_bound = 0;
      for (auto& it : unmapped_cross_edge_labels) {
        cross_edge_lower_bound += it.GetDifference();
      }
      int lb =
          vertex_lower_bound + inner_edge_lower_bound + cross_edge_lower_bound;

      // Revert everything
      unmapped_vertex_labels.update(u_label, +1);
      unmapped_vertex_labels.update(v_label, -1);
      for (int uprime : G1->GetNeighbors(u)) {
        int el = G1->GetEdgeLabel(u, uprime);
        if (state->mapping[uprime] == -1) {
          unmapped_inner_edge_labels.update(el, +1);
          unmapped_cross_edge_labels[u].update(el, -1);
          //                    flip(unmapped_inner_edge_labels, el, +1);
          //                    flip(unmapped_cross_edge_labels[u], el, -1);
        } else {
          unmapped_cross_edge_labels[uprime].update(el, +1);
          //                    flip(unmapped_cross_edge_labels[uprime], el,
          //                    +1);
        }
      }
      for (int vprime : G2->GetNeighbors(v)) {
        int el = G2->GetEdgeLabel(v, vprime);
        int fvprime = state->inverse_mapping[vprime];
        if (fvprime == -1) {
          unmapped_inner_edge_labels.update(el, -1);
          unmapped_cross_edge_labels[u].update(el, +1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[fvprime].update(el, -1);
        }
      }
      // If decided to proceed with this...
      //                fprintf(stderr, "  Consider %d-%d: cost %d + future lb
      //                %d = %d\n", u, v, child_cost, lb, child_cost+lb);
      if (child_cost + lb >= current_best) continue;
      if (threshold >= 0) {
        if (child_cost + lb > threshold) continue;
      }
      OurState* child_state = new OurState(*state);
      child_state->cost = child_cost;
      child_state->depth = state->depth + 1;
      child_state->vertex_label_bound = vertex_lower_bound;
      child_state->inner_edge_label_bound = inner_edge_lower_bound;
      child_state->cross_edge_label_bound = cross_edge_lower_bound;
      child_state->similarity = state->similarity + branch_similarities[u][v];
      child_state->ComputeLowerBound();
      child_state->mapping[u] = v;
      child_state->inverse_mapping[v] = u;
      debugger.log("Generate Child State", 1, 2);
      debugger.log(child_state->to_string(), 1, 2);
      if (depth == G1->GetNumVertices() - 1) {
        current_best = std::min(current_best, child_state->lower_bound);
        memcpy(&current_best_mapping[0], child_state->mapping,
               sizeof(int) * NumG1Vertices);
        continue;
      }
      num_pushed_nodes++;
      //                fprintf(stderr,"  Push state %d with bound %d\n",
      //                child_state->id, child_state->lower_bound);
      queue.push(child_state);
    }
  }

  void ExtendHungarianToLeaf(OurState* state) {
    if (DEBUG) {
      fprintf(stderr, "ExtendHungarianToLeaf(%d)\n", state->id);
      state->Print();
    }
    DifferenceVector diff;
    diff.init(20);
    int remaining = G2->GetNumVertices() - (state->depth + 1);
    std::vector<int> rem_left, rem_right;
    std::vector<std::vector<int>> branch_distance_matrix(
        remaining, std::vector<int>(remaining, 0));
    ComputeBranchDistanceMatrix(state, branch_distance_matrix, rem_left,
                                rem_right);
    num_hungarian[rem_left.size()]++;
    Hungarian hungarian(branch_distance_matrix);
    hungarian.Solve();
    auto& assignment = hungarian.GetAssignment();

    std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
    std::vector<int> hungarian_weight(G1->GetNumVertices(), 0);
    for (int i = 0; i < rem_left.size(); i++) {
      int u = rem_left[i];
      int v = rem_right[assignment[i]];
      hungarian_mapping[u] = v;
      hungarian_weight[u] = hungarian.AssignedWeight(i);
    }
    std::sort(rem_left.begin(), rem_left.end(), [&](int i, int j) {
      return inv_matching_order[i] < inv_matching_order[j];
      // return hungarian_weight[i] < hungarian_weight[j];
    });
    if (DEBUG and G1->GetId() == 21392 and G2->GetId() == 21392) {
      hungarian.Print();
      std::cerr << "rem_left: " << rem_left << "\n";
      std::cerr << "rem_right: " << rem_right << "\n";
      fprintf(stderr, "Hungarian ok! We will assign:\n");
      for (int i = 0; i < rem_left.size(); i++) {
        int u = rem_left[i];
        int v = rem_right[assignment[i]];
        fprintf(stderr, "%d -> %d\n", u, v);
      }
      fprintf(stderr, "\n");
    }
    int hcost = hungarian.GetTotalCost();
    int lb = (hcost + 1) / 2 + state->cost;
    state->lower_bound = lb;
    if (state->lower_bound < current_best)
      ExtendState(state, rem_left[0], hungarian_mapping[rem_left[0]]);
    for (int i = 0; i < rem_left.size(); i++) {
      int u = rem_left[i];
      int v = hungarian_mapping[u];
      state->cost = GetChildEditCost(state, u, v);
      hcost -= hungarian_weight[u];
      state->lower_bound =
          std::max(state->lower_bound, (hcost + 1) / 2 + state->cost);
      //                state->lower_bound = std::max(lb, state->cost);
      state->mapping[u] = v;
      state->inverse_mapping[v] = u;
      state->similarity = state->similarity + branch_similarities[u][v];
      state->depth++;
      //                fprintf(stderr, "State %d of %d-%d became a state with
      //                cost %d and lb %d at depth %d\n", state->id,
      //                        G1->GetId(), G2->GetId(), state->cost,
      //                        state->lower_bound, state->depth);
      if (i + 1 < rem_left.size()) {
        int nxt_u = rem_left[i + 1];
        int nxt_v = hungarian_mapping[nxt_u];
        if (state->lower_bound < current_best)
          ExtendState(state, nxt_u, nxt_v);
        else
          return;
      }
    }
    if (DEBUG) {
      fprintf(stderr, "We reached to the leaf!\n");
      state->Print();
    }
    current_best = std::min(current_best, ComputeDistance(state));
    //            if (G1->GetId() == 6707 && G2->GetId() == 5025) exit(3);
  }

  int GED() {
    PrepareGED();
    branch_similarities.clear();
    branch_similarities.resize(G1->GetNumVertices(),
                               std::vector<int>(G2->GetNumVertices(), 0));
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      for (int j = 0; j < G2->GetNumVertices(); j++) {
        branch_similarities[i][j] =
            BranchSimilarity(G1->GetBranch(i), G2->GetBranch(j));
      }
    }
    inv_matching_order.resize(G1->GetNumVertices(), -1);
    for (int i = 0; i < G1->GetNumVertices(); i++) {
      inv_matching_order[matching_order[i]] = i;
    }
    OurState* initial_state = new OurState(NULL);
    initial_state->cost = 0;
    initial_state->similarity = 0;
    initial_state->vertex_label_bound = 0;
    initial_state->vertex_label_bound = vlabel_diff->GetDifference();
    initial_state->inner_edge_label_bound = elabel_diff->GetDifference();
    initial_state->ComputeLowerBound();
    initial_state->depth = -1;
    ExtendHungarianToLeaf(initial_state);
    int64_t max_qsize = 1;
    debugger.log("I am Here", 1);
    while (!queue.empty()) {
      OurState* current_state = queue.top();
      num_nodes++;
      if (num_nodes % LOG_EVERY == 0) {
        fprintf(stderr, "GED of %d and %d\n", G1->GetId(), G2->GetId());
        fprintf(stderr,
                "%llu, Current best = %d, current cost = %d / lb = %d, depth "
                "%d, Queuesize %lu\n",
                num_nodes, current_best, current_state->cost,
                current_state->lower_bound, current_state->depth, queue.size());
      }
      queue.pop();
      if (current_state->lower_bound >= current_best) {
        continue;
      }
      if (threshold >= 0) {
        if (current_best < threshold) {
          queue = std::priority_queue<OurState*, std::vector<OurState*>,
                                      OurStateComparator>();
          break;
        }
        if (current_state->lower_bound > threshold) {
          continue;
        }
      }
      ExtendHungarianToLeaf(current_state);
      max_qsize = std::max(max_qsize, (int64_t)queue.size());
    }
    if (threshold >= 0 and current_best > threshold) {
      current_best = -1;
    }
    log.AddResult("MaxQueueSize", max_qsize, RESULT_INT64);
    log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
    log.AddResult("EditDistance", current_best, RESULT_INT);
    return current_best;
  }
};

}  // namespace GraphLib::GraphSimilarity