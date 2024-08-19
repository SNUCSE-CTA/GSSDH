#pragma once
#include "GraphSimilarity/EditDistance.h"

DEBUG LSadebugger(2);

namespace GraphLib::GraphSimilarity {
class AStarLSa : public GraphEditDistanceSolver {
 public:
  void ExtendState(State* state) {
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
      if (child_cost + lb >= current_best) continue;
      if (threshold > 0) {
        if (child_cost + lb > threshold) continue;
      }
      State* child_state = new State(state);
      //            printf("State [%p] Check Depth %d [Parent %p]! matching
      //            (%d(%d), %d(%d)) gives lb = %d, cost = %d (curbest = %d)\n",
      //            child_state, depth, state, u, u_label, v, v_label,
      //            lb+child_cost, child_cost, current_best); printf("  LB
      //            Breakdown: VLabel %d, InnerEdge %d, CrossEdge %d\n\n",
      //            vertex_lower_bound, inner_edge_lower_bound,
      //            cross_edge_lower_bound); if (depth == 2) exit(1);
      child_state->cost = child_cost;
      child_state->vertex_label_bound = vertex_lower_bound;
      child_state->inner_edge_label_bound = inner_edge_lower_bound;
      child_state->cross_edge_label_bound = cross_edge_lower_bound;
      child_state->ComputeLowerBound();
      child_state->mapping[u] = v;
      child_state->inverse_mapping[v] = u;
      LSadebugger.log("Generate Child State", 1, 2);
      LSadebugger.log(child_state->to_string(), 1, 2);
      if (depth == G1->GetNumVertices() - 1) {
        current_best = std::min(current_best, child_state->lower_bound);
        memcpy(&current_best_mapping[0], child_state->mapping,
               sizeof(int) * NumG1Vertices);
        continue;
      }
      //                fprintf(stdout,"Push state %d with bound %d\n",
      //                child_state->id, child_state->lower_bound);
      queue.push(child_state);
    }
  }

  int GED() {
    PrepareGED();
    State* initial_state = new State(NULL);
    initial_state->cost = 0;
    initial_state->vertex_label_bound = 0;
    initial_state->vertex_label_bound = vlabel_diff->GetDifference();
    initial_state->inner_edge_label_bound = elabel_diff->GetDifference();
    initial_state->ComputeLowerBound();
    initial_state->depth = -1;
    queue.push(initial_state);
    int64_t max_qsize = 1;
    while (!queue.empty()) {
      State* current_state = queue.top();
      num_nodes++;
      if (num_nodes % LOG_EVERY == 0) {
        fprintf(stderr,
                "%ld, Current best = %d, current cost = %d / lb = %d, depth "
                "%d, Queuesize %lu\n",
                num_nodes, current_best, current_state->cost,
                current_state->lower_bound, current_state->depth, queue.size());
      }
      queue.pop();
      if (current_state->lower_bound >= current_best) {
        queue =
            std::priority_queue<State*, std::vector<State*>, StateComparator>();
        break;
      }
      if (threshold >= 0) {
        if (current_best < threshold) {
          queue = std::priority_queue<State*, std::vector<State*>,
                                      StateComparator>();
          break;
        }
        if (current_state->lower_bound > threshold) {
          current_best = -1;
          queue = std::priority_queue<State*, std::vector<State*>,
                                      StateComparator>();
          break;
        }
      }
      LSadebugger.log("Current QueueTop State", 1);
      LSadebugger.log(current_state->to_string(), 1);
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
}  // namespace GraphLib::GraphSimilarity