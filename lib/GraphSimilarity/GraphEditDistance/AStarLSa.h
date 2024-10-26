#pragma once
#include "GraphSimilarity/EditDistance.h"

DEBUG LSadebugger(2);

namespace GraphLib::GraphSimilarity {
class AStarLSa : public GraphEditDistanceSolver {
public:
  int *mapping_path, *inverse_mapping_path;
  int ub;
  double hungarian_time;
  void ExtendState(State *state, GWL *gwl, ColorTree **prev_color_to_node,
                   ColorTree **curr_color_to_node) {
    if (state->cost >= current_best)
      return;
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

    if (ub < current_best) {
      current_best = ub;
      if (current_best <= threshold) {
        // std::cout << "Early Termination because of UB = " << ub << std::endl;
        return;
      }
    }

    for (int v = 0; v < G2->GetNumVertices(); v++) {
      if (state->inverse_mapping[v] != -1)
        continue;

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
        // if (verbosity) printf("Consider %d (label = %d)\n ", uprime, el);
        if (state->mapping[uprime] == -1) {
          // inner_edge_lower_bound += flip(unmapped_inner_edge_labels, el, -1);
          // cross_edge_lower_bound += flip(unmapped_cross_edge_labels[u], el,
          // +1);
          unmapped_inner_edge_labels.update(el, -1);
          unmapped_cross_edge_labels[u].update(el, +1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[uprime].update(el, -1);
          // cross_edge_lower_bound +=
          //     flip(unmapped_cross_edge_labels[uprime], el, -1);
        }
        // if (verbosity)
        //   printf("Inner is %d, Cross is %d\n", inner_edge_lower_bound,
        //          cross_edge_lower_bound);
      }
      for (int vprime : G2->GetNeighbors(v)) {
        int el = G2->GetEdgeLabel(v, vprime);
        int fvprime = state->inverse_mapping[vprime];
        // if (verbosity)
        //   printf("Consider (%d-%d) G2 (label = %d)\n", v, vprime, el);
        if (fvprime == -1) {
          unmapped_inner_edge_labels.update(el, +1);
          unmapped_cross_edge_labels[u].update(el, -1);
          // inner_edge_lower_bound += flip(unmapped_inner_edge_labels, el, +1);
          // cross_edge_lower_bound += flip(unmapped_cross_edge_labels[u], el,
          // -1);
        } else {
          // uPrime is matched: anchored on uprime ->removed
          unmapped_cross_edge_labels[fvprime].update(el, +1);
          // cross_edge_lower_bound +=
          //     flip(unmapped_cross_edge_labels[fvprime], el, +1);
        }
        // if (verbosity)
        //   printf("Inner is %d, Cross is %d\n", inner_edge_lower_bound,
        //          cross_edge_lower_bound);
      }
      vertex_lower_bound = unmapped_vertex_labels.GetDifference();
      inner_edge_lower_bound = unmapped_inner_edge_labels.GetDifference();
      cross_edge_lower_bound = 0;
      for (auto &it : unmapped_cross_edge_labels) {
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
          // flip(unmapped_inner_edge_labels, el, +1);
          // flip(unmapped_cross_edge_labels[u], el, -1);
        } else {
          unmapped_cross_edge_labels[uprime].update(el, +1);
          // flip(unmapped_cross_edge_labels[uprime], el, +1);
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
      if (child_cost + lb >= current_best)
        continue;
      if (threshold > 0) {
        if (child_cost + lb > threshold)
          continue;
      }
      State *child_state = new State(state);
      // printf(
      //     "State [%p] Check Depth %d [Parent %p]! matching (% d(% d), % d(% "
      //     "d)) gives lb = % d, cost = % d(curbest = % d)\n ",
      //     child_state, depth, state, u, u_label, v, v_label, lb + child_cost,
      //     child_cost, current_best);
      // printf("  LB Breakdown : VLabel % d, InnerEdge % d, CrossEdge % d\n\n
      // ",
      //        vertex_lower_bound, inner_edge_lower_bound,
      //        cross_edge_lower_bound);
      // if (depth == 2)
      //   exit(1);
      child_state->cost = child_cost;
      child_state->vertex_label_bound = vertex_lower_bound;
      child_state->inner_edge_label_bound = inner_edge_lower_bound;
      child_state->cross_edge_label_bound = cross_edge_lower_bound;
      child_state->ComputeLowerBound();
      child_state->mapping[u] = v;
      child_state->inverse_mapping[v] = u;

      if (mapping_path[u] == v) {
        child_state->aux_lower_bound = child_state->lower_bound;
        child_state->lower_bound = -1;
      }

      LSadebugger.log("Generate Child State", 1, 2);
      LSadebugger.log(child_state->to_string(), 1, 2);
      if (depth == G1->GetNumVertices() - 1) {
        current_best = std::min(current_best, child_state->lower_bound);
        memcpy(&current_best_mapping[0], child_state->mapping,
               sizeof(int) * NumG1Vertices);
        continue;
      }
      // fprintf(stdout, "Push state %d with bound %d\n", child_state->id,
      //         child_state->lower_bound);
      queue.push(child_state);
    }
  }

  int GED(GSSEntry *combined, GWL *gwl, ColorTree **prev_color_to_node,
          ColorTree **curr_color_to_node) {
    int initial_ub = 0;
    hungarian_time = 0;
    PrepareGED(combined);
    std::fill_n(prev_color_to_node, combined->GetNumVertices(), nullptr);
    std::fill_n(curr_color_to_node, combined->GetNumVertices(), nullptr);
    int mapping[combined->GetNumVertices()];
    int inverse_mapping[combined->GetNumVertices()];
    std::fill_n(mapping, combined->GetNumVertices(), -1);
    std::fill_n(inverse_mapping, combined->GetNumVertices(), -1);
    gwl->Init();
    gwl->GraphColoring(5, mapping, inverse_mapping, prev_color_to_node,
                       curr_color_to_node);
    Timer timer;
    timer.Start();
    // gwl->MatchByHungarian(mapping, inverse_mapping);
    gwl->VertexMatching();
    timer.Stop();
    hungarian_time += timer.GetTime();
    mapping_path = new int[combined->GetNumVertices()];
    inverse_mapping_path = new int[combined->GetNumVertices()];
    for (int i = 0; i < combined->combined_index; i++) {
      mapping_path[i] = gwl->mapping[i];
    }
    for (int i = 0; i < combined->GetNumVertices() - combined->combined_index;
         i++) {
      inverse_mapping_path[i] = gwl->inverse_mapping[i];
    }
    ub = ComputeDistance(gwl->mapping, gwl->inverse_mapping, false);
    initial_ub = ub;
    // std::cout << "Coloring gave UB = " << ub << std::endl;
    gwl->Deallocate();

    State *initial_state = new State(NULL);
    initial_state->cost = 0;
    initial_state->vertex_label_bound = 0;
    initial_state->vertex_label_bound = vlabel_diff->GetDifference();
    initial_state->inner_edge_label_bound = elabel_diff->GetDifference();
    initial_state->ComputeLowerBound();
    initial_state->depth = -1;
    queue.push(initial_state);
    int64_t max_qsize = 1;
    State *current_state = nullptr;
    while (!queue.empty()) {
      current_state = queue.top();
      num_nodes++;
      if (num_nodes % LOG_EVERY == 0) {
        fprintf(stderr,
                "%ld, Current best = %d, current cost = %d / lb = %d, depth "
                "%d, Queuesize %lu\n",
                num_nodes, current_best, current_state->cost,
                current_state->lower_bound, current_state->depth, queue.size());
      }
      queue.pop();
      if (current_state->lower_bound == -1) {
        current_state->lower_bound = current_state->aux_lower_bound;
        // Debug
        // std::cout << "Directly set lower bound to "
        //           << current_state->lower_bound << std::endl;
        // std::cout << "Matching: ";
        // for (int i = 0; i < combined->combined_index; i++) {
        //   std::cout << i << " -> " << current_state->mapping[i] << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "Mapping path: ";
        // for (int i = 0; i < combined->combined_index; i++) {
        //   std::cout << i << " -> " << mapping_path[i] << " ";
        // }
        // std::cout << std::endl;
      } else {
        // Debug
        // std::cout << "Coloring and Hungarian" << std::endl;
        std::fill_n(prev_color_to_node, combined->GetNumVertices(), nullptr);
        std::fill_n(curr_color_to_node, combined->GetNumVertices(), nullptr);
        gwl->Init();
        gwl->GraphColoring(5, current_state->mapping,
                           current_state->inverse_mapping, prev_color_to_node,
                           curr_color_to_node);
        Timer timer;
        timer.Start();
        // gwl->MatchByHungarian(current_state->mapping,
        //                       current_state->inverse_mapping);
        gwl->VertexMatching();
        timer.Stop();
        hungarian_time += timer.GetTime();
        for (int i = 0; i < combined->combined_index; i++) {
          mapping_path[i] = gwl->mapping[i];
        }
        for (int i = 0;
             i < combined->GetNumVertices() - combined->combined_index; i++) {
          inverse_mapping_path[i] = gwl->inverse_mapping[i];
        }
        ub = ComputeDistance(gwl->mapping, gwl->inverse_mapping, false);
        // std::cout << "Coloring gave UB = " << ub << std::endl;
        gwl->Deallocate();
      }
      if (current_state->lower_bound >= current_best) {
        Deallocation();
        queue = std::priority_queue<State *, std::vector<State *>,
                                    StateComparator>();
        break;
      }
      if (threshold >= 0) {
        if (current_best < threshold) {
          Deallocation();
          queue = std::priority_queue<State *, std::vector<State *>,
                                      StateComparator>();
          break;
        }
        if (current_state->lower_bound > threshold) {
          Deallocation();
          current_best = -1;
          queue = std::priority_queue<State *, std::vector<State *>,
                                      StateComparator>();
          break;
        }
      }
      LSadebugger.log("Current QueueTop State", 1);
      LSadebugger.log(current_state->to_string(), 1);
      ExtendState(current_state, gwl, prev_color_to_node, curr_color_to_node);
      if (threshold >= 0) {
        if (current_best < threshold) {
          Deallocation();
          queue = std::priority_queue<State *, std::vector<State *>,
                                      StateComparator>();
          break;
        }
      }
      max_qsize = std::max(max_qsize, (int64_t)queue.size());
      delete current_state;
    }
    std::cout << "query id: " << G1->GetId() << ", target id: " << G2->GetId()
              << std::endl;
    std::cout << "Initial ub:" << initial_ub << ", ub: " << ub
              << ", current_best: " << current_best << std::endl;
    if (threshold >= 0 and current_best > threshold) {
      current_best = -1;
    }
    log.AddResult("MaxQueueSize", max_qsize, RESULT_INT64);
    log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
    log.AddResult("EditDistance", current_best, RESULT_INT);
    return current_best;
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
    // current_best = std::min(cost, current_best);
    return cost;
  }

  void Deallocation() {
    while (!queue.empty()) {
      State *state = queue.top();
      queue.pop();
      delete state;
    }
  }
};
} // namespace GraphLib::GraphSimilarity