#pragma once
#include "Base/Hungarian.h"
#include "GraphSimilarity/EditDistance.h"

namespace GraphLib::GraphSimilarity {
struct BMaState : State {
  int next_mapping_order = -1;
  int ub = 1e9;
};

class AStarBMa : public GraphEditDistanceSolver {
  const bool DEBUG = false;
  std::priority_queue<BMaState*, std::vector<BMaState*>, StateComparator> queue;
  std::map<int, int> num_hungarian;

 public:
  void ExtendState(BMaState* state) {
    if (state->cost >= current_best) return;
    int depth = state->depth + 1;
    int u = matching_order[depth];

    for (int v = 0; v < G2->GetNumVertices(); v++) {
      if (state->inverse_mapping[v] != -1) continue;
      BMaState* child_state = new BMaState(state);
      child_state->cost = GetChildEditCost(state, u, v);
      if (DEBUG) {
        fprintf(stderr, "Parent %d has cost %d, child %d has cost %d\n",
                state->id, state->cost, child_state->id, child_state->cost);
      }
      child_state->mapping[u] = v;
      child_state->inverse_mapping[v] = u;
      auto [lb, ub] = BMaLowerBound(child_state);
      child_state->lower_bound = lb;
      child_state->ub = ub;
      if (lb == ub) continue;
      if (depth == G1->GetNumVertices() - 1) {
        if (DEBUG)
          fprintf(stderr, "Reached at leaf of cost %d\n",
                  child_state->lower_bound);
        current_best = std::min(current_best, child_state->lower_bound);
        memcpy(&current_best_mapping[0], child_state->mapping,
               sizeof(int) * NumG1Vertices);
        continue;
      }
      if (child_state->lower_bound >= current_best) continue;
      if (threshold > 0) {
        if (child_state->lower_bound > threshold) continue;
      }
      if (DEBUG) {
        fprintf(stderr, "Parent %d pushes state %d with bound %d\n", state->id,
                child_state->id, child_state->lower_bound);
        child_state->Print();
      }
      queue.push(child_state);
    }
  }

  int GED() {
    PrepareGED();
    BMaState* initial_state = new BMaState(NULL);
    initial_state->cost = 0;
    initial_state->depth = -1;
    auto [lb, ub] = BMaLowerBound(initial_state);
    initial_state->lower_bound = lb;
    queue.push(initial_state);
    int64_t max_qsize = 1;
    while (!queue.empty()) {
      BMaState* current_state = queue.top();
      if (DEBUG) current_state->Print();
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
        queue = std::priority_queue<BMaState*, std::vector<BMaState*>,
                                    StateComparator>();
        break;
      }
      if (threshold >= 0) {
        if (current_best < threshold) {
          queue = std::priority_queue<BMaState*, std::vector<BMaState*>,
                                      StateComparator>();
          break;
        }
        if (current_state->lower_bound > threshold) {
          current_best = -1;
          queue = std::priority_queue<BMaState*, std::vector<BMaState*>,
                                      StateComparator>();
          break;
        }
      }
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

  const int TARGET_ID = -2;

  std::pair<int, int> BMaLowerBound(BMaState* state) {
    DifferenceVector diff;
    diff.init(20);
    int remaining = G2->GetNumVertices() - (state->depth + 1);
    std::vector<int> rem_left, rem_right;
    std::vector<std::vector<int>> branch_distance_matrix(
        remaining, std::vector<int>(remaining, 0));
    ComputeBranchDistanceMatrix(state, branch_distance_matrix, rem_left,
                                rem_right);
    Hungarian hungarian(branch_distance_matrix);
    hungarian.Solve();
    if (DEBUG) hungarian.Print();
    auto& assignment = hungarian.GetAssignment();
    std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
    std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
    std::memcpy(hungarian_mapping.data(), state->mapping,
                sizeof(int) * G1->GetNumVertices());
    std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping,
                sizeof(int) * G2->GetNumVertices());
    for (int i = 0; i < rem_left.size(); i++) {
      int u = rem_left[i];
      int v = rem_right[assignment[i]];
      hungarian_mapping[u] = v;
      hungarian_inverse_mapping[v] = u;
    }
    int ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
    int lb = state->cost + ((hungarian.GetTotalCost() + 1) / 2);
    // int u_idx = 0, best_weight = 9999;
    // int disconnected_best = -1, disconnected_best_weight = 9999;
    // for (int u = 0; u < G1->GetNumVertices(); u++) {
    //   bool is_connected = false;
    //   for (int v : G1->GetNeighbors(u)) {
    //     if (state->mapping[v] != -1) {
    //       is_connected = true;
    //       break;
    //     }
    //   }
    //   if (state->mapping[u] == -1) {
    //     if (is_connected) {
    //       if (hungarian.AssignedWeight(u_idx) < best_weight) {
    //         state->next_mapping_order = u;
    //         best_weight = hungarian.AssignedWeight(u_idx);
    //         if (best_weight < disconnected_best_weight) {
    //           disconnected_best = u;
    //           disconnected_best_weight = best_weight;
    //         }
    //       }
    //     } else {
    //       if (hungarian.AssignedWeight(u_idx) < disconnected_best_weight) {
    //         disconnected_best = u;
    //         disconnected_best_weight = hungarian.AssignedWeight(u_idx);
    //       }
    //     }
    //     u_idx++;
    //   }
    // }
    // if (state->next_mapping_order == -1) {
    //   state->next_mapping_order = disconnected_best;
    // }
    return {lb, ub};
  }
};
}  // namespace GraphLib::GraphSimilarity
