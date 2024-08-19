#pragma once
#include "Base/Hungarian.h"
#include "GraphSimilarity/EditDistance.h"


namespace GraphLib::GraphSimilarity {
    class AStarMixed: public GraphEditDistanceSolver {
        std::map<int, int> num_hungarian;
        const bool DEBUG = false;
        std::priority_queue<State*, std::vector<State*>, StateComparator> queue;
        std::vector<int> inv_matching_order;
        long long num_pushed_nodes = 0;
    public:
        void ExtendStateLSa(State *state, int u = -1, int except = -1) {
            if (u == -1 || except == -1) {
                fprintf(stderr, "ExtendStateLSa %d: u = %d, except = %d\n", state->id, u, except);
                exit(5);
            }
            if (DEBUG) {
                fprintf(stderr, "ExtendStateLSa %d Except child %d\n", state->id, except);
                state->Print();
            }
            if (state->cost >= current_best) return;
            int depth = state->depth + 1;

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
            std::vector<DifferenceVector> unmapped_cross_edge_labels(G1->GetNumVertices(), DifferenceVector(elabel_diff->size()));
            for (int i = 0; i < G1->GetNumEdges(); i++) {
                auto [a, b] = G1->GetEdge(i);
                if (state->mapping[a] == -1) {
                    if (state->mapping[b] == -1)
                        unmapped_inner_edge_labels.update(G1->GetEdgeLabel(i), +1);
                    else
                        unmapped_cross_edge_labels[b].update(G1->GetEdgeLabel(i), +1);
                }
                else {
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
                }
                else {
                    if (state->inverse_mapping[b] == -1)
                        unmapped_cross_edge_labels[fa].update(G2->GetEdgeLabel(i), -1);
                }
            }
            for (int v = 0; v < G2->GetNumVertices(); v++) {
                if (state->inverse_mapping[v] != -1) continue;
                if (v == except) continue;
                int child_cost = GetChildEditCost(state, u, v);
                int vertex_lower_bound = state->vertex_label_bound;
                int inner_edge_lower_bound = state->inner_edge_label_bound;
                int cross_edge_lower_bound = state->cross_edge_label_bound;

                int u_label = G1->GetVertexLabel(u), v_label = G2->GetVertexLabel(v);

                unmapped_vertex_labels.update(u_label, -1);
                unmapped_vertex_labels.update(v_label, +1);
                for (int uprime : G1->GetNeighbors(u)) {
                    int el = G1->GetEdgeLabel(u, uprime);
                    if (state->mapping[uprime] == -1) {
                        unmapped_inner_edge_labels.update(el, -1);
                        unmapped_cross_edge_labels[u].update(el, +1);
                    }
                    else {
                        unmapped_cross_edge_labels[uprime].update(el, -1);
                    }
                }
                for (int vprime : G2->GetNeighbors(v)) {
                    int el = G2->GetEdgeLabel(v, vprime);
                    int fvprime = state->inverse_mapping[vprime];
                    if (fvprime == -1) {
                        unmapped_inner_edge_labels.update(el, +1);
                        unmapped_cross_edge_labels[u].update(el, -1);
                    }
                    else {
                        unmapped_cross_edge_labels[fvprime].update(el, +1);
                    }
                }
                vertex_lower_bound = unmapped_vertex_labels.GetDifference();
                inner_edge_lower_bound = unmapped_inner_edge_labels.GetDifference();
                cross_edge_lower_bound = 0;
                for (auto &it : unmapped_cross_edge_labels) {
                    cross_edge_lower_bound += it.GetDifference();
                }
                int lb = vertex_lower_bound + inner_edge_lower_bound + cross_edge_lower_bound;

                unmapped_vertex_labels.update(u_label, +1);
                unmapped_vertex_labels.update(v_label, -1);
                for (int uprime : G1->GetNeighbors(u)) {
                    int el = G1->GetEdgeLabel(u, uprime);
                    if (state->mapping[uprime] == -1) {
                        unmapped_inner_edge_labels.update(el, +1);
                        unmapped_cross_edge_labels[u].update(el, -1);
                    }
                    else {
                        unmapped_cross_edge_labels[uprime].update(el, +1);
                    }
                }
                for (int vprime : G2->GetNeighbors(v)) {
                    int el = G2->GetEdgeLabel(v, vprime);
                    int fvprime = state->inverse_mapping[vprime];
                    if (fvprime == -1) {
                        unmapped_inner_edge_labels.update(el, -1);
                        unmapped_cross_edge_labels[u].update(el, +1);
                    }
                    else {
                        unmapped_cross_edge_labels[fvprime].update(el, -1);
                    }
                }
                if (child_cost + lb >= current_best) continue;
                if (threshold > 0) {
                    if (child_cost + lb > threshold) continue;
                }
                State* child_state = new State(state);
                child_state->cost = child_cost;
                child_state->vertex_label_bound = vertex_lower_bound;
                child_state->inner_edge_label_bound = inner_edge_lower_bound;
                child_state->cross_edge_label_bound = cross_edge_lower_bound;
                child_state->ComputeLowerBound();
                child_state->mapping[u] = v;
                child_state->inverse_mapping[v] = u;
                child_state->next_mapping_order = state->next_mapping_order;
                if (depth == G1->GetNumVertices() - 1) {
                    current_best = std::min(current_best, child_state->lower_bound);
                    memcpy(&current_best_mapping[0], child_state->mapping, sizeof(int) * NumG1Vertices);
                    continue;
                }
                if (DEBUG) {
                    fprintf(stderr,"State %d pushes state %d\n", state->id, child_state->id);
                    child_state->Print();
                }
                num_pushed_nodes++;
                queue.push(child_state);
            }
        }

        void ExtendHungarianToLeaf(State *state) {
            if (DEBUG) {
                fprintf(stderr, "ExtendHungarianToLeaf(%d)\n", state->id);
                state->Print();
            }
            DifferenceVector diff; diff.init(20);
            int remaining = G2->GetNumVertices() - (state->depth + 1);
            std::vector<int> rem_left, rem_right;
            std::vector<std::vector<int>> branch_distance_matrix(remaining,
                                                                 std::vector<int>(remaining, 0));
            ComputeBranchDistanceMatrix(state, branch_distance_matrix, rem_left, rem_right);
            num_hungarian[rem_left.size()]++;
            Hungarian hungarian(branch_distance_matrix);
            hungarian.Solve();
            auto &assignment = hungarian.GetAssignment();

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
                //return hungarian_weight[i] < hungarian_weight[j];
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
                ExtendStateLSa(state, rem_left[0], hungarian_mapping[rem_left[0]]);
            for (int i = 0; i < rem_left.size(); i++) {
                int u = rem_left[i];
                int v = hungarian_mapping[u];
                state->cost = GetChildEditCost(state, u, v);
                hcost -= hungarian_weight[u];
                state->lower_bound = std::max(state->lower_bound, (hcost + 1) / 2 + state->cost);
//                state->lower_bound = std::max(lb, state->cost);
                state->mapping[u] = v;
                state->inverse_mapping[v] = u;
                state->depth++;
//                fprintf(stderr, "State %d of %d-%d became a state with cost %d and lb %d at depth %d\n", state->id,
//                        G1->GetId(), G2->GetId(), state->cost, state->lower_bound, state->depth);
                if (i + 1 < rem_left.size()) {
                    int nxt_u = rem_left[i + 1];
                    int nxt_v = hungarian_mapping[nxt_u];
                    if (state->lower_bound < current_best)
                        ExtendStateLSa(state, nxt_u, nxt_v);
                    else return;
                }
            }
            if (DEBUG) {
                fprintf(stderr, "We reached to the leaf!\n");
                state->Print();
            }
            current_best = std::min(current_best, ComputeDistance(state));
//            if (G1->GetId() == 6707 && G2->GetId() == 5025) exit(3);
        }

        int AStar() {
            global_state_id = 0;
            num_nodes = 0;
            num_pushed_nodes = 0;
            current_best = 1e9;
            inv_matching_order.resize(G1->GetNumVertices(), -1);
            for (int i = 0; i < G1->GetNumVertices(); i++) {
                inv_matching_order[matching_order[i]] = i;
            }
            State* initial_state = new State(NULL);
            initial_state->cost = 0;
            initial_state->depth = -1;
            ExtendHungarianToLeaf(initial_state);
            int64_t max_qsize = 1;
            while (!queue.empty()) {
                State* current_state = queue.top();
                num_nodes++;
                if (num_nodes % LOG_EVERY == 0) {
                    fprintf(stderr, "GED of %d and %d\n", G1->GetId(), G2->GetId());
                    fprintf(stderr,"%llu, Current best = %d, current cost = %d / lb = %d, depth %d, Queuesize %lu\n", num_nodes, current_best,
                            current_state->cost, current_state->lower_bound, current_state->depth, queue.size());
                }
                queue.pop();
                if (current_state->lower_bound >= current_best) {
                    queue = std::priority_queue<State*, std::vector<State*>, StateComparator>();
                    break;
                }
                if (threshold >= 0) {
                    if (current_best < threshold) {
                        queue = std::priority_queue<State*, std::vector<State*>, StateComparator>();
                        break;
                    }
                    if (current_state->lower_bound > threshold) {
                        current_best = -1;
                        queue = std::priority_queue<State*, std::vector<State*>, StateComparator>();
                        break;
                    }
                }
                ExtendHungarianToLeaf(current_state);
                max_qsize = std::max(max_qsize, (int64_t)queue.size());
            }
            if (threshold >= 0 and current_best > threshold) {current_best = -1;}
            log.AddResult("MaxQueueSize",max_qsize, RESULT_INT64);
            log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
            log.AddResult("AStarPushedNodes", num_pushed_nodes, RESULT_INT64);
            log.AddResult("EditDistance", current_best, RESULT_INT);
            return current_best;
        }

    };
}