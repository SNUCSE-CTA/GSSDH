#pragma once
#include "DataStructure/LabeledGraph.h"
#include "GraphSimilaritySearch.h"
#include "WeisfeilerLehman.h"
#include "Base/Logger.h"
#include "Base/Timer.h"
#include "State.h"
#include "DifferenceVector.h"


namespace GraphLib::GraphSimilarity {
bool verbosity = true;
static int32_t LOG_EVERY = 50000;

class EditDistanceSolver {
    GSSEntry *G1, *G2;
    int threshold = -1, current_best = 1e9;
    long long num_nodes = 0;
    std::vector<int> current_best_mapping;
    std::priority_queue<State*, std::vector<State*>, StateComparator> queue;
    DifferenceVector *vlabel_diff, *elabel_diff;
    std::vector<int> matching_order;
    ResultLogger log;
public:

    int GetCurrentBestGED() const {return current_best;}
    ResultLogger GetLog() {return log;}

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
                if (T[x] == 1) continue;
                T[x] = 1;
                weighted_queue.push({w[x], x});
            }
        }
        for (int i = 0; i < N; i++) {
            if (T[i] != 1) {
                matching_order.push_back(i);
            }
        }
    }


    void Initialize(GSSEntry *G1_, GSSEntry *G2_, int threshold_ = -1,
                    DifferenceVector *vlabel_diff_ = nullptr, DifferenceVector *elabel_diff_ = nullptr) {
        log.clear();
        current_best = 1e9;
        // V(G1) <= V(G2)
        if (G1_->GetNumVertices() > G2_->GetNumVertices()) {
            this->G1 = G2_;
            this->G2 = G1_;
        }
        else {
            this->G1 = G1_;
            this->G2 = G2_;
        }
        vlabel_diff = vlabel_diff_;
        elabel_diff = elabel_diff_;

        this->threshold = threshold_;
        NumG1Vertices = G1->GetNumVertices();
        NumG2Vertices = G2->GetNumVertices();
        matching_order.clear();
        ComputeMatchingOrder();
        current_best_mapping.clear();
        current_best_mapping.resize(G1->GetNumVertices(), -1);
        num_nodes = 0;
    }

    void ExtendState(State *state) {
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

//        printf("Unmapped Vertex Label Distribution: ");
//        for (int i = 0; i < vertex_labels.size(); i++) {
//            printf("%d ", unmapped_vertex_labels[i]);
//        }
//        printf("\n");

        int num_u_edges = 0;
        for (int u_nbr : G1->GetNeighbors(u)) {
            if (state->mapping[u_nbr] != -1)
                num_u_edges++;
        }

        for (int v = 0; v < G2->GetNumVertices(); v++) {
            if (state->inverse_mapping[v] != -1) continue;

            int child_cost = state->cost;
            int vertex_lower_bound = state->vertex_label_bound;
            int inner_edge_lower_bound = state->inner_edge_label_bound;
            int cross_edge_lower_bound = state->cross_edge_label_bound;

            int u_label = G1->GetVertexLabel(u), v_label = G2->GetVertexLabel(v);
            // update editorial cost
            if (u_label != v_label) {
                child_cost++;
            }
            int ec = num_u_edges;
            for (int vprime : G2->GetNeighbors(v)) {
                int uprime = state->inverse_mapping[vprime];
                if (uprime == -1) continue;
                ec++;
                int l1 = G1->GetEdgeLabel(u, uprime);
                int l2 = G2->GetEdgeLabel(v, vprime);
                if (l1 == -1) continue;
                if (l1 == l2) ec -= 2;
                else ec--;
            }
            child_cost += ec;

            // Update Vertex-label-based Lower Bound
            unmapped_vertex_labels.update(u_label, -1);
            unmapped_vertex_labels.update(v_label, +1);
//            vertex_lower_bound += flip(unmapped_vertex_labels, u_label, -1);
//            vertex_lower_bound += flip(unmapped_vertex_labels, v_label, +1);

            // Update Inner-edge-label-based Lower Bound
            for (int uprime : G1->GetNeighbors(u)) {
                int el = G1->GetEdgeLabel(u, uprime);
//                if (verbosity) printf("Consider %d (label = %d)\n",uprime, el);
                if (state->mapping[uprime] == -1) {
//                    inner_edge_lower_bound += flip(unmapped_inner_edge_labels, el, -1);
//                    cross_edge_lower_bound += flip(unmapped_cross_edge_labels[u], el, +1);
                    unmapped_inner_edge_labels.update(el, -1);
                    unmapped_cross_edge_labels[u].update(el, +1);
                }
                else {
                    // uPrime is matched: anchored on uprime ->removed
                    unmapped_cross_edge_labels[uprime].update(el, -1);
//                    cross_edge_lower_bound += flip(unmapped_cross_edge_labels[uprime], el, -1);
                }
//                if (verbosity) printf("Inner is %d, Cross is %d\n",inner_edge_lower_bound, cross_edge_lower_bound);
            }
            for (int vprime : G2->GetNeighbors(v)) {
                int el = G2->GetEdgeLabel(v, vprime);
                int fvprime = state->inverse_mapping[vprime];
//                if (verbosity) printf("Consider (%d-%d) G2 (label = %d)\n",v,vprime,el);
                if (fvprime == -1) {
                    unmapped_inner_edge_labels.update(el, +1);
                    unmapped_cross_edge_labels[u].update(el, -1);
//                    inner_edge_lower_bound += flip(unmapped_inner_edge_labels, el, +1);
//                    cross_edge_lower_bound += flip(unmapped_cross_edge_labels[u], el, -1);
                }
                else {
                    // uPrime is matched: anchored on uprime ->removed
                    unmapped_cross_edge_labels[fvprime].update(el, +1);
//                    cross_edge_lower_bound += flip(unmapped_cross_edge_labels[fvprime], el, +1);
                }
//                if (verbosity) printf("Inner is %d, Cross is %d\n",inner_edge_lower_bound, cross_edge_lower_bound);
            }
            vertex_lower_bound = unmapped_vertex_labels.GetDifference();
            inner_edge_lower_bound = unmapped_inner_edge_labels.GetDifference();
            cross_edge_lower_bound = 0;
            for (auto &it : unmapped_cross_edge_labels) {
                cross_edge_lower_bound += it.GetDifference();
            }
            int lb = vertex_lower_bound + inner_edge_lower_bound + cross_edge_lower_bound;

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
                }
                else {
                    unmapped_cross_edge_labels[uprime].update(el, +1);
//                    flip(unmapped_cross_edge_labels[uprime], el, +1);
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
//            printf("State [%p] Check Depth %d [Parent %p]! matching (%d(%d), %d(%d)) gives lb = %d, cost = %d (curbest = %d)\n", child_state, depth, state, u, u_label, v, v_label, lb+child_cost, child_cost, current_best);
//            printf("  LB Breakdown: VLabel %d, InnerEdge %d, CrossEdge %d\n\n", vertex_lower_bound, inner_edge_lower_bound, cross_edge_lower_bound);
//            if (depth == 2) exit(1);
            child_state->cost = child_cost;
            child_state->vertex_label_bound = vertex_lower_bound;
            child_state->inner_edge_label_bound = inner_edge_lower_bound;
            child_state->cross_edge_label_bound = cross_edge_lower_bound;
            child_state->ComputeLowerBound();
            child_state->mapping[u] = v;
            child_state->inverse_mapping[v] = u;
            if (depth == G1->GetNumVertices() - 1) {
                current_best = std::min(current_best, child_state->lower_bound);
                memcpy(&current_best_mapping[0], child_state->mapping, sizeof(int) * NumG1Vertices);
                continue;
            }
//            fprintf(stdout,"Push state %p with bound %d\n", child_state, child_state->lower_bound);
            queue.push(child_state);
        }
    }

    int AStar() {
        num_nodes = 0;
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
            ExtendState(current_state);
            max_qsize = std::max(max_qsize, (int64_t)queue.size());
        }
        if (threshold >= 0 and current_best > threshold) {current_best = -1;}
        log.AddResult("MaxQueueSize",max_qsize, RESULT_INT64);
        log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
        log.AddResult("EditDistance", current_best, RESULT_INT);
        return current_best;
    }



    int ComputeDistance(std::vector<int>& mapping, bool verbose=false) {
        std::vector<int> inverse_mapping(G2->GetNumVertices(), -1);
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            inverse_mapping[mapping[i]] = i;
        }
        return ComputeDistance(mapping, inverse_mapping, verbose);
    }

    int ComputeDistance(std::vector<int>& mapping, std::vector<int>& inverse_mapping, bool verbose = false) {
        int cost = 0;
        // vertex re-labeling cost
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            int l1 = G1->GetVertexLabel(i);
            int l2 = G2->GetVertexLabel(mapping[i]);
            if (l1 != l2) {
                if (verbose)
                    printf("Vertex %d(%d)-%d(%d) re-labeling cost!\n", i, l1, mapping[i], l2);
                cost++;
            }
        }
        if (verbose)
            printf("#Missing Vertices: %d\n",(G2->GetNumVertices() - G1->GetNumVertices()));
        cost += (G2->GetNumVertices() - G1->GetNumVertices());
        for (auto &[u, v] : G1->GetEdges()) {
            int fu = mapping[u], fv = mapping[v];
            int l1 = G1->GetEdgeLabel(u, v);
            int l2 = G2->GetEdgeLabel(fu, fv);
            if (l1 != l2) {
                if (verbose)
                    printf("Edge (%d, %d)(%d)-(%d, %d)(%d) re-labeling cost!\n",
                           u,v,l1,fu,fv,l2);
                cost++;
            }
        }
        for (auto &[u, v] : G2->GetEdges()) {
            int inv_u = inverse_mapping[u], inv_v = inverse_mapping[v];
            if (inv_u == -1 || inv_v == -1) {
                if (verbose)
                    printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n",
                           u, v, inv_u, inv_v);
                cost++;
            }
            else {
                int l = G1->GetEdgeLabel(inv_u, inv_v);
                if (l == -1) {
                    if (verbose)
                        printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n",
                               u, v, inv_u, inv_v);
                    cost++;
                }
            }
        }
        if (verbose)
            printf("Total ED Cost: %d\n",cost);
        current_best = std::min(cost, current_best);
        return cost;
    }

};

}