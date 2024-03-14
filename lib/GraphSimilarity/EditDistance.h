#pragma once 
#include "DataStructure/LabeledGraph.h"
#include "DataStructure/LabeledGraphDatabase.h"
#include "WeisfeilerLehman.h"
#include "Base/Logger.h"
#include "Base/Timer.h"
#include "State.h"

struct DifferenceVector {
    std::vector<int> value;
    int pos = 0, neg = 0;
    DifferenceVector(int sz = 0) {value.resize(sz, 0);}

    void init(int sz) {value.resize(sz, 0); pos = neg = 0;}
    void reset() { std::fill(value.begin(), value.end(), 0); pos = neg = 0;}

    unsigned long size() const {return value.size();}

    int& operator[](int index) {
        return value[index];
    }

    // val should be +1 or -1
    void update(int idx, int val) {
        if (val > 1) {
            for (int i = 0; i < val; i++) update(idx, 1);
            return;
        }
        if (val < -1) {
            for (int i = 0; i < -val; i++) update(idx, -1);
            return;
        }
        if (value[idx] == 0) {
            (val > 0 ? pos : neg)++;
        }
        else if (value[idx] < 0) {
            neg += (val > 0 ? -1 : 1);
        }
        else if (value[idx] > 0) {
            pos += (val > 0 ? 1 : -1);
        }
        value[idx] += val;
    }

    int GetSetDifferenceDistance() {
        return std::max(pos, neg);
    }
};

namespace GraphLib {
bool verbosity = true;
static int32_t LOG_EVERY = 50000;

static int num_wl_iterations = 0;
class EditDistanceSolver {
    WeisfeilerLehman *WL = nullptr;
    LabeledGraphDatabaseEntry *G1, *G2;
    int threshold = -1, current_best = 1e9;
    long long num_nodes = 0;
    std::vector<int> current_best_mapping;
    std::priority_queue<State*, std::vector<State*>, StateComparator> queue;
    DifferenceVector vertex_labels, edge_labels;
    std::vector<int> matching_order;

    ResultLogger log;
public:

    int GetCurrentBestGED() const {return current_best;}
    ResultLogger GetLog() {return log;}

    LabeledGraphDatabaseEntry *GCombined = nullptr;

    EditDistanceSolver(LabeledGraphDatabase &DB) {
        vertex_labels.init(DB.GetNumGlobalVertexLabels()+1);
        edge_labels.init(DB.GetNumGlobalEdgeLabels()+1);
        printf("Intitializing edit distance solver...#Labels = (%lu, %lu)\n",vertex_labels.size(), edge_labels.size());
    }

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

    bool Initialize(LabeledGraphDatabaseEntry *G1_, LabeledGraphDatabaseEntry *G2_, int threshold_ = -1) {
        log.clear();
        Timer timer; timer.Start();
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

        this->threshold = threshold_;
        if (this->threshold >= 0) {
            if (!Verify()) {
                timer.Stop();
                log.AddResult("FilteringTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
                return false;
            }
        }
        NumG1Vertices = G1->GetNumVertices();
        NumG2Vertices = G2->GetNumVertices();
        matching_order.clear();
        ComputeMatchingOrder();
        // fprintf(stderr,"Compute Edit Distance: %d-%d vertex\n", this->G1->GetNumVertices(), this->G2->GetNumVertices());
        //
        // fprintf(stderr,"MatchingOrder: ");
        // for (auto x : matching_order) {
        //     fprintf(stderr,"%d ", x);
        // }
        // fprintf(stderr,"\n");
        current_best_mapping.clear();
        current_best_mapping.resize(G1->GetNumVertices(), -1);
        vertex_labels.reset();
        edge_labels.reset();
        num_nodes = 0;
        if (num_wl_iterations > 0) {
            delete WL;
            delete GCombined;
            WL = nullptr;
            GCombined = nullptr;

            GCombined = new LabeledGraphDatabaseEntry();
            GCombined->CombineGraph(G1, G2);
            WL = new WeisfeilerLehman(this->GCombined);
            WL->Reset(current_best_mapping);
//            for (int i = 0; i < G1->GetNumVertices(); i++) {
//                printf("We start from %d->%d\n", i, current_best_mapping[i]);
//            }
            WL->Refine(2);
            WL->Match(current_best_mapping);
//            for (int x : current_best_mapping) {
//                fprintf(stderr, "%d ", x);
//            }
//            fprintf(stderr, "\n");
            current_best = ComputeDistance(current_best_mapping, false);
//            fprintf(stderr, "color-distance: %d\n", current_best);
        }
        else {
            std::iota(current_best_mapping.begin(), current_best_mapping.end(), 0);
            current_best = ComputeDistance(current_best_mapping, false);
        }
        timer.Stop();
        log.AddResult("FilteringTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
        return true;
//        exit(0);
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


    int DegreeSequenceBound() {
        int pos = 0, neg = 0;
        auto &g1_deg = G1->GetDegreeSequence();
        auto &g2_deg = G2->GetDegreeSequence();
        for (int i = 0; i < G2->GetNumVertices(); i++) {
            int a = (i >= G1->GetNumVertices()) ? 0 : g1_deg[i];
            int b = g2_deg[i];
            if (a > b) pos += (a - b);
            if (a < b) neg += (b - a);
        }
        return (pos + 1) / 2 + (neg + 1) / 2;
    }

    int LabelSetLowerBound() {
        vertex_labels.reset();
        edge_labels.reset();
        for (auto &[l, t] : G1->GetVertexLabelFrequency())
            vertex_labels.update(l, t);
        for (auto &[l, t] : G2->GetVertexLabelFrequency())
            vertex_labels.update(l, -t);
        for (auto &[l, t] : G1->GetEdgeLabelFrequency())
            edge_labels.update(l, t);
        for (auto &[l, t] : G2->GetEdgeLabelFrequency())
            edge_labels.update(l, -t);
        int cost = 0;
        cost = vertex_labels.GetSetDifferenceDistance() + edge_labels.GetSetDifferenceDistance();
        // fprintf(stderr, "LS Bound: %d\n", cost);
        return cost;
    }

    int NaiveLowerBound() {
        int bound = abs(G1->GetNumVertices() - G2->GetNumVertices()) + abs(G1->GetNumEdges() - G2->GetNumEdges());
        // fprintf(stderr, "Naive Bound: %d\n", bound);
        return bound;
    }

    bool Verify() {
        if (NaiveLowerBound() > threshold) return false;
        if (DegreeSequenceBound() > threshold) return false;
        if (LabelSetLowerBound() > threshold) return false;
        return true;
    }

    int flip(std::vector<int>& v, int idx, int dir) {
        int t = 0;
        t -= abs(v[idx]);
        v[idx] += dir;
        t += abs(v[idx]);
//        return 2 * (dir * v[idx] > 0) - 1;
        return t;
    }

    void ExtendState(State *state) {
        if (state->cost >= current_best) return;
        int depth = state->depth + 1;
        int u = matching_order[depth];

        if (num_wl_iterations > 0) {
            std::vector<int> current_mapping(G1->GetNumVertices(), -1);
            for (int i = 0; i < G1->GetNumVertices(); i++) {
                current_mapping[i] = state->mapping[i];
            }
            WL->Reset(current_mapping);
            WL->Refine(num_wl_iterations);
            WL->Match(current_mapping);
            int ub = ComputeDistance(current_mapping, false);
            if (ub == state->lower_bound)
                return;
            if (ub < current_best) {
                current_best = ub;
                current_best_mapping = current_mapping;
            }
        }


        DifferenceVector unmapped_vertex_labels(vertex_labels.size());
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            if (state->mapping[i] == -1)
                unmapped_vertex_labels.update(G1->GetVertexLabel(i), +1);
        }
        for (int i = 0; i < G2->GetNumVertices(); i++) {
            if (state->inverse_mapping[i] == -1)
                unmapped_vertex_labels.update(G2->GetVertexLabel(i), -1);
        }

        DifferenceVector unmapped_inner_edge_labels(edge_labels.size());
        std::vector<DifferenceVector> unmapped_cross_edge_labels(G1->GetNumVertices(), DifferenceVector(edge_labels.size()));
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
            vertex_lower_bound = unmapped_vertex_labels.GetSetDifferenceDistance();
            inner_edge_lower_bound = unmapped_inner_edge_labels.GetSetDifferenceDistance();
            cross_edge_lower_bound = 0;
            for (auto &it : unmapped_cross_edge_labels) {
                cross_edge_lower_bound += it.GetSetDifferenceDistance();
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
        Timer timer; timer.Start();
        State* initial_state = new State(NULL);
        initial_state->cost = 0;
        initial_state->vertex_label_bound = 0;
        initial_state->vertex_label_bound = vertex_labels.GetSetDifferenceDistance();
        initial_state->inner_edge_label_bound = edge_labels.GetSetDifferenceDistance();
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
        timer.Stop();
        log.AddResult("MaxQueueSize",max_qsize, RESULT_INT64);
        log.AddResult("AStarTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
        log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
        log.AddResult("EditDistance", current_best, RESULT_INT);
        return current_best;
    }
};

}