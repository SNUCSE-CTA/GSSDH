#pragma once
/*
 * LabeledGraph database for managing collection of small to medium-sized graphs.
 */

#include "Base/BasicAlgorithms.h"
#include "DataStructure/LabeledGraph.h"
#include "Branch.h"
#include "DifferenceVector.h"
#include "Base/Hungarian.h"
#include "GraphSimilarity/GSSEntry.h"
#include "GraphSimilarity/EditDistance.h"
#include "GraphSimilarity/PartitionFilter.h"

namespace GraphLib::GraphSimilarity {

    class GraphSimilaritySearch {
        double total_filtering_time = 0.0, total_verifying_time = 0.0;
        std::vector<GSSEntry*> data_graphs;
        GSSEntry *query = nullptr;
        std::vector<std::vector<int>> indexed_branch_edit_distance;
        std::map<Branch, int> seen_branch_ids;
        std::vector<Branch> seen_branch_structures;
        std::unordered_map<std::string, int> global_vertex_labels, global_edge_labels;
        int num_global_vertex_labels = 1, num_global_edge_labels = 1, tau = 0, num_indexed_branches;
        ResultLogger log;
        DifferenceVector vlabel_diff, elabel_diff;
        Hungarian *hungariansolver = nullptr;

        EditDistanceSolver GEDSolver;

        int num_answer = 0;
        std::vector<ResultLogger> ged_logs;
    public:
        int total_candidates = 0, total_filtered = 0;
        ResultLogger GetLog() {return log;}
        void LoadGraphDatabase(std::string &filename, int which);
        void RetrieveSimilarGraphs(GSSEntry* query, int tau);
        std::vector<GSSEntry*>& GetData() {return data_graphs;}
        int GetNumGlobalVertexLabels() {return num_global_vertex_labels;}
        int GetNumGlobalEdgeLabels() {return num_global_edge_labels;}
        void BuildBranches();
        GraphSimilaritySearch() {
            hungariansolver = new Hungarian(100);
        }
        ~GraphSimilaritySearch() {
            for (auto entry : data_graphs) {
                delete entry;
            }
            data_graphs.clear();
        }

        /**
         * Functions for GED Filtering
         */

        int NaiveCountBound(int data_idx) {
            GSSEntry *data = data_graphs[data_idx];
            int bound = abs(data->GetNumVertices() - query->GetNumVertices()) + abs(data->GetNumEdges() - query->GetNumEdges());
            return bound;
        }

        int DegreeSequenceBound(int data_idx) {
            GSSEntry *data = data_graphs[data_idx];
            int pos = 0, neg = 0;
            auto &g1_deg = data->GetDegreeSequence();
            auto &g2_deg = query->GetDegreeSequence();
            for (int i = 0; i < query->GetNumVertices(); i++) {
                int a = (i >= data->GetNumVertices()) ? 0 : g1_deg[i];
                int b = g2_deg[i];
                if (a > b) pos += (a - b);
                if (a < b) neg += (b - a);
            }
            return (pos + 1) / 2 + (neg + 1) / 2;
        }

        int LabelSetDifferenceBound(int data_idx) {
            GSSEntry *data = data_graphs[data_idx];
            vlabel_diff.reset();
            elabel_diff.reset();
            for (auto &[l, t] : data->GetVertexLabelFrequency())
                vlabel_diff.update(l, t);
            for (auto &[l, t] : query->GetVertexLabelFrequency())
                vlabel_diff.update(l, -t);
            for (auto &[l, t] : data->GetEdgeLabelFrequency())
                elabel_diff.update(l, t);
            for (auto &[l, t] : query->GetEdgeLabelFrequency())
                elabel_diff.update(l, -t);
            int cost = 0;
            cost = vlabel_diff.GetDifference() + elabel_diff.GetDifference();
            return cost;
        }

        int BranchBound(int data_idx);
        int PartitionBound(int data_idx);
    };

    void GraphSimilaritySearch::LoadGraphDatabase(std::string &filename, int opt) {
        std::ifstream fin(filename);
        std::vector<int> graph_ids;
        std::vector<std::vector<int>> graph_vertices, graph_vertex_labels, graph_edge_labels;
        std::vector<std::vector<std::pair<int, int>>> graph_edges;
        while (!fin.eof()) {
            std::string raw;
            std::getline(fin, raw);
            auto line = parse(raw, " ");
            if (line[0] == "t") {
                graph_ids.push_back(std::stoi(line[2]));
                graph_vertices.emplace_back();
                graph_vertex_labels.emplace_back();
                graph_edge_labels.emplace_back();
                graph_edges.emplace_back();
            }
            else if (line[0] == "v") {
                int id, label;
                id = std::stoi(line[1]);
                string raw_label = line[2];
                if (global_vertex_labels.find(raw_label) == global_vertex_labels.end()) {
                    global_vertex_labels[raw_label] = num_global_vertex_labels++;
                }
                label = global_vertex_labels[raw_label];
                graph_vertices.back().push_back(id);
                graph_vertex_labels.back().push_back(label);
            }
            else if (line[0] == "e") {
                int u, v, label;
                u = std::stoi(line[1]);
                v = std::stoi(line[2]);
                string raw_label = line[3];
                if (global_edge_labels.find(raw_label) == global_edge_labels.end()) {
                    global_edge_labels[raw_label] = num_global_edge_labels++;
                }
                label = global_edge_labels[raw_label];
                graph_edges.back().emplace_back(u, v);
                graph_edge_labels.back().push_back(label);
            }
        }
        for (int i = 0; i < graph_ids.size(); i++) {
            GSSEntry *G = new GSSEntry;
            G->LoadFromDatabase(graph_ids[i], graph_vertices[i], graph_vertex_labels[i],
                                 graph_edges[i], graph_edge_labels[i]);
            if (opt == -1)
                data_graphs.emplace_back(G);
            else {
                RetrieveSimilarGraphs(G, opt);
            }
        }
        if (opt >= 0) {
            log.AddResult("NUM_CANDIDATES", total_candidates, RESULT_INT);
            log.AddResult("NUM_FILTERED", total_filtered, RESULT_INT);
            log.AddResult("FILTERING_TIME", total_filtering_time, RESULT_DOUBLE_FIXED);
            log.AddResult("VERIFY_TIME", total_verifying_time, RESULT_DOUBLE_FIXED);
            log.AddResult("Ans", num_answer, RESULT_INT);
            log.AddResult("TotalSearchSpace", (int64_t)Total(ged_logs, "AStarNodes") , RESULT_INT64);
            log.AddResult("TotalMaxQueueSize", (int64_t)Total(ged_logs, "MaxQueueSize"), RESULT_INT64);
        }
    }

    void GraphSimilaritySearch::BuildBranches() {
        for (int i = 0; i < data_graphs.size(); i++) {
            data_graphs[i]->BuildBranches(seen_branch_ids, seen_branch_structures);
//            fprintf(stderr, "Datagraph [%d] has %lu branches\n", i, data_graphs[i]->GetBranchIDs().size());
        }
        log.AddResult("DATA_NUM_BRANCHES", (int)seen_branch_structures.size(), RESULT_INT);
        indexed_branch_edit_distance.resize(seen_branch_structures.size(), std::vector<int>(seen_branch_structures.size(), 0));
        for (int i = 0; i < seen_branch_structures.size(); i++) {
            for (int j = i+1; j < seen_branch_structures.size(); j++) {
                indexed_branch_edit_distance[i][j] = (BranchEditDistance(seen_branch_structures[i], seen_branch_structures[j]));
                indexed_branch_edit_distance[j][i] = indexed_branch_edit_distance[i][j];
            }
        }
        num_indexed_branches = indexed_branch_edit_distance.size();
    }

    void GraphSimilaritySearch::RetrieveSimilarGraphs(GSSEntry* query_, int tau_) {
        this->query = query_;
        this->tau = tau_;
        int num_filtered = 0, num_candidates = 0;
        vlabel_diff.init(100); elabel_diff.init(20);
        for (int data_idx = 0; data_idx < data_graphs.size(); data_idx++) {
            GSSEntry *data = data_graphs[data_idx];
            vlabel_diff.reset(); elabel_diff.reset();
            Timer timer; timer.Start();
            Timer verification_timer;
            int naive_count_bound = -1, degree_seq_bound = -1, label_set_bound = -1, branch_bound = -1, partition_bound = -1, ed;
            naive_count_bound = NaiveCountBound(data_idx);
            if (naive_count_bound > tau) goto filtered;
            degree_seq_bound  = DegreeSequenceBound(data_idx);
            if (degree_seq_bound > tau) goto filtered;
            label_set_bound   = LabelSetDifferenceBound(data_idx);
            if (label_set_bound > tau) goto filtered;
            branch_bound      = BranchBound(data_idx);
            if (branch_bound > tau) goto filtered;
            timer.Stop(); total_filtering_time+=timer.GetTime();

            verify:
            num_candidates++;
            verification_timer.Start();
//            fprintf(stdout, "CandPair %d %d\n", query->GetId(), data->GetId());
            GEDSolver.Initialize(query, data, this->tau, &vlabel_diff, &elabel_diff);
            if (GEDSolver.AStar() != -1) {
                num_answer++;
            }
            ged_logs.push_back(GEDSolver.GetLog());
            verification_timer.Stop(); total_verifying_time+=verification_timer.GetTime();
            continue;

            filtered:
            timer.Stop(); total_filtering_time+=timer.GetTime();
            num_filtered++;
            continue;

        }
        total_candidates += num_candidates;
        total_filtered += num_filtered;
        fprintf(stderr, "Filtered %d out of %lu graphs\n", num_filtered, data_graphs.size());
    }

    int GraphSimilaritySearch::BranchBound(int data_idx) {
        GSSEntry *data = data_graphs[data_idx];
        query->BuildBranches(seen_branch_ids, seen_branch_structures);
        int max_num_vertices = std::max(data->GetNumVertices(), query->GetNumVertices());
        std::vector<std::vector<int>> cost_matrix(max_num_vertices, std::vector<int>(max_num_vertices, 0));
        for (int i = 0; i < query->GetNumVertices(); i++) {
            int query_branch_id = query->GetBranchID(i);
            for (int j = 0; j < data->GetNumVertices(); j++) {
                int data_branch_id = data->GetBranchID(j);
                if (query_branch_id < num_indexed_branches and data_branch_id < num_indexed_branches) {
                    cost_matrix[i][j] = indexed_branch_edit_distance[query_branch_id][data_branch_id];
                }
                else {
                    cost_matrix[i][j] = BranchEditDistance(seen_branch_structures[query_branch_id],
                                                           seen_branch_structures[data_branch_id]);
                }
            }
            if (data->GetNumVertices() < max_num_vertices) {
                int null_distance = BranchEditDistanceFromNull(seen_branch_structures[query_branch_id]);
                for (int j = data->GetNumVertices(); j < max_num_vertices; j++) {
                    cost_matrix[i][j] = null_distance;
                }
            }
        }
        if (query->GetNumVertices() < max_num_vertices) {
            for (int j = 0; j < data->GetNumVertices(); j++) {
                int data_branch_id = data->GetBranchID(j);
                int null_distance = BranchEditDistanceFromNull(seen_branch_structures[data_branch_id]);
                for (int i = query->GetNumVertices(); i < max_num_vertices; i++) {
                    cost_matrix[i][j] = null_distance;
                }
            }
        }
        Hungarian hungarian(cost_matrix);
        hungarian.Solve();
        int cost = (hungarian.GetTotalCost() + 1) / 2;
        return cost;
    }

    int GraphSimilaritySearch::PartitionBound(int data_idx) {
        GSSEntry *data = data_graphs[data_idx];
        PartitionFilter partition_filter(query, data);
        return partition_filter.GetPartitionBound();
    }
}