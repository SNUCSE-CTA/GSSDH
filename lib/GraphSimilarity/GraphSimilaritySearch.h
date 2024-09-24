#pragma once
/*
 * LabeledGraph database for managing collection of small to medium-sized
 * graphs.
 */

// #include "Base/BasicAlgorithms.h"
// #include "Base/Hungarian.h"
#include "Branch.h"
#include "DataStructure/LabeledGraph.h"
#include "DifferenceVector.h"
#include "GraphSimilarity/EditDistance.h"
// #include "GraphSimilarity/GSSEntry.h"
// #include "GraphSimilarity/GraphEditDistance/AStarBMa.h"
#include "GraphSimilarity/PartitionFilter.h"
// #include "GraphSimilarity/GraphEditDistance/AStarMixed.h"
// #include "GraphSimilarity/GraphEditDistance/OurGED.h"
#include "Base/DynamicHungarian.h"
#ifdef CC
// #include "GraphSimilarity/GraphEditDistance/AStarDH.h"
// #include "GraphSimilarity/GraphEditDistance/AStarDHo.h"
#include "GraphSimilarity/GraphEditDistance/DFSDH.h"
#endif
#ifdef COLORING
#include "GraphSimilarity/GraphEditDistance/AStarLSa.h"
#endif

namespace GraphLib::GraphSimilarity {
class GraphSimilaritySearch {
  double total_filtering_time = 0.0, total_verifying_time = 0.0;
  std::vector<GSSEntry *> data_graphs, query_graphs;
  GSSEntry *query = nullptr;
  std::vector<std::vector<int>> indexed_branch_edit_distance;
  std::map<Branch, int> seen_branch_ids;
  std::vector<Branch> seen_branch_structures;
  std::unordered_map<std::string, int> global_vertex_labels, global_edge_labels;
  int num_global_vertex_labels = 1, num_global_edge_labels = 1, tau = 0,
      num_indexed_branches;
  ResultLogger log;
#ifdef CC
  Hungarian *hungariansolver = nullptr;
  // AStarDHo GEDSolver;
  // AStarDH GEDSolver;
  DFSDH GEDSolver;
#endif
#ifdef COLORING
  AStarLSa GEDSolver;
#endif
  // AStarBMa GEDSolver;
  double total_hg_time = 0.0, total_bd_time = 0.0;
  int64_t total_hungarian_vertex_num = 0;
  int64_t total_dfs_cnt = 0;
  int num_answer = 0;
  std::vector<ResultLogger> ged_logs;

public:
  int total_candidates = 0, total_filtered = 0;
  ResultLogger GetLog() { return log; }
  void LoadGraphDatabase(std::string &filename, int which);
  void RetrieveSimilarGraphs(GSSEntry *query, int tau, GWL *gwl,
                             ColorTree **prev_color_to_node,
                             ColorTree **curr_color_to_node);
  std::vector<GSSEntry *> &GetData() { return data_graphs; }
  int GetNumGlobalVertexLabels() { return num_global_vertex_labels; }
  int GetNumGlobalEdgeLabels() { return num_global_edge_labels; }
  void BuildBranches();
  GraphSimilaritySearch() {
#ifdef CC
    hungariansolver = new Hungarian(100);
#endif
  }
  ~GraphSimilaritySearch() {
    for (auto entry : data_graphs) {
      delete entry;
    }
    data_graphs.clear();
  }

  int BranchBound(int data_idx);
  int PartitionBound(int data_idx);

  void ProcessSimilaritySearch(int tau_) {
    this->tau = tau_;
    int q_index = 0;
#ifdef COLORING
    GWL *gwl = new GWL(nullptr);
    ColorTree **prev_color_to_node = new ColorTree *[1000];
    ColorTree **curr_color_to_node = new ColorTree *[1000];
    for (auto &q : query_graphs) {
      // fprintf(stderr, "Processing query %d\n", q_index++);
      // fprintf(stderr, "Query id: %d\n", q->GetId());
      RetrieveSimilarGraphs(q, tau, gwl, prev_color_to_node,
                            curr_color_to_node);
    }
    delete gwl;
    delete[] prev_color_to_node;
    delete[] curr_color_to_node;
#endif
#ifdef CC
    for (auto &q : query_graphs) {
      // fprintf(stderr, "Processing query %d\n", q_index++);
      // fprintf(stderr, "Query id: %d\n", q->GetId());
      RetrieveSimilarGraphs(q, tau, nullptr, nullptr, nullptr);
    }
#endif
    log.AddResult("NUM_CANDIDATES", total_candidates, RESULT_INT);
    log.AddResult("NUM_FILTERED", total_filtered, RESULT_INT);
    log.AddResult("FILTERING_TIME", total_filtering_time, RESULT_DOUBLE_FIXED);
    log.AddResult("VERIFY_TIME", total_verifying_time, RESULT_DOUBLE_FIXED);
    log.AddResult("Ans", num_answer, RESULT_INT);
    log.AddResult("TotalSearchSpace", (int64_t)Total(ged_logs, "AStarNodes"),
                  RESULT_INT64);
    log.AddResult("TotalMaxQueueSize", (int64_t)Total(ged_logs, "MaxQueueSize"),
                  RESULT_INT64);
    /*BMa time*/
    log.AddResult("HUNGARIAN_TIME", total_hg_time, RESULT_DOUBLE_FIXED);
    log.AddResult("BranchDistance_TIME", total_bd_time, RESULT_DOUBLE_FIXED);
    log.AddResult("Hungarian_Vertices", total_hungarian_vertex_num,
                  RESULT_INT64);
    log.AddResult("DFS_depth_cnt", total_dfs_cnt, RESULT_INT64);
  }

  void CombineGraphs(GSSEntry *g1, GSSEntry *g2, GSSEntry *combined) {
    if (g1->GetNumVertices() <= g2->GetNumVertices())
      combined->CombineGraph(g1, g2);
    else
      combined->CombineGraph(g2, g1);
  }
};

void GraphSimilaritySearch::LoadGraphDatabase(std::string &filename, int opt) {
  std::ifstream fin(filename);
  std::vector<int> graph_ids;
  std::vector<std::vector<int>> graph_vertices, graph_vertex_labels,
      graph_edge_labels;
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
    } else if (line[0] == "v") {
      int id, label;
      id = std::stoi(line[1]);
      string raw_label = line[2];
      if (global_vertex_labels.find(raw_label) == global_vertex_labels.end()) {
        global_vertex_labels[raw_label] = num_global_vertex_labels++;
      }
      label = global_vertex_labels[raw_label];
      graph_vertices.back().push_back(id);
      graph_vertex_labels.back().push_back(label);
    } else if (line[0] == "e") {
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
    else
      query_graphs.emplace_back(G);
  }
}

void GraphSimilaritySearch::BuildBranches() {
  for (int i = 0; i < data_graphs.size(); i++) {
    data_graphs[i]->BuildBranches(seen_branch_ids, seen_branch_structures);
    //            fprintf(stderr, "Datagraph [%d] has %lu branches\n", i,
    //            data_graphs[i]->GetBranchIDs().size());
  }
  log.AddResult("DATA_NUM_BRANCHES", (int)seen_branch_structures.size(),
                RESULT_INT);
  indexed_branch_edit_distance.resize(
      seen_branch_structures.size(),
      std::vector<int>(seen_branch_structures.size(), 0));
  for (int i = 0; i < seen_branch_structures.size(); i++) {
    for (int j = i + 1; j < seen_branch_structures.size(); j++) {
      indexed_branch_edit_distance[i][j] = (BranchEditDistance(
          seen_branch_structures[i], seen_branch_structures[j]));
      indexed_branch_edit_distance[j][i] = indexed_branch_edit_distance[i][j];
    }
  }
  num_indexed_branches = indexed_branch_edit_distance.size();
  for (auto &q : query_graphs) {
    q->BuildBranches(seen_branch_ids, seen_branch_structures);
  }
}

void GraphSimilaritySearch::RetrieveSimilarGraphs(
    GSSEntry *query_, int tau_, GWL *gwl, ColorTree **prev_color_to_node,
    ColorTree **curr_color_to_node) {
  this->query = query_;
  this->tau = tau_;
  int num_filtered = 0, num_candidates = 0;
  double coloring_time = 0, filtering_time = 0, verifying_time = 0;
  for (int data_idx = 0; data_idx < data_graphs.size(); data_idx++) {
    // std::cout << "Processing data graph " << data_idx << std::endl;
    GSSEntry *data = data_graphs[data_idx];
    // if (data->GetId() != query->GetId())
    //   continue;
    GEDSolver.InitializeSolver(query, data, this->tau);
    Timer filtering_timer;
    filtering_timer.Start();
    bool filtering_result = GEDSolver.GEDVerificiationFiltering();
    filtering_timer.Stop();
    total_filtering_time += filtering_timer.GetTime();
    filtering_time += filtering_timer.GetTime();
    if (filtering_result) {
      num_candidates++;
      Timer verification_timer;
      verification_timer.Start();
#ifdef CC
      int ged = GEDSolver.GED();
#endif
#ifdef COLORING
      GSSEntry *combined = new GSSEntry;
      CombineGraphs(query, data, combined);
      gwl->SetGraph(combined);
      int ged =
          GEDSolver.GED(combined, gwl, prev_color_to_node, curr_color_to_node);
      delete combined;
#endif
      if (ged != -1) {
        num_answer++;
        // std::cout << "Query : "<<query->GetId() << "\n";
        // std::cout << data_idx << "\n";
      }
#ifdef CC
      /*AStarBMa time*/
      // total_hungarian_vertex_num += GEDSolver.GetVertNum();
      total_dfs_cnt += GEDSolver.GetCnt();
      total_hg_time += GEDSolver.Gethgtime();
      total_bd_time += GEDSolver.Getbdtime();
#endif
      verification_timer.Stop();
      total_verifying_time += verification_timer.GetTime();
      verifying_time += verification_timer.GetTime();
    } else {
      num_filtered++;
    }
    ged_logs.push_back(GEDSolver.GetLog());
    continue;
  }
  total_candidates += num_candidates;
  total_filtered += num_filtered;
  // fprintf(stderr, "Filtered %d out of %lu graphs\n", num_filtered,
  //         data_graphs.size());
  // fprintf(stderr, "Coloring time = %.6lf\n", coloring_time);
  // fprintf(stderr, "Filtering time = %.6lf\n", filtering_time);
  // fprintf(stderr, "Verifying time = %.6lf\n", verifying_time);
}

int GraphSimilaritySearch::BranchBound(int data_idx) {
  GSSEntry *data = data_graphs[data_idx];
  query->BuildBranches(seen_branch_ids, seen_branch_structures);
  int max_num_vertices =
      std::max(data->GetNumVertices(), query->GetNumVertices());
  std::vector<std::vector<int>> cost_matrix(
      max_num_vertices, std::vector<int>(max_num_vertices, 0));
  for (int i = 0; i < query->GetNumVertices(); i++) {
    int query_branch_id = query->GetBranchID(i);
    for (int j = 0; j < data->GetNumVertices(); j++) {
      int data_branch_id = data->GetBranchID(j);
      if (query_branch_id < num_indexed_branches and
          data_branch_id < num_indexed_branches) {
        cost_matrix[i][j] =
            indexed_branch_edit_distance[query_branch_id][data_branch_id];
      } else {
        cost_matrix[i][j] =
            BranchEditDistance(seen_branch_structures[query_branch_id],
                               seen_branch_structures[data_branch_id]);
      }
    }
    if (data->GetNumVertices() < max_num_vertices) {
      int null_distance =
          BranchEditDistanceFromNull(seen_branch_structures[query_branch_id]);
      for (int j = data->GetNumVertices(); j < max_num_vertices; j++) {
        cost_matrix[i][j] = null_distance;
      }
    }
  }
  if (query->GetNumVertices() < max_num_vertices) {
    for (int j = 0; j < data->GetNumVertices(); j++) {
      int data_branch_id = data->GetBranchID(j);
      int null_distance =
          BranchEditDistanceFromNull(seen_branch_structures[data_branch_id]);
      for (int i = query->GetNumVertices(); i < max_num_vertices; i++) {
        cost_matrix[i][j] = null_distance;
      }
    }
  }
  Hungarian hungarian(cost_matrix);
  hungarian.Solve();
  //        hungarian.Print();
  int cost = (hungarian.GetTotalCost() + 1) / 2;
  return cost;
}

int GraphSimilaritySearch::PartitionBound(int data_idx) {
  GSSEntry *data = data_graphs[data_idx];
  PartitionFilter partition_filter(query, data);
  return partition_filter.GetPartitionBound();
}
} // namespace GraphLib::GraphSimilarity
