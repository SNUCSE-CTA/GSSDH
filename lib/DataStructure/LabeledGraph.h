#pragma once
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "Base/Base.h"
#include "Graph.h"

namespace GraphLib {
// Vertex and edge labeled graph
class LabeledGraph : public Graph {
public:
  int num_vertex_labels = 0, num_edge_labels = 0;
  int combined_index = 0;
  std::vector<int> vertex_label, edge_label;
  std::vector<int> vertex_color;
  std::unordered_map<int, int> vertex_label_frequency, edge_label_frequency;

  LabeledGraph() {}
  ~LabeledGraph() {}
  LabeledGraph &operator=(const LabeledGraph &) = delete;

  inline int GetVertexLabel(int v) const { return vertex_label[v]; }

  virtual inline int GetEdgeLabel(int edge_id) const {
    return edge_label[edge_id];
  }
  inline int GetNumVertexLabels() const { return num_vertex_labels; }
  inline int GetNumEdgeLabels() const { return num_edge_labels; }

  std::unordered_map<int, int> &GetVertexLabelFrequency() {
    return vertex_label_frequency;
  }
  std::unordered_map<int, int> &GetEdgeLabelFrequency() {
    return edge_label_frequency;
  }
  int GetVertexLabelFrequency(int l) { return vertex_label_frequency[l]; }
  int GetEdgeLabelFrequency(int l) { return edge_label_frequency[l]; }
  bool HaveVertexLabel(int l) {
    return vertex_label_frequency.find(l) != vertex_label_frequency.end();
  }
  bool HaveEdgeLabel(int l) {
    return edge_label_frequency.find(l) != edge_label_frequency.end();
  }
  void CombineGraph(LabeledGraph *g1, LabeledGraph *g2) {
    combined_index = g1->GetNumVertices();
    num_vertex_labels =
        std::max(g1->GetNumVertexLabels(), g2->GetNumVertexLabels());
    num_vertex = g1->GetNumVertices() + g2->GetNumVertices();
    num_edge = g1->GetNumEdges() + g2->GetNumEdges();
    adj_list.resize(num_vertex);
    vertex_label.resize(num_vertex);
    edge_label.resize(num_edge);
    for (int i = 0; i < g1->GetNumVertices(); i++)
      vertex_label[i] = g1->GetVertexLabel(i);
    for (int i = 0; i < g2->GetNumVertices(); i++)
      vertex_label[i + g1->GetNumVertices()] = g2->GetVertexLabel(i);
    for (int i = 0; i < g1->GetNumEdges(); i++) {
      edge_label[i] = g1->GetEdgeLabel(i);
      auto [v1, v2] = g1->GetEdge(i);
      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
    }
    for (int i = 0; i < g2->GetNumEdges(); i++) {
      edge_label[i + g1->GetNumEdges()] = g2->GetEdgeLabel(i);
      auto [v1, v2] = g2->GetEdge(i);
      v1 += g1->GetNumVertices();
      v2 += g1->GetNumVertices();
      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
    }
  }

  // Read/Write graph from/to file.
  void LoadIGraph(const std::string &filename);
  void LoadCFLGraph(const std::string &filename);
  void WriteToFile(const std::string &filename);
};

void LabeledGraph::LoadIGraph(const std::string &filename) {
  std::ifstream fin(filename);
  int v;
  std::string ignore, type, line;
  fin >> ignore >> ignore >> v;
  num_vertex = v;
  adj_list.resize(num_vertex);
  vertex_label.resize(num_vertex);
  int num_lines = 0;
  while (getline(fin, line)) {
    auto tok = parse(line, " ");
    type = tok[0];
    tok.pop_front();
    if (type[0] == 'v') {
      int id = std::stoi(tok.front());
      tok.pop_front();
      int l;
      if (tok.empty())
        l = 0;
      else {
        l = std::stoi(tok.front());
        tok.pop_front();
      }
      vertex_label[id] = l;
      num_vertex_labels = std::max(num_vertex_labels, l + 1);
    } else if (type[0] == 'e') {
      int v1, v2;
      v1 = std::stoi(tok.front());
      tok.pop_front();
      v2 = std::stoi(tok.front());
      tok.pop_front();
      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
      edge_list.emplace_back(v1, v2);
      edge_list.emplace_back(v2, v1);
      int el = tok.empty() ? 0 : std::stoi(tok.front());
      edge_label.emplace_back(el);
      edge_label.emplace_back(el);
      num_edge_labels = std::max(num_edge_labels, el + 1);
      max_degree = std::max(
          max_degree, (int)std::max(adj_list[v1].size(), adj_list[v2].size()));
    }
    num_lines++;
  }
  num_edge = edge_list.size();
}

void LabeledGraph::LoadCFLGraph(const std::string &filename) {
  std::ifstream fin(filename);
  int v, e;
  std::string ignore, type, line;
  fin >> ignore >> v >> e;
  num_vertex = v;
  num_edge = e * 2;
  adj_list.resize(num_vertex);
  vertex_label.resize(num_vertex);
  edge_label.resize(num_edge);
  int num_lines = 0;
  while (getline(fin, line)) {
    auto tok = parse(line, " ");
    type = tok[0];
    tok.pop_front();
    if (type[0] == 'v') {
      int id = std::stoi(tok.front());
      tok.pop_front();
      int l;
      if (tok.empty())
        l = 0;
      else {
        l = std::stoi(tok.front());
        tok.pop_front();
      }
      vertex_label[id] = l;
      num_vertex_labels = std::max(num_vertex_labels, l + 1);
    } else if (type[0] == 'e') {
      int v1, v2;
      v1 = std::stoi(tok.front());
      tok.pop_front();
      v2 = std::stoi(tok.front());
      tok.pop_front();
      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
      edge_list.push_back({v1, v2});
      edge_list.push_back({v2, v1});
      int el = tok.empty() ? 0 : std::stoi(tok.front());
      edge_label[edge_list.size() - 2] = edge_label[edge_list.size() - 1] = el;
      max_degree = std::max(
          max_degree, (int)std::max(adj_list[v1].size(), adj_list[v2].size()));
      num_edge_labels = std::max(num_edge_labels, el + 1);
    }
    num_lines++;
  }
}

void LabeledGraph::WriteToFile(const std::string &filename) {
  std::filesystem::path filepath = filename;
  create_directories(filepath.parent_path());
  std::ofstream out(filename);
  out << "t " << GetNumVertices() << ' ' << GetNumEdges() / 2 << '\n';
  for (int i = 0; i < GetNumVertices(); i++) {
    out << "v " << i << ' ' << GetVertexLabel(i) << ' ' << GetDegree(i) << '\n';
  }
  int idx = 0;
  for (auto &e : edge_list) {
    if (e.first < e.second) {
      out << "e " << e.first << ' ' << e.second << ' ' << GetEdgeLabel(idx)
          << '\n';
    }
    idx++;
  }
}

} // namespace GraphLib