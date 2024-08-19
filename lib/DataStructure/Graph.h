#pragma once
#include <vector>
/*
 * Basic graph class. All graph class should inherit from this.
 * Assume unlabeled graph with only adjacency list and edge list.
 */
namespace GraphLib {
class Graph {
 protected:
  std::vector<std::vector<int>> adj_list;
  std::vector<std::pair<int, int>> edge_list;
  int num_vertex, num_edge, max_degree, id;

  // Degeneracy and core-number
  std::vector<int> core_num;
  int degeneracy;

 public:
  ~Graph() = default;
  Graph() {};

  // Basic functions to access graph properties
  std::vector<int> &GetNeighbors(int v) { return adj_list[v]; }

  std::vector<std::pair<int, int>> &GetEdges() { return edge_list; }

  inline std::pair<int, int> GetEdge(int i) { return edge_list[i]; }

  inline int GetNumVertices() const { return num_vertex; }

  inline int GetNumEdges() const { return num_edge; }

  inline int GetMaxDegree() const { return max_degree; }

  inline int GetDegree(int v) const { return adj_list[v].size(); }

  inline int GetId() const { return id; }

  // Basic functions to analyze graph data
  void ComputeCoreNum();

  inline int GetCoreNum(int v) const { return core_num[v]; }

  inline int GetDegeneracy() const { return degeneracy; }
};

/**
 * @brief Compute the core number of each vertex
 * @date Oct 21, 2022
 */
void Graph::ComputeCoreNum() {
  core_num.resize(num_vertex, 0);
  int *bin = new int[GetMaxDegree() + 1];
  int *pos = new int[GetNumVertices()];
  int *vert = new int[GetNumVertices()];

  std::fill(bin, bin + (GetMaxDegree() + 1), 0);

  for (int v = 0; v < GetNumVertices(); v++) {
    core_num[v] = adj_list[v].size();
    bin[core_num[v]] += 1;
  }

  int start = 0;
  int num;

  for (int d = 0; d <= GetMaxDegree(); d++) {
    num = bin[d];
    bin[d] = start;
    start += num;
  }

  for (int v = 0; v < GetNumVertices(); v++) {
    pos[v] = bin[core_num[v]];
    vert[pos[v]] = v;
    bin[core_num[v]] += 1;
  }

  for (int d = GetMaxDegree(); d--;) bin[d + 1] = bin[d];
  bin[0] = 0;

  for (int i = 0; i < GetNumVertices(); i++) {
    int v = vert[i];

    for (int u : GetNeighbors(v)) {
      if (core_num[u] > core_num[v]) {
        int du = core_num[u];
        int pu = pos[u];
        int pw = bin[du];
        int w = vert[pw];

        if (u != w) {
          pos[u] = pw;
          pos[w] = pu;
          vert[pu] = w;
          vert[pw] = u;
        }

        bin[du]++;
        core_num[u]--;
      }
    }
  }

  degeneracy = 0;
  for (int i = 0; i < GetNumVertices(); i++) {
    degeneracy = std::max(core_num[i], degeneracy);
  }

  delete[] bin;
  delete[] pos;
  delete[] vert;
}

}  // namespace GraphLib