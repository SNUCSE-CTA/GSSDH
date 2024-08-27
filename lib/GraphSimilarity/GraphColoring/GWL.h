#pragma once
#include <stack>
#include <unordered_set>

#include "DataStructure/LabeledGraph.h"

using GraphLib::LabeledGraph;
namespace GraphLib {
class GWL {
 protected:
  LabeledGraph *G;
  std::unordered_map<int, std::vector<int>> color_tree_map;
  std::unordered_map<int, int> child_to_parent_map;

 public:
  void GraphColoring(const int iter,
                     const std::unordered_set<int> &fixed_vertices);
};

struct MapEqual {
  bool operator()(const std::map<int, int> &m1,
                  const std::map<int, int> &m2) const {
    return m1 == m2;
  }
};

void GWL::GraphColoring(const int iter,
                        const std::unordered_set<int> &fixed_vertices) {
  std::unordered_map<std::map<int, int>, int, MapEqual> frequency_to_color_map;
  G->vertex_color.resize(G->GetNumVertices());
  G->vertex_color.assign(G->GetNumVertices(), 0);
  int next_color = 1;
  color_tree_map[0];

  for (int i = 0; i < iter; ++i) {
    for (int u = 0; u < G->GetNumVertices(); ++u) {
      if (fixed_vertices.find(u) != fixed_vertices.end()) {
        continue;
      }

      std::vector<int> neighbors = G->GetNeighbors(u);

      std::map<int, int> color_frequency;
      for (int v : neighbors) {
        color_frequency[G->vertex_color[v]]++;
      }

      auto it = frequency_to_color_map.find(color_frequency);
      int old_color = G->vertex_color[u];

      if (it != frequency_to_color_map.end()) {
        G->vertex_color[u] = it->second;
      } else {
        frequency_to_color_map[color_frequency] = next_color;
        G->vertex_color[u] = next_color;
        color_tree_map[old_color].push_back(next_color);
        color_tree_map[next_color];
        child_to_parent_map[next_color] = old_color;
        next_color++;
      }
    }
  }
}
}  // namespace GraphLib
