#pragma once
#include <stack>
#include <unordered_set>

#include "GraphSimilarity/GSSEntry.h"

using GraphLib::LabeledGraph;
namespace GraphLib {
class GWL {
protected:
  // The input graph (combined graph)
  GraphSimilarity::GSSEntry *G;
  // The color tree
  std::unordered_map<int, std::vector<int>> color_tree_map;
  // The structure to store the vertices of each color
  std::unordered_map<int, std::vector<int>> color_node_to_vertex;
  // The structure to store the parent of each color
  std::unordered_map<int, int> child_to_parent_map;

public:
  GWL(GraphSimilarity::GSSEntry *g) : G(g) {}
  void GraphColoring(const int iter,
                     const std::unordered_set<int> &fixed_vertices);
  void debug() {
    std::cout << "Vertex color: ";
    for (int i = 0; i < G->GetNumVertices(); ++i) {
      std::cout << G->vertex_color[i] << " ";
    }
    std::cout << "Child to parent map: " << std::endl;
    for (auto &it : child_to_parent_map) {
      std::cout << it.first << ": " << it.second << std::endl;
    }
    std::cout << "Color node to vertex: " << std::endl;
    for (auto &it : color_node_to_vertex) {
      std::cout << it.first << ": ";
      for (int v : it.second) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "Color tree map: " << std::endl;
    for (auto &it : color_tree_map) {
      std::cout << it.first << ": ";
      for (int v : it.second) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
    }
  }
};

struct MapHash {
  std::size_t operator()(const std::map<int, int> &m) const {
    std::size_t seed = 0;
    for (const auto &pair : m) {
      seed ^= std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second);
    }
    return seed;
  }
};

struct MapEqual {
  bool operator()(const std::map<int, int> &m1,
                  const std::map<int, int> &m2) const {
    return m1 == m2;
  }
};

void GWL::GraphColoring(const int iter,
                        const std::unordered_set<int> &fixed_vertices) {
  // The hash map to store the frequency of the color
  std::unordered_map<std::map<int, int>, int, MapHash, MapEqual>
      frequency_to_color_map;
  // The hash map to store the hash value of the (l(v), l(e)) pair
  std::unordered_map<std::pair<int, int>, int> handle_hash_map;
  int handle_hash = 0;

  // Initialize the color of the vertices
  int next_color = G->num_vertex_labels;
  G->vertex_color.resize(G->GetNumVertices());
  for (int i = 0; i < G->GetNumVertices(); ++i) {
    G->vertex_color[i] = G->vertex_label[i];
  }
  color_tree_map[0];
  for (int i = 1; i < next_color; ++i) {
    color_tree_map[0].push_back(i);
  }

  // Temporary color for each vertex
  std::vector<int> temp_color;
  temp_color.resize(G->GetNumVertices());
  for (auto &v : fixed_vertices) {
    temp_color[v] = G->vertex_color[v];
  }

  // Iterate the coloring process
  for (int i = 0; i < iter; ++i) {
    for (int u = 0; u < G->GetNumVertices(); ++u) {
      // Skip the fixed vertices
      if (fixed_vertices.find(u) != fixed_vertices.end()) {
        continue;
      }

      // Get the neighbors of the vertex
      std::vector<int> neighbors = G->GetNeighbors(u);

#ifdef DEBUG
      std::cout << "Vertex: " << u << ", label: " << G->GetVertexLabel(u)
                << std::endl;
      std::cout << "Neighbors: ";
      for (int v : neighbors) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
#endif

      // Get the color frequency of the neighbor (l(v), l(e)) pairs
      std::map<int, int> color_frequency;
      for (int v : neighbors) {
        const int color = G->vertex_color[v];
        auto it = handle_hash_map.find({color, G->GetEdgeLabel(u, v)});
        if (it == handle_hash_map.end()) {
          handle_hash_map[{color, G->GetEdgeLabel(u, v)}] = handle_hash;
          color_frequency[handle_hash]++;
          handle_hash++;
        } else {
          color_frequency[it->second]++;
        }
      }

      // Get the color frequency of the vertex (l(u), -1) pair
      auto handle_iter = handle_hash_map.find({G->vertex_color[u], -1});
      if (handle_iter == handle_hash_map.end()) {
        handle_hash_map[{G->vertex_color[u], -1}] = handle_hash;
        color_frequency[handle_hash]++;
        handle_hash++;
      } else {
        color_frequency[handle_iter->second]++;
      }

      // Find the color frequency in the map
      auto freq_iter = frequency_to_color_map.find(color_frequency);
      int old_color = G->vertex_color[u];

      // If the color frequency is already in the map, assign the color to the
      // vertex. Otherwise, assign a new color to the vertex.
      if (freq_iter != frequency_to_color_map.end()) {
        temp_color[u] = freq_iter->second;
        color_node_to_vertex[freq_iter->second].push_back(u);
      } else {
        frequency_to_color_map[color_frequency] = next_color;
        temp_color[u] = next_color;
        color_tree_map[old_color].push_back(next_color);
        color_tree_map[next_color];
        child_to_parent_map[next_color] = old_color;
        color_node_to_vertex[next_color];
        color_node_to_vertex[next_color].push_back(u);
        next_color++;
      }
    }
    for (int u = 0; u < G->GetNumVertices(); ++u) {
      G->vertex_color[u] = temp_color[u];
    }
  }
}
} // namespace GraphLib
