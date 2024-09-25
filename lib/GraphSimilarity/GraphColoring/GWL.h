#pragma once
#include <stack>
#include <queue>
#include <unordered_set>

#include "GraphSimilarity/GSSEntry.h"
#include "Base/Hungarian.h"

using GraphLib::LabeledGraph;
namespace GraphLib {
class ColorTree {
public:
  ColorTree *parent;
  // int color;
  std::vector<int> vertices;
  std::vector<ColorTree *> children;
  int height;
  bool was_in_queue = false;
  ColorTree() : parent(nullptr) { was_in_queue = false; }
};

class GWL {
protected:
  // The input graph (combined graph)
  GraphSimilarity::GSSEntry *G;
  // The color tree
  ColorTree *root;
  // The leaf nodes
  std::queue<ColorTree *> leaf_nodes;
  // The mapping from vertex to leaf node
  std::vector<ColorTree *> vertex_to_leaf_node;
  const int INF = 1e9;

public:
  // The matching information
  std::vector<int> mapping;
  std::vector<int> inverse_mapping;

  GWL(GraphSimilarity::GSSEntry *g) : G(g) {}

  void SetGraph(GraphSimilarity::GSSEntry *g) { G = g; }

  void GraphColoring(const int iter, int *mapping_array,
                     int *inverse_mapping_array, ColorTree **prev_color_to_node,
                     ColorTree **curr_color_to_node);

  void VertexMatching();

  void ComputeMatchingCost();

  void MakeBipartiteGraph();

  void Init() {
    root = new ColorTree();
    root->height = 0;
    vertex_to_leaf_node.resize(G->GetNumVertices());
  }

  void Deallocate() {
    std::stack<ColorTree *> stack;
    stack.push(root);
    while (!stack.empty()) {
      ColorTree *node = stack.top();
      stack.pop();
      for (ColorTree *child : node->children) {
        stack.push(child);
      }
      delete node;
    }
  }

  void debug() {
    std::cout << "Vertex coloring: ";
    for (int i = 0; i < G->GetNumVertices(); ++i) {
      std::cout << G->vertex_color[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Color tree: " << std::endl;
    std::stack<ColorTree *> stack;
    stack.push(root);
    while (!stack.empty()) {
      ColorTree *node = stack.top();
      stack.pop();
      std::cout << "Vertices: ";
      for (int v : node->vertices) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
    }
    int i = 0, num_leaf_nodes = leaf_nodes.size();
    std::cout << "Leaf nodes: ";
    while (i < num_leaf_nodes) {
      ColorTree *node = leaf_nodes.front();
      leaf_nodes.pop();
      std::cout << "Parent: ";
      for (int v : node->parent->vertices) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
      std::cout << "Vertices: ";
      for (int v : node->vertices) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
      leaf_nodes.push(node);
      i++;
    }
    // Print matching
    std::cout << "Matching: " << std::endl;
    for (int i = 0; i < G->combined_index; ++i) {
      std::cout << i << " -> " << mapping[i] + G->combined_index << std::endl;
    }
    for (int i = 0; i < G->GetNumVertices() - G->combined_index; ++i) {
      std::cout << i + G->combined_index << " -> " << inverse_mapping[i]
                << std::endl;
    }
  }
};

void GWL::GraphColoring(const int iter, int *mapping_array,
                        int *inverse_mapping_array,
                        ColorTree **prev_color_to_node,
                        ColorTree **curr_color_to_node) {
  // The hash map to store the frequency of the color
  std::unordered_map<long long int, int> frequency_to_color_map;

  // Initialize the color of each vertex
  int next_color = 0;
  int num_vertices = G->GetNumVertices();
  int relaxed_color[1000];
  for (int i = 0; i < 1000; ++i) {
    relaxed_color[i] = -1;
  }
  int combined_index = G->combined_index;
  for (int i = 0; i < num_vertices; ++i) {
    if ((i < combined_index && mapping_array[i] == -1) ||
        (i >= combined_index &&
         inverse_mapping_array[i - combined_index] == -1)) {
      int curr_color = G->GetVertexLabel(i);
      if (relaxed_color[curr_color] == -1) {
        relaxed_color[curr_color] = next_color;
        G->vertex_color[i] = next_color++;
      } else {
        G->vertex_color[i] = relaxed_color[curr_color];
      }
      curr_color = G->vertex_color[i];
      if (prev_color_to_node[curr_color] == nullptr) {
        ColorTree *node = new ColorTree();
        node->parent = root;
        node->height = root->height + 1;
        root->children.push_back(node);
        prev_color_to_node[curr_color] = node;
      }
    } else {
      if (i < combined_index) {
        G->vertex_color[i] = next_color;
        ColorTree *node = new ColorTree();
        node->parent = root;
        root->children.push_back(node);
        node->vertices.push_back(i);
        node->height = root->height + 1;
        prev_color_to_node[next_color] = node;
        leaf_nodes.push(node);
        vertex_to_leaf_node[i] = node;
        next_color++;
      } else {
        G->vertex_color[i] =
            G->vertex_color[inverse_mapping_array[i - combined_index]];
        prev_color_to_node[G->vertex_color[i]]->vertices.push_back(i);
        vertex_to_leaf_node[i] = prev_color_to_node[G->vertex_color[i]];
      }
    }
  }

  // Temporary color for each vertex
  std::vector<int> temp_color;
  temp_color.resize(num_vertices);

  // The number of edge labels
  const int num_edge_labels = G->GetNumEdgeLabels();

  // Iterate the coloring process
  for (int i = 0; i < iter; ++i) {
    const int num_vertex_colors = next_color;
    const int digit_value = (num_vertex_colors - 1) * num_edge_labels;
    next_color = 0;
    for (int u = 0; u < num_vertices; ++u) {
      // Skip the matched vertices
      if ((u < combined_index && mapping_array[u] != -1) ||
          (u >= combined_index &&
           inverse_mapping_array[u - combined_index] != -1)) {
        if (u < combined_index) {
          temp_color[u] = next_color++;
        } else {
          temp_color[u] = temp_color[inverse_mapping_array[u - combined_index]];
        }
        continue;
      }
      // Get the color of current vertex
      const int curr_color = G->vertex_color[u];

      // Get the neighbors of the vertex
      std::vector<int> neighbors = G->GetNeighbors(u);
      // Get the color frequency of the neighbor (l(v), l(e)) pairs
      std::vector<int> hash_values;
      for (int v : neighbors) {
        const int color = G->vertex_color[v];
        const int edge_label = G->GetEdgeLabel(u, v);
        int hash_value = edge_label * num_vertex_colors + color;
        hash_values.push_back(hash_value);
      }
      hash_values.push_back(curr_color);

      // Sort the hash values
      std::sort(hash_values.begin(), hash_values.end());
      long long int color_frequency = 0;
      for (int hash_value : hash_values) {
        color_frequency = color_frequency * digit_value + hash_value;
      }

      // Find the color frequency in the map
      auto freq_iter = frequency_to_color_map.find(color_frequency);

      // If the color frequency is already in the map, assign the color to the
      // vertex. Otherwise, assign a new color to the vertex.
      if (freq_iter != frequency_to_color_map.end()) {
        temp_color[u] = freq_iter->second;
        curr_color_to_node[freq_iter->second]->vertices.push_back(u);
        if (i == iter - 1) {
          vertex_to_leaf_node[u] = curr_color_to_node[freq_iter->second];
        }
      } else {
        frequency_to_color_map[color_frequency] = next_color;
        temp_color[u] = next_color;
        ColorTree *node = new ColorTree();
        node->parent = prev_color_to_node[curr_color];
        node->parent->vertices.clear();
        prev_color_to_node[curr_color]->children.push_back(node);
        curr_color_to_node[next_color] = node;
        node->vertices.push_back(u);
        node->height = node->parent->height + 1;
        next_color++;
        if (i == iter - 1) {
          leaf_nodes.push(node);
          vertex_to_leaf_node[u] = node;
        }
      }
    }
    for (int u = 0; u < num_vertices; ++u) {
      G->vertex_color[u] = temp_color[u];
      prev_color_to_node[u] = curr_color_to_node[u];
    }
    frequency_to_color_map.clear();
  }
}

void GWL::VertexMatching() {
  root->was_in_queue = false;
  const int combined_index = G->combined_index;
  mapping = std::vector<int>(combined_index, -1);
  inverse_mapping = std::vector<int>(G->GetNumVertices() - combined_index, -1);
  while (!leaf_nodes.empty()) {
    ColorTree *leaf_node = leaf_nodes.front();
    leaf_nodes.pop();

    // Split the leaf node into two groups by the combined index
    const int num_vertices = leaf_node->vertices.size();
    int G1_index = 0, G2_index = 0, num_G1_vertices = 0, num_G2_vertices = 0;
    std::sort(leaf_node->vertices.begin(), leaf_node->vertices.end());
    for (int i = 0; i < num_vertices; ++i) {
      if (leaf_node->vertices[i] < combined_index) {
        num_G1_vertices++;
      } else {
        break;
      }
    }
    G2_index = num_G1_vertices;
    num_G2_vertices = num_vertices - num_G1_vertices;

    // Match the vertices in the two groups
    while (G1_index < num_G1_vertices && G2_index < num_vertices) {
      const int u = leaf_node->vertices[G1_index];
      const int v = leaf_node->vertices[G2_index];
      mapping[u] = v - combined_index;
      inverse_mapping[v - combined_index] = u;
      G1_index++;
      G2_index++;
    }

    // Add the remaining vertices to the parent node
    if (num_G1_vertices != num_G2_vertices && leaf_node != root) {
      int start_index = num_G1_vertices > num_G2_vertices ? G1_index : G2_index;
      int end_index =
          num_G1_vertices > num_G2_vertices ? num_G1_vertices : num_vertices;
      ColorTree *parent = leaf_node->parent;
      for (int i = start_index; i < end_index; ++i) {
        parent->vertices.push_back(leaf_node->vertices[i]);
      }
      if (parent->was_in_queue == false) {
        parent->was_in_queue = true;
        leaf_nodes.push(parent);
      }
    }
  }
}

void GWL::MakeBipartiteGraph() {
  // Debug coloring
  std::cout << "Vertex coloring: ";
  for (int i = 0; i < G->GetNumVertices(); ++i) {
    std::cout << i + 1 << " -> " << G->vertex_color[i] << std::endl;
  }
  std::cout << std::endl;
  // Make the bipartite graph based on the distance on the color tree
  std::vector<std::vector<int>> cost_matrix;
  int combined_index = G->combined_index;
  int matrix_size =
      std::max(combined_index, G->GetNumVertices() - combined_index);
  cost_matrix.resize(matrix_size);
  for (int i = 0; i < matrix_size; ++i) {
    cost_matrix[i].resize(matrix_size);
  }
  for (int i = 0; i < matrix_size; ++i) {
    for (int j = 0; j < matrix_size; ++j) {
      cost_matrix[i][j] = INF;
    }
  }
  for (int i = 0; i < G->GetNumVertices(); ++i) {
    std::cout << "Vertex " << i << ", leaf node: " << vertex_to_leaf_node[i]
              << std::endl;
  }
  for (int i = 0; i < combined_index; ++i) {
    for (int j = combined_index; j < G->GetNumVertices(); ++j) {
      ColorTree *node = vertex_to_leaf_node[i];
      ColorTree *leaf_node = vertex_to_leaf_node[j];
      // The initial distance is the diff of neighbor color multiset
      std::vector<int> node_neighbor_colors;
      std::vector<int> leaf_node_neighbor_colors;
      for (int u : G->GetNeighbors(i)) {
        node_neighbor_colors.push_back(G->vertex_color[u]);
      }
      for (int u : G->GetNeighbors(j)) {
        leaf_node_neighbor_colors.push_back(G->vertex_color[u]);
      }
      int distance =
          node_neighbor_colors.size() - leaf_node_neighbor_colors.size();
      if (distance < 0) {
        distance = -distance;
      }
      std::sort(node_neighbor_colors.begin(), node_neighbor_colors.end());
      std::sort(leaf_node_neighbor_colors.begin(),
                leaf_node_neighbor_colors.end());
      std::vector<int> diff;
      std::set_difference(
          node_neighbor_colors.begin(), node_neighbor_colors.end(),
          leaf_node_neighbor_colors.begin(), leaf_node_neighbor_colors.end(),
          std::inserter(diff, diff.begin()));
      distance += diff.size();
      // Debug
      if (i == 4 && j == 8) {
        std::cout << "node_vector: ";
        for (int u : node_neighbor_colors) {
          std::cout << u << " ";
        }
        std::cout << std::endl;
        std::cout << "leaf_node_vector: ";
        for (int u : leaf_node_neighbor_colors) {
          std::cout << u << " ";
        }
        std::cout << std::endl;
        std::cout << "Diff: ";
        for (int u : diff) {
          std::cout << u << " ";
        }
        std::cout << std::endl;
        std::cout << "Distance: " << distance << std::endl;
      }
      while (node != leaf_node) {
        if (node->height > leaf_node->height) {
          node = node->parent;
          distance += 2;
        } else {
          leaf_node = leaf_node->parent;
          distance += 2;
        }
      }
      cost_matrix[i][j - combined_index] = distance;
    }
  }
  // Debug
  for (int i = 0; i < matrix_size; ++i) {
    for (int j = 0; j < matrix_size; ++j) {
      if (cost_matrix[i][j] == INF) {
        std::cout << "INF ";
      } else {
        std::cout << cost_matrix[i][j] << " ";
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Hungarian algorithm
  Hungarian hungarian(cost_matrix);
  hungarian.Solve();
  std::vector<int> assignment = hungarian.GetAssignment();
  mapping = std::vector<int>(combined_index, -1);
  inverse_mapping = std::vector<int>(G->GetNumVertices() - combined_index, -1);
  for (int i = 0; i < combined_index; ++i) {
    if (assignment[i] != -1) {
      mapping[i] = assignment[i];
      inverse_mapping[assignment[i]] = i;
    }
  }

  // Debug
  for (int i = 0; i < combined_index; ++i) {
    std::cout << i + 1 << " -> " << mapping[i] + 1 << std::endl;
  }
  std::cout << std::endl;
}
} // namespace GraphLib
