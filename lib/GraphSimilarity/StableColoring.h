#pragma once
#include <stack>

#include "DataStructure/LabeledGraph.h"

using GraphLib::LabeledGraph;
namespace GraphLib {
//@TODO  Implement this as canonical coloring (or better, take as parameter)
/**
 @brief BBG Stable Coloring Algorithm

 Efficient color refinement (also known as 1-WL test, Naive vertex
 classification) by Berkholz, Bonsma and Grohe [1] Time complexity: O((V + E)
 log V), Space complexity: O(V + E)
 @note Current implementation only works for undirected graphs, and does not
 produce canonical coloring (i.e., it is possible that the coloring for two
 isomorphic graphs differ up to permutation) This is to avoid sorting, thereby
 making the algorithm more efficient for pairwise color refinement (used for
 Graph Isomorphism [2]).

       Tested with tree isomorphism problem from Estonian Informatics Olympiad
 2018, Iso testcases.

 References
   1. C. Berkholz, P. Bonsma, and M. Grohe.
        Tight Lower and Upper Bounds for the Complexity of Canonical Colour
 Refinement. Theory of Computing Systems 2017.
   2. G. Gu, Y. Nam, K. Park, Z. Galil, G. F. Italiano, and W.-S. Han,
        Scalable Graph Isomorphism: Combining Pairwise Color Refinement and
 Backtracking via Compressed Candidate Space, ICDE 2021
 */
class BBGColorTree {
public:
  BBGColorTree *parent;
  // int color;
  std::vector<int> vertices;
  std::vector<BBGColorTree *> children;
  int height;
  bool was_in_queue = false;
  BBGColorTree(int h = 0) : parent(nullptr) {
    was_in_queue = false;
    height = h;
  }

  void Debug() { DebugHelper(); }

  void DebugHelper(int indent = 0) {
    std::cout << std::string(indent, ' ') << "Height: " << height << std::endl;
    std::cout << std::string(indent, ' ') << "Vertices: ";
    for (int v : vertices) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << std::string(indent, ' ') << "Children: " << std::endl;
    for (BBGColorTree *child : children) {
      child->DebugHelper(indent + 2);
    }
  }
};

class StableColoring {
protected:
  LabeledGraph *G;
  int *S, stack_size;
  bool *in_stack, *used;
  int *color, *cdeg;
  int *color_cand, *color_split;
  int *max_cdeg, *min_cdeg, *num_cdeg;
  int *color_mapping;
  int *inv_idx;
  int num_color_cand = 0, num_color_split = 0;
  int num_colors = 0;
  int max_color_degree = 0;
  std::vector<std::vector<int>> color_partition, aux;
  BBGColorTree *root;
  std::vector<BBGColorTree *> leaf_nodes;
  std::vector<BBGColorTree *> vertex_to_leaf_node;
  std::vector<BBGColorTree *> color_to_node;
  virtual void SplitUpColor(int s) {
    max_color_degree = max_cdeg[s];
    memset(num_cdeg, 0, sizeof(int) * (max_color_degree + 1));
    num_cdeg[0] = color_partition[s].size() - aux[s].size();
    int b = 0;
    for (int v : aux[s]) {
      num_cdeg[cdeg[v]]++;
      if (num_cdeg[cdeg[v]] > num_cdeg[b])
        b = cdeg[v];
    }
    for (int i = 0; i <= max_color_degree; i++) {
      if (num_cdeg[i] >= 1) {
        int c = color_mapping[i] = (i == min_cdeg[s]) ? s : num_colors++;
        if (in_stack[s]) {
          if (c != s) {
            S[stack_size++] = c;
            in_stack[c] = true;
          }
        } else {
          if (i != b) {
            S[stack_size++] = c;
            in_stack[c] = true;
          }
        }
      }
    }
    BBGColorTree *prev_node = color_to_node[s];
    BBGColorTree *aux_node = new BBGColorTree(prev_node->height + 1);
    // find prev_node in leaf_nodes and remove it
    leaf_nodes.erase(
        std::find(leaf_nodes.begin(), leaf_nodes.end(), prev_node));
    aux_node->parent = prev_node;
    aux_node->vertices = prev_node->vertices;
    prev_node->children.push_back(aux_node);
    prev_node->vertices.clear();
    color_to_node[s] = aux_node;
    leaf_nodes.push_back(aux_node);
    for (int v : aux[s]) {
      int c = color_mapping[cdeg[v]];
      if (c != s) {
        ChangeColor(v, c);
        auto it =
            std::find(aux_node->vertices.begin(), aux_node->vertices.end(), v);
        if (it != aux_node->vertices.end()) {
          aux_node->vertices.erase(it);
        }
        if (color_to_node[c] == nullptr) {
          BBGColorTree *new_node = new BBGColorTree(prev_node->height + 1);
          new_node->parent = prev_node;
          prev_node->children.push_back(new_node);
          color_to_node[c] = new_node;
          vertex_to_leaf_node[v] = new_node;
          leaf_nodes.push_back(new_node);
        } else {
          BBGColorTree *node = color_to_node[c];
          node->vertices.push_back(v);
          vertex_to_leaf_node[v] = node;
        }
      }
    }
  }
  virtual void ChangeColor(int v, int new_color) {
    int s = color[v];
    int v_idx = inv_idx[v];
    inv_idx[color_partition[s].back()] = v_idx;
    color_partition[s][v_idx] = color_partition[s].back();
    color_partition[s].pop_back();
    color_partition[new_color].push_back(v);
    inv_idx[v] = color_partition[new_color].size() - 1;
    color[v] = new_color;
  }

  int *perm;

public:
  StableColoring(LabeledGraph *g) {
    G = g;
    num_color_cand = num_color_split = stack_size = 0;
    color = new int[G->GetNumVertices()]();
    cdeg = new int[G->GetNumVertices()]();
    used = new bool[G->GetNumVertices()]();
    max_cdeg = new int[G->GetNumVertices()]();
    min_cdeg = new int[G->GetNumVertices()]();
    num_cdeg = new int[G->GetNumVertices()]();
    in_stack = new bool[G->GetNumVertices()]();
    color_mapping = new int[G->GetNumVertices()]();
    color_cand = new int[G->GetNumVertices()]();
    color_split = new int[G->GetNumVertices()]();
    inv_idx = new int[G->GetNumVertices()]();
    S = new int[G->GetNumVertices()]();
    color_partition.resize(G->GetNumVertices());
    aux.resize(G->GetNumVertices());
    perm = new int[G->GetNumVertices()]();
    std::iota(perm, perm + G->GetNumVertices(), 0);
    std::sort(perm, perm + G->GetNumVertices(), [&](int i, int j) -> bool {
      int di = G->GetDegree(i), dj = G->GetDegree(j);
      int li = G->GetVertexLabel(i), lj = G->GetVertexLabel(j);
      return (di < dj) or (di == dj && li < lj);
    });
    // give initial coloring according to (label, degree) pair
    int current_color = 0;
    color[perm[0]] = current_color;
    S[stack_size++] = current_color;
    in_stack[current_color] = true;
    inv_idx[perm[0]] = color_partition[current_color].size();
    color_partition[current_color].push_back(perm[0]);
    for (int i = 1; i < G->GetNumVertices(); i++) {
      int u = perm[i], v = perm[i - 1];
      int du = G->GetDegree(u), dv = G->GetDegree(v);
      int lu = G->GetVertexLabel(u), lv = G->GetVertexLabel(v);
      if (du != dv or lu != lv) {
        current_color++;
        S[stack_size++] = current_color;
        in_stack[current_color] = true;
      }
      color[u] = current_color;
      inv_idx[u] = color_partition[current_color].size();
      color_partition[current_color].push_back(u);
      node->vertices.push_back(u);
      vertex_to_leaf_node[u] = node;
    }
    num_colors = current_color + 1;

    // Debug color tree
    // root->Debug();
  }

  void RefineStep() {
    int r = S[--stack_size];
    in_stack[r] = false;
    for (int v : color_partition[r]) {
      for (int w : G->GetNeighbors(v)) {
        if (color_partition[color[w]].size() == 1)
          continue;
        int c = color[w];
        cdeg[w]++;
        if (cdeg[w] == 1)
          aux[c].push_back(w);
        if (!used[c]) {
          color_cand[num_color_cand++] = c;
          used[c] = true;
        }
        if (cdeg[w] > max_cdeg[c])
          max_cdeg[c] = cdeg[w];
      }
    }
    for (int i = 0; i < num_color_cand; i++) {
      int c = color_cand[i];
      if (color_partition[c].size() != aux[c].size())
        min_cdeg[c] = 0;
      else {
        min_cdeg[c] = max_cdeg[c];
        for (int v : aux[c])
          min_cdeg[c] = std::min(min_cdeg[c], cdeg[v]);
      }
      if (min_cdeg[c] < max_cdeg[c])
        color_split[num_color_split++] = c;
    }
    for (int i = 0; i < num_color_split; i++)
      SplitUpColor(color_split[i]);
    num_color_split = 0;
    for (int i = 0; i < num_color_cand; i++) {
      int c = color_cand[i];
      for (int v : aux[c])
        cdeg[v] = 0;
      max_cdeg[c] = 0;
      min_cdeg[c] = 0;
      aux[c].clear();
      used[c] = false;
    }
    num_color_cand = 0;
  }

  virtual void Refine() {
    while (stack_size > 0) {
      RefineStep();
    }
  }
  inline int GetColor(int v) const { return color[v]; }
  inline std::vector<int> &GetVerticesByColor(int c) {
    return color_partition[c];
  }
  void Individualize(int v) {
    if (color_partition[color[v]].size() == 1)
      return;
    int new_color = num_colors++;
    ChangeColor(v, new_color);
    S[stack_size++] = new_color;
    in_stack[new_color] = true;
    Refine();
  }
  inline int GetNumColors() const { return num_colors; }
};

class PairwiseStableColoring : public StableColoring {
private:
  std::vector<std::vector<std::pair<int, int>>> history;
  int VL, VR;
  int *left_color_count, *right_color_count;
  bool color_count_equal = true;

  bool record_history = false;
  int num_colors_before = 0;

  void ChangeColor(int v, int new_color) {
    int s = color[v];
    int v_idx = inv_idx[v];
    (v < VL ? right_color_count : left_color_count)[s]--;
    (v < VL ? right_color_count : left_color_count)[new_color]++;
    inv_idx[color_partition[s].back()] = v_idx;
    color_partition[s][v_idx] = color_partition[s].back();
    color_partition[s].pop_back();
    color_partition[new_color].push_back(v);
    inv_idx[v] = color_partition[new_color].size() - 1;
    color[v] = new_color;
    if (record_history) {
      history.back().emplace_back(v, s);
    }
  }

public:
  std::vector<int> mapping;
  std::vector<int> inverse_mapping;

  PairwiseStableColoring(LabeledGraph *g1, LabeledGraph *g2,
                         LabeledGraph *combined)
      : StableColoring(combined) {
    VL = g1->GetNumVertices();
    VR = g2->GetNumVertices();
    left_color_count = new int[VL + VR]();
    right_color_count = new int[VL + VR]();
    for (int i = 0; i < num_colors; i++) {
      for (auto v : color_partition[i]) {
        (v < VL ? right_color_count : left_color_count)[i]++;
      }
      color_count_equal &= (left_color_count[i] == right_color_count[i]);
    }
  }

  bool CheckColorCount() { return color_count_equal; }

  void Refine() {
    while (stack_size > 0) {
      RefineStep();
      if (!color_count_equal) {
        stack_size = 0;
        return;
      }
    }
  }

  bool PairwiseIndividualize(int u, int v) {
    if ((color[u] == color[v]) and (color_partition[color[u]].size() == 2))
      return true;
    record_history = true;
    history.emplace_back();
    int c = color[u];
    int new_color = num_colors++;
    ChangeColor(u, new_color);
    ChangeColor(v, new_color);
    S[stack_size++] = new_color;
    in_stack[new_color] = true;
    Refine();
#ifdef DEBUG
    printf("Have done %lu work\n", history.back().size());
#endif
    record_history = false;
    return CheckColorCount();
  }

  void SplitUpColor(int s) {
    StableColoring::SplitUpColor(s);
    for (int i = 0; i <= max_color_degree; i++) {
      int c = color_mapping[i];
      color_count_equal &= (left_color_count[c] == right_color_count[c]);
      if (left_color_count[c] == 1 and right_color_count[c] == 1) {
      }
    }
  }

  inline int GetLeftColorCount(int c) const { return left_color_count[c]; }
  inline int GetRightColorCount(int c) const { return right_color_count[c]; }

  inline int GetNumVertexClass(int u) { return left_color_count[color[u]]; }

  void RevertColoring() {
    record_history = false;
#ifdef DEBUG
    printf("Reverting %d bits of history...\n", history.back().size());
#endif
    auto &last_history = history.back();
    for (int i = last_history.size() - 1; i >= 0; i--) {
      std::pair<int, int> c = last_history[i];
      ChangeColor(c.first, c.second);
    }
    history.pop_back();
    while (color_partition[num_colors - 1].size() == 0) {
      num_colors--;
    }
    color_count_equal = true;
    record_history = true;
  }

  inline int BinaryMapping(int u) {
    if (GetNumVertexClass(u) != 1)
      return -1;
    int v = -u;
    for (int x : color_partition[color[u]])
      v += x;
    return v;
  }

  void PrintEntireColorPartition(int upto = 5) {
    fprintf(stdout, "Color info:\n");
    for (int i = 0; i < std::min(GetNumColors(), GetNumColors()); i++) {
      fprintf(stdout, "  Color %d: %d left, %d right\n", i,
              GetLeftColorCount(i), GetRightColorCount(i));
      fprintf(stdout, "    Vertices: ");
      for (int v : GetVerticesByColor(i)) {
        fprintf(stdout, "%d ", v);
      }
      fprintf(stdout, "\n");
    }
    root->Debug();
  }
};
} // namespace GraphLib
