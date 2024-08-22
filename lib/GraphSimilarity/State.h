#pragma once
#include "DataStructure/LabeledGraph.h"

static int NumG1Vertices = 0, NumG2Vertices = 0;
struct State {
  State *parent = NULL;
  int id = 0, next_mapping_order = -1;
  int cost = 0, depth = -1, lower_bound = 0;
  int vertex_label_bound = 0, inner_edge_label_bound = 0,
      cross_edge_label_bound = 0;
  int *mapping;
  int *inverse_mapping;
  static int global_state_id;
  State(State *p, int c = -1, int d = -1) {
    id = global_state_id++;
    parent = p;
    cost = c;
    depth = d;
    lower_bound = 0;
    vertex_label_bound = inner_edge_label_bound = cross_edge_label_bound = 0;
    mapping = new int[NumG1Vertices];
    inverse_mapping = new int[NumG2Vertices];
    if (p != NULL) {
      std::memcpy(mapping, p->mapping, sizeof(int) * NumG1Vertices);
      std::memcpy(inverse_mapping, p->inverse_mapping,
                  sizeof(int) * NumG2Vertices);
      depth = p->depth + 1;
      cost = p->cost;
    } else {
      memset(mapping, -1, sizeof(int) * NumG1Vertices);
      memset(inverse_mapping, -1, sizeof(int) * NumG2Vertices);
    }
  };

  // Copy constructor
  State(const State &other) {
    id = global_state_id++;
    parent = other.parent;
    next_mapping_order = other.next_mapping_order;
    cost = other.cost;
    depth = other.depth;
    lower_bound = other.lower_bound;
    vertex_label_bound = other.vertex_label_bound;
    inner_edge_label_bound = other.inner_edge_label_bound;
    cross_edge_label_bound = other.cross_edge_label_bound;
    mapping = new int[NumG1Vertices];
    inverse_mapping = new int[NumG2Vertices];
    std::memcpy(mapping, other.mapping, sizeof(int) * NumG1Vertices);
    std::memcpy(inverse_mapping, other.inverse_mapping,
                sizeof(int) * NumG2Vertices);
  }
};

int State::global_state_id = 0;
