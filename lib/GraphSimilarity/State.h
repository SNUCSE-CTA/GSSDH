#pragma once
#include "DataStructure/LabeledGraph.h"

static int NumG1Vertices = 0, NumG2Vertices = 0;
struct State {
  State* parent = NULL;
  int id = 0, next_mapping_order = -1;
  int cost = 0, depth = -1, lower_bound = 0;
  int vertex_label_bound = 0, inner_edge_label_bound = 0,
      cross_edge_label_bound = 0;
  int* mapping;
  int* inverse_mapping;
  static int global_state_id;
  State(State* p, int c = -1, int d = -1) {
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
  State(const State& other) {
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

  // Destructor
  ~State() {
    delete[] mapping;
    delete[] inverse_mapping;
  }

  // Assignment operator
  State& operator=(const State& other) {
    if (this == &other) return *this;

    // Free existing memory
    delete[] mapping;
    delete[] inverse_mapping;

    // Copy data
    id = other.id;
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

    return *this;
  }
  void ComputeLowerBound() {
    lower_bound = cost + vertex_label_bound + inner_edge_label_bound +
                  cross_edge_label_bound;
  }

  std::string to_string() const {
    std::ostringstream oss;
    oss << "State " << id << ": ";
    for (int i = 0; i < NumG1Vertices; i++) {
      oss << mapping[i] << " ";
    }
    oss << "| Cost = " << cost << " | Lb = " << lower_bound << "";
    return oss.str();
  }

  void Print() {
    fprintf(stderr, "State[%d]: ", id);
    for (int i = 0; i < NumG1Vertices; i++) {
      fprintf(stderr, "%d ", mapping[i]);
    }
    fprintf(stderr, "| Cost = %d | Lb = %d\n\n", cost, lower_bound);
  }
};

struct StateComparator {
  bool operator()(const State* a, const State* b) const {
    if (a->lower_bound == b->lower_bound) {
      if (a->depth == b->depth) {
        return a->id < b->id;
      }
      return a->depth < b->depth;
    }
    return a->lower_bound > b->lower_bound;
  }
};

int State::global_state_id = 0;
