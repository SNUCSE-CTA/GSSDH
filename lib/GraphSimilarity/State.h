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
};

int State::global_state_id = 0;
