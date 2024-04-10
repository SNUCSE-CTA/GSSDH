#pragma once
#include "DataStructure/LabeledGraph.h"

static int NumG1Vertices = 0, NumG2Vertices = 0;
struct State {
    State* parent = NULL;
    int cost = 0, depth = -1, lower_bound = 0;
    int vertex_label_bound = 0, inner_edge_label_bound = 0, cross_edge_label_bound = 0;
    int* mapping;
    int* inverse_mapping;
    State(State* p, int c = -1, int d = -1) {
        parent = p;
        cost = c; depth = d;
        lower_bound = 0;
        vertex_label_bound = inner_edge_label_bound = cross_edge_label_bound = 0;
        mapping = new int[NumG1Vertices];
        inverse_mapping = new int[NumG2Vertices];
        if (p != NULL) {
            std::memcpy(mapping, p->mapping, sizeof(int) * NumG1Vertices);
            std::memcpy(inverse_mapping, p->inverse_mapping, sizeof(int) * NumG2Vertices);
            depth = p->depth + 1;
            cost = p->cost;
        }
        else {
            memset(mapping, -1, sizeof(int) * NumG1Vertices);
            memset(inverse_mapping, -1, sizeof(int) * NumG2Vertices);
        }
    };
    void ComputeLowerBound() {
        lower_bound = cost + vertex_label_bound + inner_edge_label_bound + cross_edge_label_bound;
    }
};

struct StateComparator {
    bool operator()(const State* a, const State* b) const {
        if (a->lower_bound == b->lower_bound) {
            return a->depth < b->depth;
        }
        return a->lower_bound > b->lower_bound;
    }
};