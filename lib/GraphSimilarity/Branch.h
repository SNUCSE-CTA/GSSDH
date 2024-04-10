#pragma once
#include "Base/BasicAlgorithms.h"
#include <vector>

namespace GraphLib::GraphSimilarity {
    struct Branch {
        int vertex_label;
        std::vector<int> edge_labels;

        int GetVertexLabel() const { return vertex_label; }
        const std::vector<int>& GetEdgeLabels() const { return edge_labels; }

        const bool operator<(const Branch& other) const {
            if (vertex_label != other.vertex_label)
                return vertex_label  < other.vertex_label;
            return edge_labels < other.edge_labels;
        }
    };

    int BranchEditDistance(const Branch& a, const Branch& b) {
        int ed = 0;
        if (a.vertex_label!= b.vertex_label) ed += 2;
        ed += std::max(a.edge_labels.size(), b.edge_labels.size());
        ed -= VectorIntersectionSize(a.edge_labels, b.edge_labels);
        return ed;
    }

    int BranchEditDistanceFromNull(const Branch& a) {
        return 2 + a.edge_labels.size();
    }
}