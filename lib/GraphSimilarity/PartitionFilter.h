#pragma once
#include <unordered_set>
#include "GraphSimilarity/GSSEntry.h"

namespace GraphLib::GraphSimilarity {
    struct MismatchingStructure {
        std::vector<int> vertices, edges;
        MismatchingStructure(std::vector<int> vertices_, std::vector<int> edges_) : vertices(vertices_), edges(edges_) {
//            fprintf(stderr, "Mismatch found)");
//            fprintf(stderr, " Vertices: ");
//            for (auto v : vertices) {
//                fprintf(stderr, "%d ", v);
//            }
//            fprintf(stderr, "| Edges: ");
//            for (auto e : edges) {
//                fprintf(stderr, "%d ", e);
//            }
//            fprintf(stderr, "\n");
        };
    };

    class PartitionFilter {
        GSSEntry *query = nullptr, *data = nullptr;
        std::vector<MismatchingStructure> mismatches;
        std::vector<int> used_vertex, used_edge;
    public:
        PartitionFilter(GSSEntry *query_ = nullptr, GSSEntry *data_ = nullptr) : query(query_), data(data_) {
            BuildPartition();
        };

        int GetPartitionBound() {
            return mismatches.size();
        }

        void BuildPartition() {
//            fprintf(stderr, "Building partition for %d %d\n",query->GetId(), data->GetId());
            used_vertex.resize(query->GetNumVertices(), 0);
            used_edge.resize(query->GetNumEdges(), 0);

            // Size 1 structures of query that are not in data
            for (int i = 0; i < query->GetNumVertices(); i++) {
                // If label is not in data
                int vl = query->GetVertexLabel(i);
                if (!data->HaveVertexLabel(vl)) {
                    mismatches.emplace_back(std::vector<int>{i}, std::vector<int>());
                    used_vertex[i] = 1;
                }
            }
            for (int i = 0; i < query->GetNumEdges(); i++) {
                int el = query->GetEdgeLabel(i);
                if (!data->HaveEdgeLabel(el)) {
                    mismatches.emplace_back(std::vector<int>(), std::vector<int>{i});
                    used_edge[i] = 1;
                }
            }

            // Size 2 structures of query that are not in data
            std::unordered_set<std::pair<int, int>> size_2_label_freq;
            for (int i = 0; i < data->GetNumVertices(); i++) {
                for (int e : data->GetIncidentEdges(i)) {
                    size_2_label_freq.emplace(data->GetVertexLabel(i), data->GetEdgeLabel(e));
                }
            }
            for (int i = 0; i < query->GetNumVertices(); i++) {
                if (used_vertex[i]) continue;
                for (int e : query->GetIncidentEdges(i)) {
                    if (used_edge[i]) continue;
                    if (size_2_label_freq.find({query->GetVertexLabel(i), query->GetEdgeLabel(e)}) == size_2_label_freq.end()) {
                        mismatches.emplace_back(std::vector<int>{i}, std::vector<int>{e});
                        used_vertex[i] = 1;
                        used_edge[e] = 1;
                        break;
                    }
                }
            }

            // Size 3 structures: V-E-V
            for (int i = 0; i < query->GetNumEdges(); i++) {
                if (used_edge[i]) continue;
                int u, v; std::tie(u, v) = query->GetEdge(i);
                if (used_vertex[u] or used_vertex[v]) continue;
                int ul = query->GetVertexLabel(u), vl = query->GetVertexLabel(v);
                if (ul > vl) std::swap(ul, vl);
                int el = query->GetEdgeLabel(i);
                for (int de : data->GetEdgesByLabel(el)) {
                    int du, dv; std::tie(du, dv) = data->GetEdge(de);
                    int dul = data->GetVertexLabel(du), dvl = data->GetVertexLabel(dv);
                    if (dul > dvl) std::swap(dul, dvl);
                    if (ul == dul && vl == dvl) {
                        goto nxt_edge;
                    }
                }
                mismatches.emplace_back(std::vector<int>{u, v}, std::vector<int>{i});
                used_vertex[u] = used_vertex[v] = 1;
                used_edge[i] = 1;
                nxt_edge:continue;
            }

            // Size 3 structures: E-V-E
            for (int v = 0; v < query->GetNumVertices(); v++) {
                if (used_vertex[v]) continue;
                int vl = query->GetVertexLabel(v);
                auto &inc = query->GetIncidentEdges(v);
                for (int i = 0; i < inc.size(); i++) {
                    if (used_edge[inc[i]]) continue;
                    for (int j = i+1; j < inc.size(); j++) {
                        if (used_edge[inc[j]]) continue;
                        if (query->GetEdgeLabel(inc[i]) == query->GetEdgeLabel(inc[j])) {
                            int el = query->GetEdgeLabel(inc[i]);
                            bool found = false;
                            for (int dv : data->GetVerticesByLabel(vl)) {
                                if (data->NumIncidentEdgesWithLabel(dv, el) >= 2) {
                                    found = true; break;
                                }
                            }
                            if (!found) {
                                mismatches.emplace_back(std::vector<int>{v}, std::vector<int>{inc[i], inc[j]});
                                used_vertex[v] = 1;
                                used_edge[inc[i]] = used_edge[inc[j]] = 1;
                                goto nxt_vertex;
                            }
                        }
                    }
                }
                nxt_vertex:continue;
            }


            return;
        }
    };
}