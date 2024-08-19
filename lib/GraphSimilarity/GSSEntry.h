#pragma once
#include "Branch.h"
#include "DataStructure/LabeledGraph.h"

namespace GraphLib::GraphSimilarity {
    /**
 * Since the DB is a collection of small graphs,
 * a seperate class inheriting LabeledGraph is needed.
 * @note Unlike Labeledgraph, this class does not store edges in bi-directional way.
 */
    class GSSEntry : public LabeledGraph {
        std::vector<std::vector<int>> adj_matrix;
        std::vector<std::vector<int>> incidence_list;
        std::vector<std::vector<int>> num_incident_edge_by_label;
        std::vector<std::vector<int>> vertices_by_label;
        std::vector<std::vector<int>> edges_by_label;
        std::vector<int> degree_sequence;
        std::vector<int> branch_ids;
        std::vector<Branch> branches;

    public:
        // Works only for small graphs
        inline int GetEdgeLabel(int edge_id) const {return edge_label[edge_id];}
        inline int GetEdgeLabel(int u, int v) const {return adj_matrix[u][v];}
        inline int NumIncidentEdgesWithLabel(int v, int label) {
            if (label >= num_edge_labels) return 0;
            return num_incident_edge_by_label[v][label];
        }
        inline std::vector<int>& GetIncidentEdges(int v) {return incidence_list[v];}
        inline std::vector<int>& GetVerticesByLabel(int label) {return vertices_by_label[label];}
        inline std::vector<int>& GetEdgesByLabel(int label) {return edges_by_label[label];}
        inline std::vector<int>& GetDegreeSequence() {return degree_sequence;}
        const Branch& GetBranch(int v) {return branches[v];}
        void BuildDegreeSequence() {
            degree_sequence.resize(num_vertex);
            for (int i = 0; i < num_vertex; i++) {
                degree_sequence[i] = adj_list[i].size();
            }
            std::sort(degree_sequence.begin(), degree_sequence.end(), std::greater<int>());
        }

        void LoadFromDatabase(int graph_id, std::vector<int> &vertices, std::vector<int> &vertex_labels,
                              std::vector<std::pair<int, int>> &edges, std::vector<int> &edge_labels) {
            id = graph_id;
            num_vertex = vertices.size();
            num_vertex_labels = 0;
            num_edge = edges.size();
            num_edge_labels = 0;
            vertex_label.resize(num_vertex);
            adj_list.resize(num_vertex);
            adj_matrix.resize(num_vertex);
            incidence_list.resize(num_vertex);
            num_incident_edge_by_label.resize(num_vertex);
            for (int i = 0; i < num_vertex; i++) {
                num_vertex_labels = std::max(num_vertex_labels, vertex_labels[i] + 1);
                vertex_label[vertices[i]] = vertex_labels[i];
                vertex_label_frequency[vertex_labels[i]]++;
                adj_matrix[i].resize(num_vertex, -1);
            }
            vertices_by_label.resize(num_vertex_labels);
            edge_list = edges;
            edge_label = edge_labels;
            for (int i = 0; i < edge_list.size(); i++) {
                auto &[u, v] = edge_list[i];
                num_edge_labels = std::max(num_edge_labels, edge_label[i] + 1);
                adj_matrix[u][v] = edge_label[i];
                adj_matrix[v][u] = edge_label[i];
                adj_list[u].push_back(v);
                adj_list[v].push_back(u);
                incidence_list[u].push_back(i);
                incidence_list[v].push_back(i);
                edge_label_frequency[edge_label[i]]++;
            }
            edges_by_label.resize(num_edge_labels);
            for (int i = 0; i < num_vertex; i++) {
                vertices_by_label[GetVertexLabel(i)].push_back(i);
                num_incident_edge_by_label[i].resize(num_edge_labels, 0);
                for (int e : incidence_list[i]) {
                    num_incident_edge_by_label[i][GetEdgeLabel(e)]++;
                }
            }
            for (int i = 0; i < edge_list.size(); i++) {
                edges_by_label[GetEdgeLabel(i)].push_back(i);
            }
            BuildDegreeSequence();
        };

        void CombineGraph(GSSEntry *g1, GSSEntry *g2) {
            num_vertex = g1->GetNumVertices() + g2->GetNumVertices();
            adj_matrix.resize(num_vertex, std::vector<int>(num_vertex, -1));
            num_edge = g1->GetNumEdges() + g2->GetNumEdges();
            adj_list.resize(num_vertex);
            vertex_label.resize(num_vertex);
            edge_label.resize(num_edge);
            for (int i = 0; i < g1->GetNumVertices(); i++) vertex_label[i] = g1->GetVertexLabel(i);
            for (int i = 0; i < g2->GetNumVertices(); i++) vertex_label[i + g1->GetNumVertices()] = g2->GetVertexLabel(i);
            for (int i = 0; i < g1->GetNumEdges(); i++) {
                edge_label[i] = g1->GetEdgeLabel(i);
                auto [v1, v2] = g1->GetEdge(i);
                adj_list[v1].push_back(v2);
                adj_list[v2].push_back(v1);
                adj_matrix[v1][v2] = adj_matrix[v2][v1] = edge_label[i];
            }
            for (int i = 0; i < g2->GetNumEdges(); i++) {
                edge_label[i + g1->GetNumEdges()] = g2->GetEdgeLabel(i);
                auto [v1, v2] = g2->GetEdge(i);
                v1 += g1->GetNumVertices();
                v2 += g1->GetNumVertices();
                adj_list[v1].push_back(v2);
                adj_list[v2].push_back(v1);
                adj_matrix[v1][v2] = adj_matrix[v2][v1] = edge_label[i + g1->GetNumEdges()];
            }
        }

        void BuildBranches(std::map<Branch, int>& seen_branch_ids, std::vector<Branch>& seen_branch_structures) {
            branches.resize(num_vertex);
            branch_ids.resize(num_vertex);
            for (int i = 0; i < num_vertex; i++) {
                std::sort(adj_list[i].begin(), adj_list[i].end(), [this, i](int u, int v) {
                    return adj_matrix[i][u] < adj_matrix[i][v];
                });
                branches[i].vertex_label = vertex_label[i];
                branches[i].edge_labels.resize(adj_list[i].size());
                for (int j = 0; j < adj_list[i].size(); j++) {
                    branches[i].edge_labels[j] = adj_matrix[i][adj_list[i][j]];
                }
                auto it = seen_branch_ids.find(branches[i]);
                if (it == seen_branch_ids.end()) {
                    branch_ids[i] = seen_branch_structures.size();
                    seen_branch_ids[branches[i]] = seen_branch_structures.size();
                    seen_branch_structures.push_back(branches[i]);
                }
                else {
                    branch_ids[i] = it->second;
                }
            }
        }

        const int GetBranchID(int v) const {return branch_ids[v];}
        std::vector<int>& GetBranchIDs() {return branch_ids;}
    };

}