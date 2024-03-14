#pragma once
/*
 * LabeledGraph database for managing collection of small to medium-sized graphs.
 */

#include "LabeledGraph.h"

namespace GraphLib {

    /**
     * Since the DB is a collection of small graphs,
     * a seperate class inheriting LabeledGraph is needed.
     * @note Unlike Labeledgraph, this class does not store edges in bi-directional way.
     */
    class LabeledGraphDatabaseEntry : public LabeledGraph {
        //@TODO Build adjacency matrix only if V <= 500?
        std::vector<std::vector<int>> adj_matrix;
        std::vector<int> degree_sequence;
    public:
        // Works only for small graphs
        inline int GetEdgeLabel(int edge_id) const {return edge_label[edge_id];}

        inline int GetEdgeLabel(int u, int v) const {return adj_matrix[u][v];}

        inline std::vector<int>& GetDegreeSequence() {return degree_sequence;}

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
            num_edge = edges.size();
            vertex_label.resize(num_vertex);
            adj_list.resize(num_vertex);
            adj_matrix.resize(num_vertex);
            for (int i = 0; i < num_vertex; i++) {
                vertex_label[vertices[i]] = vertex_labels[i];
                vertex_label_frequency[vertex_labels[i]]++;
                adj_matrix[i].resize(num_vertex, -1);
            }
            edge_list = edges;
            edge_label = edge_labels;
            for (int i = 0; i < edge_list.size(); i++) {
                auto &[u, v] = edge_list[i];
                adj_matrix[u][v] = edge_label[i];
                adj_matrix[v][u] = edge_label[i];
                adj_list[u].push_back(v);
                adj_list[v].push_back(u);
                edge_label_frequency[edge_label[i]]++;
            }
            BuildDegreeSequence();
        };

        void CombineGraph(LabeledGraphDatabaseEntry *g1, LabeledGraphDatabaseEntry *g2) {
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

    };

    class LabeledGraphDatabase {
        std::vector<LabeledGraphDatabaseEntry*> data;
        std::vector<LabeledGraphDatabaseEntry*> queries;
        std::unordered_map<std::string, int> global_vertex_labels, global_edge_labels;
        int num_global_vertex_labels = 1, num_global_edge_labels = 1;
    public:
        void LoadGraphDatabase(std::string &filename, int which);
        std::vector<LabeledGraphDatabaseEntry*>& GetQueries() {return queries;}
        std::vector<LabeledGraphDatabaseEntry*>& GetData() {return data;}
        int GetNumGlobalVertexLabels() {return num_global_vertex_labels;}
        int GetNumGlobalEdgeLabels() {return num_global_edge_labels;}
        void WriteGraphDatabase(std::string path, std::string dataset);
    };


    void LabeledGraphDatabase::LoadGraphDatabase(std::string &filename, int which) {
        std::ifstream fin(filename);
    //    std::cerr << "Loading LabeledGraph database " << fileSize(filename) << " bytes" << std::endl;
        std::vector<int> graph_ids;
        std::vector<std::vector<int>> graph_vertices, graph_vertex_labels, graph_edge_labels;
        std::vector<std::vector<std::pair<int, int>>> graph_edges;
        while (!fin.eof()) {
            std::string raw;
            std::getline(fin, raw);
            auto line = parse(raw, " ");
            if (line[0] == "t") {
                graph_ids.push_back(std::stoi(line[2]));
                graph_vertices.emplace_back();
                graph_vertex_labels.emplace_back();
                graph_edge_labels.emplace_back();
                graph_edges.emplace_back();
            }
            else if (line[0] == "v") {
                int id, label;
                id = std::stoi(line[1]);
                string raw_label = line[2];
                if (global_vertex_labels.find(raw_label) == global_vertex_labels.end()) {
                    global_vertex_labels[raw_label] = num_global_vertex_labels++;
                }
                label = global_vertex_labels[raw_label];
                graph_vertices.back().push_back(id);
                graph_vertex_labels.back().push_back(label);
            }
            else if (line[0] == "e") {
                int u, v, label;
                u = std::stoi(line[1]);
                v = std::stoi(line[2]);
                string raw_label = line[3];
                if (global_edge_labels.find(raw_label) == global_edge_labels.end()) {
                    global_edge_labels[raw_label] = num_global_edge_labels++;
                }
                label = global_edge_labels[raw_label];
                graph_edges.back().emplace_back(u, v);
                graph_edge_labels.back().push_back(label);
            }
        }
    //    std::cerr << "Loaded " << graph_ids.size() << " graphs." << std::endl;
        for (int i = 0; i < graph_ids.size(); i++) {
            LabeledGraphDatabaseEntry *G = new LabeledGraphDatabaseEntry;
            G->LoadFromDatabase(graph_ids[i], graph_vertices[i], graph_vertex_labels[i],
                                 graph_edges[i], graph_edge_labels[i]);
            (which?data:queries).emplace_back(G);
        }
    //    printf("#Vertex Labels: %d\n", num_global_vertex_labels);
    //    printf("#Edge Labels: %d\n", num_global_edge_labels);
    }

    void LabeledGraphDatabase::WriteGraphDatabase(std::string path, std::string dataset) {
        fprintf(stderr, "Writing LabeledGraph database to %s\n", path.c_str());
        std::vector<LabeledGraphDatabaseEntry*> all_data = data;
        for (int i = 0; i < queries.size(); i++) all_data.emplace_back(queries[i]);
        std::ofstream edge_out(path+"/"+dataset+"_A.txt");
        std::ofstream edge_label_out(path+"/"+dataset+"_edge_labels.txt");
        std::ofstream vertex_label_out(path+"/"+dataset+"_node_labels.txt");
        std::ofstream vertex_out(path+"/"+dataset+"_graph_indicator.txt");
        int v_offset = 0;
        for (int i = 0; i < all_data.size(); i++) {
            LabeledGraph* g = all_data[i];
            for (int j = 0; j < g->GetNumEdges(); j++) {
                edge_out << g->GetEdge(j).first + 1 + v_offset << ", " << g->GetEdge(j).second + 1 + v_offset << std::endl;
                edge_label_out << g->GetEdgeLabel(j) - 1 << std::endl;
            }
            for (int j = 0; j < g->GetNumVertices(); j++) {
                vertex_label_out << g->GetVertexLabel(j) - 1 << std::endl;
                vertex_out << i + 1 << std::endl;
            }
            v_offset += g->GetNumVertices();
        }
    }
}