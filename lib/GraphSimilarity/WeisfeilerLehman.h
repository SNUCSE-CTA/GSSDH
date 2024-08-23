// #pragma once
// #include <map>
// #include <set>
// #include <unordered_map>

// #include "GraphSimilaritySearch.h"
// const bool DEBUG = 0;
// const int OFFSET = 1000000;
// const int NUM_MAX_COLOR = 1000;
// namespace GraphLib::GraphSimilarity {
// class WeisfeilerLehman {
//  public:
//   GSSEntry *G;
//   int *color, *tmp_color, *perm, *num_color_childs, *num_left, *num_right;
//   int *degree;
//   int *ignore;
//   int num_colors = 0;
//   int *color_parent;

//   void InitializeColor(int num_fixed_color = 0) {
//     std::iota(perm, perm + G->GetNumVertices(), 0);
//     std::sort(perm, perm + G->GetNumVertices(), [&](int i, int j) -> bool {
//       if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
//       int di = degree[i], dj = degree[j];
//       int li = G->GetVertexLabel(i), lj = G->GetVertexLabel(j);
//       // return (di < dj) or (di == dj && li < lj);
//       return (li < lj);
//     });
//     int current_color = num_fixed_color + 1;
//     if (DEBUG) {
//       printf("Initial perm order: ");
//       for (int i = 0; i < G->GetNumVertices(); i++) {
//         printf("%d ", perm[i]);
//       }
//       printf("\n");
//     }
//     color[perm[0]] = current_color;
//     for (int i = 1; i < G->GetNumVertices(); i++) {
//       if (ignore[perm[i]]) {
//         color[perm[i]] = ignore[perm[i]];
//         continue;
//       }
//       int u = perm[i], v = perm[i - 1];
//       int du = degree[u], dv = degree[v];
//       int lu = G->GetVertexLabel(u), lv = G->GetVertexLabel(v);
//       if (du != dv or lu != lv) {
//         // printf("%d is assigned different color (%d) from/%d\n", u,
//         //        current_color + 1, v);
//         current_color++;
//       }
//       if (lu != lv) {
//         current_color++;
//       }
//       color[u] = current_color;
//     }
//     num_colors = current_color + 1;
//     color_parent[0] = -1;
//     for (int i = 1; i < num_colors; i++) {
//       color_parent[i] = 0;
//       num_color_childs[0]++;
//     }
//   }
//   WeisfeilerLehman(GSSEntry *g) {
//     G = g;
//     degree = new int[G->GetNumVertices()]();
//     for (int i = 0; i < G->GetNumVertices(); i++) {
//       degree[i] = G->GetDegree(i);
//     }
//     ignore = new int[G->GetNumVertices()]();
//     color = new int[G->GetNumVertices()]();
//     tmp_color = new int[G->GetNumVertices()]();
//     // Give initial color with (degree and label) information
//     perm = new int[G->GetNumVertices()]();
//     color_parent = new int[NUM_MAX_COLOR]();
//     num_color_childs = new int[NUM_MAX_COLOR]();
//     num_left = new int[NUM_MAX_COLOR]();
//     num_right = new int[NUM_MAX_COLOR]();
//     memset(color_parent, -1, sizeof(int) * NUM_MAX_COLOR);
//     InitializeColor();
//   }

//   void Reset(std::vector<int> &current_mapping) {
//     for (int i = 0; i < G->GetNumVertices(); i++) {
//       degree[i] = G->GetDegree(i);
//     }
//     std::memset(ignore, false, sizeof(int) * G->GetNumVertices());
//     int VL = current_mapping.size();
//     int num_ignore = 0;
//     for (int i = 0; i < VL; i++) {
//       if (current_mapping[i] != -1) {
//         ignore[i] = ignore[current_mapping[i] + VL] = (++num_ignore);
//       }
//     }
//     std::memset(color, 0, sizeof(int) * G->GetNumVertices());
//     std::memset(tmp_color, 0, sizeof(int) * G->GetNumVertices());
//     std::memset(color_parent, -1, sizeof(int) * num_colors);
//     std::memset(num_color_childs, 0, sizeof(int) * num_colors);
//     std::memset(num_left, 0, sizeof(int) * num_colors);
//     std::memset(num_right, 0, sizeof(int) * num_colors);
//     InitializeColor(num_ignore - 1);
//   }

//   void RefineStep() {
//     int num_diff_neighbors = 0;
//     memcpy(tmp_color, color, sizeof(int) * G->GetNumVertices());
//     std::map<std::vector<int>, int> neighbor_hash;
//     for (int u = 0; u < G->GetNumVertices(); u++) {
//       if (ignore[u]) continue;
//       std::vector<int> h = {color[u]};
//       for (int v : G->GetNeighbors(u)) {
//         h.push_back(color[v] + G->GetEdgeLabel(u, v) * OFFSET);
//       }
//       std::sort(h.begin(), h.end());
//       if (neighbor_hash.find(h) == neighbor_hash.end()) {
//         neighbor_hash[h] = num_diff_neighbors++;
//       }
//       tmp_color[u] = num_colors + neighbor_hash[h];
//       if (color_parent[tmp_color[u]] == -1) {
//         color_parent[tmp_color[u]] = color[u];
//         num_color_childs[color[u]]++;
//       }
//     }
//     memcpy(color, tmp_color, sizeof(int) * G->GetNumVertices());
//     num_colors = num_colors + num_diff_neighbors;
//   }

//   void Refine(int iterations = 2) {
//     for (int t = 0; t < iterations; t++) {
//       RefineStep();
//     }
//     if (DEBUG) {
//       for (int i = 0; i < G->GetNumVertices(); i++) {
//         printf("Vertex %d: Color = %d (before = %d)\n", i, color[i],
//                color_parent[color[i]]);
//       }
//     }
//   }

//   void Match(std::vector<int> &mapping) {
//     int VL = mapping.size();
//     int VR = G->GetNumVertices() - VL;
//     int mapped_idx = 0;
//     std::queue<int> color_queue;
//     std::vector<std::deque<int>> color_partition(num_colors);
//     std::vector<std::vector<int>> mapped_neighbors(G->GetNumVertices());
//     if (DEBUG) {
//       printf("Start mapping...");
//       for (int i = 0; i < VL; i++) {
//         printf("%d ", mapping[i]);
//       }
//       printf("\n");
//     }
//     for (int i = 0; i < G->GetNumVertices(); i++) {
//       if (ignore[i]) continue;
//       if (DEBUG)
//         fprintf(stdout, "Vertex %d: Start with color = %d\n", i, color[i]);
//       if (i >= VL) {
//         color_partition[color[i]].push_back(i);
//         num_right[color[i]]++;
//       } else {
//         color_partition[color[i]].push_front(i);
//         num_left[color[i]]++;
//       }
//     }
//     for (int i = 0; i < num_colors; i++) {
//       if (DEBUG) printf("Color %d: Parent %d\n", i, color_parent[i]);
//       if (num_color_childs[i] == 0) {
//         color_queue.push(i);
//       }
//     }
//     while (!color_queue.empty()) {
//       int c = color_queue.front();
//       if (DEBUG)
//         printf("Check color %d (L, R = %d, %d)\n", c, num_left[c],
//                num_right[c]);
//       color_queue.pop();
//       if (num_left[c] > 0 and num_right[c] > 0) {
//         std::sort(color_partition[c].begin(), color_partition[c].end(),
//                   [&](int a, int b) -> bool {
//                     if (a < VL and b < VL)
//                       return mapped_neighbors[a] > mapped_neighbors[b];
//                     else if (a < VL and b >= VL)
//                       return true;
//                     else if (a >= VL and b < VL)
//                       return false;
//                     else
//                       return mapped_neighbors[a] < mapped_neighbors[b];
//                   });
//         while (num_left[c] > 0 and num_right[c] > 0) {
//           int l = color_partition[c].front();
//           int r = color_partition[c].back();
//           color_partition[c].pop_front();
//           color_partition[c].pop_back();
//           mapping[l] = r - VL;
//           // mapped_idx++;
//           // for (int lnbr : G->GetNeighbors(l)) {
//           //   mapped_neighbors[lnbr].push_back(mapped_idx);
//           // }
//           // for (int rnbr : G->GetNeighbors(r)) {
//           //   mapped_neighbors[rnbr].push_back(mapped_idx);
//           // }
//           if (DEBUG) printf("!!Map %d to %d (%d)\n", l, r, r - VL);
//           num_left[c]--;
//           num_right[c]--;
//         }
//       }

//       int p = color_parent[c];
//       if (p == -1) continue;
//       for (auto x : color_partition[c]) {
//         if (x >= VL) {
//           color_partition[p].push_back(x);
//           num_right[p]++;
//         } else {
//           color_partition[p].push_front(x);
//           num_left[p]++;
//         }
//         if (DEBUG) printf("move vertex %d to color %d\n", x, p);
//       }
//       num_left[c] = num_right[c] = 0;
//       num_color_childs[p]--;
//       if (num_color_childs[p] == 0) {
//         color_queue.push(p);
//       }
//     }
//   }
// };
// }  // namespace GraphLib::GraphSimilarity
