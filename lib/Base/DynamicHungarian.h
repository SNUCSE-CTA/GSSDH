#pragma once
#include "Base.h"
#include "PrettyPrint.h"
#include "Timer.h"

/**
 * @brief Hungarian algorithm for Assignment problem (Minimum-weighted bipartite
 * matching).
 * @param N number of vertices
 * @param cost_matrix cost matrix (edge weight)
 * @date April 09, 2024
 * @todo implement row-weight changes
 * @todo improve performance (not super-optimized)
 */

class Hungarian {
  const int INF = 1e9;
  std::vector<std::vector<int>> cost_matrix;
  std::vector<int> assignment, inverse_assignment, alpha, beta;
  std::vector<int> left_visited, right_visited;

  std::vector<int> changed_vertex; //for dynamic
  int total_cost, N, theta;
  bool dbg = 0;
  int64_t functioncall = 0;
  double hgtime = 0.0;

 public:
  void InitializeSolver() {
    total_cost = INF;
    assignment = std::vector<int>(N, -1);
    inverse_assignment = std::vector<int>(N, -1);
    left_visited = std::vector<int>(N, 0);
    right_visited = std::vector<int>(N, 0);
    alpha = std::vector<int>(N, 0);
    beta = std::vector<int>(N, 0);
    theta = 0;
  }

  Hungarian(int N_) {
    this->N = N_;
    InitializeSolver();
    cost_matrix = std::vector<std::vector<int>>(N, std::vector<int>(N, 0));
  }

  Hungarian(const std::vector<std::vector<int>>& cost_matrix_) {
    this->N = cost_matrix_.size();
    cost_matrix = cost_matrix_;
    InitializeSolver();
  }

  void FillCostMatrix(std::vector<std::vector<int>>& cost_matrix_) {
    N = cost_matrix_.size();
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        cost_matrix[i][j] = cost_matrix_[i][j];
      }
    }
  }

  void InitializeVariables() {
    //        fprintf(stderr, "Invoked Hungarian::InitializeVariables()\n");
    std::fill(alpha.begin(), alpha.end(), 0);
    std::fill(beta.begin(), beta.end(), INF);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        beta[j] = std::min(beta[j], cost_matrix[i][j]);
      }
    }
  }

  bool FindAugmentingPath(int i) {
    functioncall++;
    left_visited[i] = true;
    for (int j = 0; j < N; j++) {
      if (right_visited[j]) continue;
      if (alpha[i] + beta[j] != cost_matrix[i][j]) continue;
      //functioncall++;
      right_visited[j] = true;
      if (inverse_assignment[j] == -1 ||
          FindAugmentingPath(inverse_assignment[j])) {
        inverse_assignment[j] = i;
        assignment[i] = j;
        // functioncall++;
        return true;
      }
    }
    return false;
  }

  void RecalculatePotential() {
    // functioncall++;
    theta = INF;
    for (int i = 0; i < N; i++) {
      if (left_visited[i]) {
        for (int j = 0; j < N; j++) {
          if (!right_visited[j]) {
            theta = std::min(theta, cost_matrix[i][j] - alpha[i] - beta[j]);
            // functioncall++;
          }
        }
      }
    }
    for (int i = 0; i < N; i++) {
      if (left_visited[i]) {
        alpha[i] += theta;
      }
    }
    for (int j = 0; j < N; j++) {
      if (right_visited[j]) {
        beta[j] -= theta;
      }
    }
  }

 void Solve() {
    // double pt = 0.0;
    InitializeVariables();
    for(int i = 0; i < N; i++){
      // functioncall++;
      while (true) {

        std::fill(left_visited.begin(), left_visited.end(), 0);
        std::fill(right_visited.begin(), right_visited.end(), 0);
        Timer p;
        p.Start();
        bool flag = FindAugmentingPath(i);
        p.Stop();
        hgtime += p.GetTime(); 
        if (flag) break;                      
  
        // Timer p;
        // p.Start();
        RecalculatePotential();

        //     p.Stop();
        // hgtime += p.GetTime(); 
      }
    }
    total_cost = 0;
    for (int i = 0; i < N; i++) {
      // std::cout <<i << " " << assignment[i]  <<  " " <<cost_matrix[i][assignment[i]] <<"\n";
      total_cost += cost_matrix[i][assignment[i]];
    }
    // return pt;
  }

// void SolvePartial(int G1size, int *statemapping){
//   // std::cout << vertex << "\n";
//   InitializeVariables();
//   for(int i = 0 ; i < N;i++){
//     if(assignment[i] != -1) continue;
  //   while(true){
  //     std::fill(left_visited.begin(), left_visited.end(), 0);
  //     std::fill(right_visited.begin(), right_visited.end(), 0);
  //     if (FindAugmentingPath(i)) break;
  //     RecalculatePotential();
  //   }
  // }
//   total_cost = 0;
//     for (int i = 0; i < N; i++) {
//         if(i < G1size && statemapping[i] != -1) continue;
//         // std::cout << "u : "<< cost_matrix[i][assignment[i]] << " ";
//         total_cost += cost_matrix[i][assignment[i]];
//     }
//     // std::cout << "\n";
// }

void SolvePartial(){
  for(int i = 0 ; i < N; i++){
    if(assignment[i] != -1) continue;
      while(true){
        std::fill(left_visited.begin(), left_visited.end(), 0);
        std::fill(right_visited.begin(), right_visited.end(), 0);
        if (FindAugmentingPath(i)) break;
        RecalculatePotential();
    }
  }
  total_cost = 0;
    for (int i = 0; i < N; i++) {
      // std::cout <<i << " " << assignment[i]  <<  " " <<cost_matrix[i][assignment[i]] <<"\n";
      total_cost += cost_matrix[i][assignment[i]];
  }
}

// void Match(int u, int v){ //state->matrix //need to optimze
//   for(int j = 0 ; j < N; j++){
//     if(j != v){
//       ChangeCost(u, j, INF);
//     }
//   }
//   for(int i = 0 ; i < N; i++){
//     if(i != u){
//       ChangeCost(i, v, INF);
//     }
//   }
// }
void ChangeCost(int i, int j, int newCost){
    int vertex = i;
    int oldCost = cost_matrix[i][j];
    cost_matrix[i][j] = newCost;
    if(newCost > oldCost && assignment[i] == j){
        assignment[i] = -1;
        inverse_assignment[j] = -1;
    }
    else if(newCost < oldCost && alpha[i] + beta[j] > newCost){
        alpha[i] = INF;
        for(int j = 0 ; j < N; j++){
            alpha[i] = std::min(alpha[i], cost_matrix[i][j] - beta[j]);
        }
        if(assignment[i] != j){
            inverse_assignment[assignment[i]] = -1;
            assignment[i] = -1;
        }
    }
}

  void ChangeRow(int i, std::vector<int>& newRow){
    int vertex = i;
    inverse_assignment[assignment[i]] = -1;
    assignment[i] = -1;
    alpha[i] = INF;
    for(int j = 0; j < N; j++){
        cost_matrix[i][j] = newRow[j];
    }
    for(int j = 0 ; j < N ; j++){
        alpha[i] = std::min(alpha[i], cost_matrix[i][j] - beta[j]);
    }
    while (true) {
        std::fill(left_visited.begin(), left_visited.end(), 0);
        std::fill(right_visited.begin(), right_visited.end(), 0);
        if (FindAugmentingPath(vertex)) break;
        RecalculatePotential();
    }
}

  void ChangeCol(int j, std::vector<int>& newCol){
    int vertex = inverse_assignment[j];
    assignment[inverse_assignment[j]] = -1;
    inverse_assignment[j] = -1;
    beta[j] = INF;
    for(int i = 0;i < N; i++){
        cost_matrix[i][j] = newCol[i];
    }
    for(int i = 0;i < N; i++){
        beta[j] = std::min(beta[j], cost_matrix[i][j] - alpha[i]);
    }
    while (true) {
        std::fill(left_visited.begin(), left_visited.end(), 0);
        std::fill(right_visited.begin(), right_visited.end(), 0);
        if (FindAugmentingPath(vertex)) break;
        RecalculatePotential();
    }
  }

  void CopyAssignmentAlphaBeta(std::vector<int>& parent_assignment, std::vector<int>&parent_inverse_assignment, std::vector<int>&parent_alpha, std::vector<int>&parent_beta, std::vector<int>& rem_left, std::vector<int>& rem_right){
    for(int u_idx = 0; u_idx < rem_left.size(); u_idx++){
      int u_curr = rem_left[u_idx];
      assignment[u_idx] = parent_assignment[u_curr];
      alpha[u_idx] = parent_alpha[u_curr];
    }
    for(int v_idx = 0 ; v_idx < rem_right.size();v_idx++){
      int v_curr = rem_right[v_idx];
      inverse_assignment[v_idx] = parent_inverse_assignment[v_curr];
      beta[v_idx] = parent_beta[v_curr];
    }
    // std::cout << parent_assignment << "\n";
    // std::cout << assignment << "\n";
    std::cout << parent_beta << "\n";
    std::cout << beta << "\n";
  }
  void CopyAssignmentAlphaBetaReverse(std::vector<int>& parent_assignment, std::vector<int>&parent_inverse_assignment, std::vector<int>&parent_alpha, std::vector<int>&parent_beta, std::vector<int>& rem_left, std::vector<int>& rem_right){
    for(int u_idx = 0; u_idx < rem_left.size(); u_idx++){
      int u_curr = rem_left[u_idx];
      parent_assignment[u_curr] = assignment[u_curr];
      parent_alpha[u_curr] = alpha[u_idx];
    }
    for(int v_idx = 0 ; v_idx < rem_right.size();v_idx++){
      int v_curr = rem_right[v_idx];
      parent_inverse_assignment[v_curr] = inverse_assignment[v_idx];
      parent_beta[v_curr] = beta[v_curr];
    }

  }

  void Print() {
    fprintf(stderr, "Invoked Hungarian::Print() N = %d\n", N);
    for (int i = 0; i < N; i++) {
      std::cerr << cost_matrix[i] << '\n';
    }
    fprintf(stderr, "Assignment:\n");
    for (int i = 0; i < N; i++) {
      fprintf(stderr, "  %d(cost=%d)\n", assignment[i],
              cost_matrix[i][assignment[i]]);
    }
    fprintf(stderr, "TotalCost = %d\n", total_cost);
  }

  std::vector<int>& GetAssignment() { return assignment; }

  std::vector<std::vector<int>>& GetMatrix(){return cost_matrix; }
  int GetTotalCost() { return total_cost; }
  int64_t GetCnt(){return functioncall; }
  double GetTime(){return hgtime; }

  int AssignedWeight(int i) { return cost_matrix[i][assignment[i]]; }

  void PrintAssignment(){
    std::cout <<"Assignment\n";
    for(int i = 0 ; i < N; i++){
      std::cout << assignment[i] << " ";
    }
    std::cout << "\nInverse Assignment\n";
    for(int i = 0 ; i < N;i++){
      std::cout << inverse_assignment[i] << " ";
    }
    std::cout << "\n";
  }
  void PrintMatrix(){
    std::cout <<"hungarian class matrix\n";
    for(int i = 0 ; i < N; i++){
      for(int j = 0 ; j < N; j++){
        std::cout << cost_matrix[i][j] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "-------------------------\n";
  }
};