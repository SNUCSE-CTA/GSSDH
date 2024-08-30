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
    left_visited[i] = true;
    for (int j = 0; j < N; j++) {
      if (right_visited[j]) continue;
      if (alpha[i] + beta[j] != cost_matrix[i][j]) continue;
      right_visited[j] = true;
      if (inverse_assignment[j] == -1 ||
          FindAugmentingPath(inverse_assignment[j])) {
        inverse_assignment[j] = i;
        assignment[i] = j;
        return true;
      }
    }
    return false;
  }

  void RecalculatePotential() {
    theta = INF;
    for (int i = 0; i < N; i++) {
      if (left_visited[i]) {
        for (int j = 0; j < N; j++) {
          if (!right_visited[j]) {
            theta = std::min(theta, cost_matrix[i][j] - alpha[i] - beta[j]);
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
    InitializeVariables();
    for (int i = 0; i < N; i++) {
      while (true) {
        std::fill(left_visited.begin(), left_visited.end(), 0);
        std::fill(right_visited.begin(), right_visited.end(), 0);
        if (FindAugmentingPath(i)) break;
        RecalculatePotential();
      }
    }
    total_cost = 0;
    for (int i = 0; i < N; i++) {
      total_cost += cost_matrix[i][assignment[i]];
    }
  }

void SolvePartial(int vertex){
  // std::cout << vertex << "\n";
  while(true){
      std::fill(left_visited.begin(), left_visited.end(), 0);
      std::fill(right_visited.begin(), right_visited.end(), 0);
      if (FindAugmentingPath(vertex)) break;
      RecalculatePotential();
  }
}

void TotalCost(int*isMapping){
  total_cost = 0;
  for(int i = 0 ; i < N; i++){
    if(isMapping[i] == -1){
    total_cost += cost_matrix[i][assignment[i]];
    }
  }
}
int ChangeCost(int i, int j, int newCost){
    int vertex = i;
    int oldCost = cost_matrix[i][j];
    cost_matrix[i][j] = newCost;
    if(newCost > oldCost && assignment[i] == j){
        assignment[i] = -1;
        inverse_assignment[j] = -1;
      //   while (true) {
      //       std::fill(left_visited.begin(), left_visited.end(), 0);
      //       std::fill(right_visited.begin(), right_visited.end(), 0);
      //       if (FindAugmentingPath(vertex)) break;
      //       RecalculatePotential();
      // }
        return vertex;
    }
    else if(newCost < oldCost && alpha[i] + beta[j] > newCost){
        alpha[i] = INF;
        for(int j = 0 ; j < N; j++){
            alpha[i] = std::min(alpha[i], cost_matrix[i][j] - beta[j]);
        }
        if(assignment[i] != j){
            inverse_assignment[assignment[i]] = -1;
            assignment[i] = -1;
        // while (true) {
        //     std::fill(left_visited.begin(), left_visited.end(), 0);
        //     std::fill(right_visited.begin(), right_visited.end(), 0);
        //     if (FindAugmentingPath(vertex)) break;
        //     RecalculatePotential();
        // }
        return vertex;
        }
        else return -1;
    }
    else{
      return -1;
    }
}

int Match(int u, int v, std::vector<std::vector<int>>& matrix){ //state->matrix

  int u_prime = inverse_assignment[v];
  int v_prime = assignment[u];

  if(u_prime != u){
    inverse_assignment[v_prime] = -1;
    assignment[u_prime] = -1;
  }
    assignment[u] = v;
    inverse_assignment[v] = u;
  for(int j = 0 ; j < N; j++){
    if(j != v){
      cost_matrix[u][j] = INF;
      matrix[u][j] = INF;
    }
  }
  for(int i = 0 ; i < N; i++){
    if(i != u){
      cost_matrix[i][v] = INF;
      matrix[i][v] = INF;
    }
  }
  return u_prime;
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

  void CopyParent(std::vector<int>& parent_assignment){
    assignment = parent_assignment;
    inverse_assignment.resize(parent_assignment.size(), -1);
    for (int i = 0; i < assignment.size(); i++) {
        inverse_assignment[assignment[i]] = i;
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

  int GetTotalCost() { return total_cost; }

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