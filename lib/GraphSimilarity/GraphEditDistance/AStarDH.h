#pragma once
// #include "Base/Hungarian.h"
// #include "Base/DynamicHungarian.h"
#include "GraphSimilarity/EditDistance.h"
#include "Base/Timer.h"
#include "Base/DynamicHungarian.h"

namespace GraphLib::GraphSimilarity {
struct DHState : State {
  int next_mapping_order = -1;
  int ub = 1e9;
  std::vector<std::vector<int>> matrix;
  std::vector<int> hungarian_assignment;
  std::vector<int> alpha, beta;
};


class AStarDH : public GraphEditDistanceSolver {
  const bool DEBUG = false;
  const bool dbg = true;
  int cnt = 0;
  std::priority_queue<DHState*, std::vector<DHState*>, StateComparator> queue;
  std::map<int, int> num_hungarian;

  // Timer hg_timer;
  double hungarian_time = 0.0, branchdistance_time = 0.0;
  int64_t hungarian_vertex_num = 0;
/*from hungarian*/
  const int INF = 1e9;
  std::vector<int> assignment, inverse_assignment;
  std::vector<int> left_visited, right_visited;
  int total_cost, N, theta;

 public:

/*from hungarian*/
void InitializeSolve_(std::vector<int>& alpha, std::vector<int>& beta) {
    total_cost = INF;
    assignment = std::vector<int>(N, -1);
    inverse_assignment = std::vector<int>(N, -1);
    // std::cerr<<__LINE__<<"\n";
    left_visited = std::vector<int>(N, 0);
    right_visited = std::vector<int>(N, 0);
    alpha = std::vector<int>(N, 0);
    beta = std::vector<int>(N, 0);
    theta = 0;
        // std::cerr<<__LINE__<<"\n";
}

void InitializeVariables_(std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta) {
    //        fprintf(stderr, "Invoked Hungarian::InitializeVariables()\n");
    std::fill(alpha.begin(), alpha.end(), 0);
    std::fill(beta.begin(), beta.end(), INF);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        beta[j] = std::min(beta[j], cost_matrix[i][j]);
      }
    }
  }
  // bool FindAugmentingPath(int i, std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta) {
  //   left_visited[i] = true;
  //   for (int j = 0; j < N; j++) {
  //     if (right_visited[j]) continue;
  //     if (alpha[i] + beta[j] != cost_matrix[i][j]) continue;
  //     right_visited[j] = true;
  //     if (inverse_assignment[j] == -1 ||
  //         FindAugmentingPath(inverse_assignment[j], cost_matrix, alpha, beta)) {
  //       inverse_assignment[j] = i;
  //       assignment[i] = j;
  //       return true;
  //     }
  //   }
  //   return false;
  // }

  bool FindAugmentingPath_(int i, std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta){
    left_visited[i] = true;
    for(int j = 0 ; j < N; j++){
        if(!right_visited[j] && (alpha[i] + beta[j] == cost_matrix[i][j])){
            right_visited[j] = true;
            if(inverse_assignment[j] == -1 || FindAugmentingPath_(inverse_assignment[j], cost_matrix, alpha, beta)){
                inverse_assignment[j] = i;
                assignment[i] = j;
                return true;
            }
        }
    }
    return false;
}
  void RecalculatePotential_(std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta) {
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

 bool FindAugmentingPath_p(int i, std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta, std::vector<int> &remain_right){
  auto &matrix = cost_matrix;
    left_visited[i] = true;
    for(auto j : remain_right){
        if(!right_visited[j] && (alpha[i] + beta[j] == matrix[i][j])){
            right_visited[j] = true;
            if(inverse_assignment[j] == -1 || FindAugmentingPath_p(inverse_assignment[j], cost_matrix, alpha, beta, remain_right)){
                inverse_assignment[j] = i;
                assignment[i] = j;
                return true;
            }
        }
    }
    return false;
}
  void RecalculatePotential_p(std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta, std::vector<int> &remain_left, std::vector<int> &remain_right) {
    auto &matrix = cost_matrix; 
    theta = INF;
    for (auto i : remain_left) {
      if (left_visited[i]) {
        int alpha_i = alpha[i];
        for (auto j : remain_right) {
          if (!right_visited[j]) {
            theta = std::min(theta, matrix[i][j] - alpha_i - beta[j]);
          }
        }
      }
    }
    for (auto i : remain_left) {
      if (left_visited[i]) {
        alpha[i] += theta;
      }
    }
    for (auto j : remain_right) {
      if (right_visited[j]) {
        beta[j] -= theta;
      }
    }
  }

void SolvePartial_p(DHState *state, std::vector<int> &remain_left, std::vector<int> &remain_right){
  auto &matrix_ = state->matrix;
  auto &alpha_ = state->alpha;
  auto &beta_ = state->beta;
  for(auto i : remain_left){
    if(assignment[i] != -1) continue;
    hungarian_vertex_num++;
    while(true){
      std::fill(left_visited.begin(), left_visited.end(), 0);
      std::fill(right_visited.begin(), right_visited.end(), 0);
      if (FindAugmentingPath_p(i, matrix_, alpha_, beta_, remain_right)) break;
      RecalculatePotential_p(matrix_, alpha_, beta_, remain_left, remain_right);
    }
  }


  total_cost = 0;
  for (int i = 0; i < G1->GetNumVertices(); i++){
        if(state->mapping[i] == -1){
            total_cost += state->matrix[i][assignment[i]];
        }
    }
  for(int i = G1->GetNumVertices() ; i < G2->GetNumVertices();i++){
        total_cost += state->matrix[i][assignment[i]];
    }
}


void Solve_(std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha,std::vector<int>& beta) {
    InitializeVariables_(cost_matrix, alpha, beta);
    for (int i = 0; i < N; i++) {
      hungarian_vertex_num++;
      while (true) {
        std::fill(left_visited.begin(), left_visited.end(), 0);
        std::fill(right_visited.begin(), right_visited.end(), 0);
        if (FindAugmentingPath_(i, cost_matrix, alpha, beta)) break;
        RecalculatePotential_(cost_matrix, alpha, beta);
      }
    }
    total_cost = 0;
    for (int i = 0; i < N; i++) {
      total_cost += cost_matrix[i][assignment[i]];
    }
  }


void SolvePartial_(DHState *state){
  for(int i = 0 ; i < N; i ++){
    if(assignment[i] != -1) continue;
    hungarian_vertex_num++;
    while(true){
      
      std::fill(left_visited.begin(), left_visited.end(), 0);
      std::fill(right_visited.begin(), right_visited.end(), 0);
      if (FindAugmentingPath_(i, state->matrix, state->alpha, state->beta)) break;
      RecalculatePotential_(state->matrix,state->alpha, state->beta);
    }
  }
  total_cost = 0;
  for (int i = 0; i < G1->GetNumVertices(); i++){
        if(state->mapping[i] == -1){
            total_cost += state->matrix[i][assignment[i]];
            // std::cout << i << " " << assignment[i] << " " <<state->matrix[i][assignment[i]] << "\n";
        }
    }
  for(int i = G1->GetNumVertices() ; i < G2->GetNumVertices();i++){
        total_cost += state->matrix[i][assignment[i]];
    }
    // std::cout << "\n";
}


void Match_(DHState *state){
    int u = matching_order[state->depth];
    int v = state->mapping[u];
    for(int j = 0; j < N; j++){
      if(j == v) continue;
        ChangeCost(u, j, INF, state->matrix, state->alpha, state->beta, false);
    }
    for(int i = 0 ; i < N; i++){
      if(i == u) continue;
      ChangeCost(i, v, INF, state->matrix,state->alpha, state->beta, false);
    }
      inverse_assignment[v] = u;
      assignment[u] = v;
}

void ChangeCost(int i, int j, int newCost, std::vector<std::vector<int>>& cost_matrix, std::vector<int>& alpha, std::vector<int>& beta, bool rowFlag){
    int vertex = i;
    int oldCost = cost_matrix[i][j];
    cost_matrix[i][j] = newCost;
    if(newCost > oldCost && assignment[i] == j){  
        assignment[i] = -1;
        inverse_assignment[j] = -1;
    }
    // else if(newCost < oldCost && alpha[i] + beta[j] > newCost){
    //     alpha[i] = INF;
    //     if(rowFlag == true){
    //       for(int j = 0 ; j < N; j++){
    //         alpha[i] = std::min(alpha[i], cost_matrix[i][j] - beta[j]);
    //       }
    //     }
    //     if(assignment[i] != j){
    //         inverse_assignment[assignment[i]] = -1;
    //         assignment[i] = -1;
    //     }
    // }
}

  void ExtendState(DHState* state) {
    if (state->cost >= current_best) return;
    int depth = state->depth + 1;
    int u = matching_order[depth];
    for (int v = 0; v < G2->GetNumVertices(); v++) {
      if (state->inverse_mapping[v] != -1) continue;
      DHState* child_state = new DHState(state);
      child_state->matrix = state->matrix;//copy matrix from parent

      child_state->cost = GetChildEditCost(state, u, v);
      if (DEBUG) {
        fprintf(stderr, "Parent %d has cost %d, child %d has cost %d\n",
                state->id, state->cost, child_state->id, child_state->cost);
      }
      child_state->mapping[u] = v;
      child_state->inverse_mapping[v] = u;
      child_state->depth = depth;
      auto [lb, ub] = DHLowerBound(child_state);
      child_state->lower_bound = lb;
      child_state->ub = ub;
      if (lb == ub) continue;
      if (depth == G1->GetNumVertices() - 1) {
        if (DEBUG)
          fprintf(stderr, "Reached at leaf of cost %d\n",
                  child_state->lower_bound);
        current_best = std::min(current_best, child_state->lower_bound);
        memcpy(&current_best_mapping[0], child_state->mapping,
               sizeof(int) * NumG1Vertices);
        continue;
      }
      if (child_state->lower_bound >= current_best) continue;
      if (threshold > 0) {
        if (child_state->lower_bound > threshold) continue;
      }

      if (DEBUG) {
        fprintf(stderr, "Parent %d pushes state %d with bound %d\n", state->id,
                child_state->id, child_state->lower_bound);
        child_state->Print();
      }
      queue.push(child_state);
    }
  }

  int GED() {
    PrepareGED();
    DHState* initial_state = new DHState(NULL);
    initial_state->cost = 0;
    initial_state->depth = -1;

    hungarian_time = 0.0;
    branchdistance_time = 0.0;
    auto [lb, ub] = DHLowerBound(initial_state);

    initial_state->lower_bound = lb;
    queue.push(initial_state);
    int64_t max_qsize = 1;
    while (!queue.empty()) {
      DHState* current_state = queue.top();
      if (DEBUG) current_state->Print();
      num_nodes++;
      if (num_nodes % LOG_EVERY == 0) {
        fprintf(stderr, "GED of %d and %d\n", G1->GetId(), G2->GetId());
        fprintf(stderr,
                "%llu, Current best = %d, current cost = %d / lb = %d, depth "
                "%d, Queuesize %lu\n",
                num_nodes, current_best, current_state->cost,
                current_state->lower_bound, current_state->depth, queue.size());
      }
      queue.pop();
      if (current_state->lower_bound >= current_best) {
        queue = std::priority_queue<DHState*, std::vector<DHState*>,
                                    StateComparator>();
        break;
      }
      if (threshold >= 0) {
        if (current_best < threshold) {
          queue = std::priority_queue<DHState*, std::vector<DHState*>,
                                      StateComparator>();
          break;
        }
        if (current_state->lower_bound > threshold) {
          current_best = -1;
          queue = std::priority_queue<DHState*, std::vector<DHState*>,
                                      StateComparator>();
          break;
        }
      }
      ExtendState(current_state);
      max_qsize = std::max(max_qsize, (int64_t)queue.size());
    }
    if (threshold >= 0 and current_best > threshold) {
      current_best = -1;
    }
    log.AddResult("MaxQueueSize", max_qsize, RESULT_INT64);
    log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
    log.AddResult("EditDistance", current_best, RESULT_INT);
    return current_best;
  }

  const int TARGET_ID = -2;

  void ComputeBranchDistanceMatrixInitial(DHState *state){
    DifferenceVector diff;
    diff.init(20);
    int N = state->matrix.size();
    int colSize = G2->GetNumVertices(); //right
    int rowSize = G1->GetNumVertices(); //left
    for(int v = 0 ; v < colSize; v++){
        auto &v_nbrs = G2->GetNeighbors(v);
        int u = 0;
        for(u = 0 ; u < rowSize; u++){
            auto &u_nbrs = G1->GetNeighbors(u);
            diff.reset();
            if(G1->GetVertexLabel(u) != G2->GetVertexLabel(v)){
                state->matrix[u][v] += 2;
            }
            for(int l = 0 ; l < u_nbrs.size();l++){
                int u_nbr = u_nbrs[l];
                int u_el = G1->GetEdgeLabel(u, u_nbr);
                if(state->mapping[u_nbr] == -1){
                    diff.update(u_el, 1);
                }
                else{
                    int v_mapping_el = G2->GetEdgeLabel(v, state->mapping[u_nbr]);
                    if(v_mapping_el != u_el){
                        state->matrix[u][v] += 2;
                    }
                }
            }
            for(int r = 0 ; r < v_nbrs.size(); r++){
                int v_nbr = v_nbrs[r];
                int v_el = G2->GetEdgeLabel(v, v_nbr);
                if(state->inverse_mapping[v_nbr] == -1){
                    diff.update(v_el, -1);
                }
                else{
                    int u_mapping_el = G1->GetEdgeLabel(u, state->inverse_mapping[v_nbr]);
                    if(u_mapping_el == -1){
                        state->matrix[u][v] += 2;
                    }
                }
            }
            int inner_distance = diff.GetDifference();
            state->matrix[u][v] += inner_distance;
        }
        int from_null = BranchEditDistanceFromNull(G2->GetBranch(v));
        for(int v_nbr : G2->GetNeighbors(v)){
            if(state->inverse_mapping[v_nbr] != -1){
                from_null++;
            }
        }
        for(; u < colSize; u++){
            state->matrix[u][v] = from_null;
        }
    }
  }

//   std::vector<std::tuple<int,int,int>> ComputeBranchDistanceMatrixDynamic(DHState *state){ //return u, v, cost
    void ComputeBranchDistanceMatrixDynamic(DHState *state){ 
        DifferenceVector diff;
        diff.init(20);
        // std::vector<std::tuple<int,int,int>> vertex;
        int u = matching_order[state->depth];
        int v = state->mapping[u];
        auto &u_nbrs = G1->GetNeighbors(u);
        auto &v_nbrs = G2->GetNeighbors(v);
        int newCost = 0;
        std::vector<bool> u_visited(G1->GetNumVertices(), 0);
        for(int i = 0; i < u_nbrs.size(); i++){
            int u_curr = u_nbrs[i];
            u_visited[u_curr] = true; //do not compute twice
            if(state->mapping[u_curr] != -1) continue; //already mapped vertex shoud not change weight
            for(int j = 0 ; j < G2->GetNumVertices(); j++){
                int v_curr = j;
                if(state->inverse_mapping[v_curr] != -1) continue;
                
                newCost = 0;

                auto &u_curr_nbrs = G1->GetNeighbors(u_curr);
                diff.reset();
                if(G1->GetVertexLabel(u_curr) != G2->GetVertexLabel(v_curr)){
                    newCost += 2;
                }
                for(int l = 0; l < u_curr_nbrs.size(); l++){
                    int u_curr_nbr = u_curr_nbrs[l];
                    int u_curr_el = G1->GetEdgeLabel(u_curr, u_curr_nbr);
                    if(state->mapping[u_curr_nbr] == -1){
                        diff.update(u_curr_el, 1);
                    }
                    else{
                        int v_curr_mapping_el = G2->GetEdgeLabel(v_curr, state->mapping[u_curr_nbr]);
                        if(v_curr_mapping_el != u_curr_el){
                            newCost += 2;
                        }
                    }
                }
                auto &v_curr_nbrs = G2->GetNeighbors(v_curr);
                for(int r = 0; r < v_curr_nbrs.size(); r++){
                    int v_curr_nbr = v_curr_nbrs[r];
                    int v_curr_el = G2->GetEdgeLabel(v_curr, v_curr_nbr);
                    if(state->inverse_mapping[v_curr_nbr] == -1){
                        diff.update(v_curr_el, -1);
                    }
                    else{
                        int u_curr_mapping_el = G1->GetEdgeLabel(u_curr, state->inverse_mapping[v_curr_nbr]);
                        if(u_curr_mapping_el == -1){
                            newCost+=2;
                        }
                    }
                }
                int inner_distance = diff.GetDifference();
                newCost += inner_distance;

                // state->matrix[u_curr][v_curr] = newCost;
                if(j == G2->GetNumVertices() -1){
                  ChangeCost(u_curr, v_curr, newCost, state->matrix, state->alpha, state->beta, true);
                }
                else{
                  ChangeCost(u_curr, v_curr, newCost, state->matrix, state->alpha, state->beta, false);
                }
                // vertex.push_back(std::make_tuple(u_curr, v_curr, newCost));
            }
            // ComputeAlpha(state, u_curr);
        }
        newCost = 0;
        for(int j = 0 ; j < v_nbrs.size(); j++){
            int v_curr = v_nbrs[j];
            if(state->inverse_mapping[v_curr] != -1) continue;
            for(int i = 0 ; i < G1->GetNumVertices(); i++){
                int u_curr = i;
                if(state->mapping[u_curr] != -1 || u_visited[u_curr] == true) continue;
                newCost = 0;
                auto v_curr_nbrs = G2->GetNeighbors(v_curr);
                diff.reset();
                if(G1->GetVertexLabel(u_curr) != G2->GetVertexLabel(v_curr)){
                    newCost +=2;
                }
                auto &u_curr_nbrs = G1->GetNeighbors(u_curr);
                for(int l = 0; l < u_curr_nbrs.size(); l++){
                    int u_curr_nbr = u_curr_nbrs[l];
                    int u_curr_el = G1->GetEdgeLabel(u_curr, u_curr_nbr);
                    if(state->mapping[u_curr_nbr] == -1){
                        diff.update(u_curr_el, 1);
                    }
                    else{
                        int v_curr_mapping_el = G2->GetEdgeLabel(v_curr, state->mapping[u_curr_nbr]);
                        if(v_curr_mapping_el != u_curr_el){
                            newCost += 2;
                        }
                    }
                }
                for(int r = 0; r < v_curr_nbrs.size(); r++){
                    int v_curr_nbr = v_curr_nbrs[r];
                    int v_curr_el = G2->GetEdgeLabel(v_curr, v_curr_nbr);
                    if(state->inverse_mapping[v_curr_nbr] == -1){
                        diff.update(v_curr_el, -1);
                    }
                    else{
                        int u_curr_mapping_el = G1->GetEdgeLabel(u_curr, state->inverse_mapping[v_curr_nbr]);
                        if(u_curr_mapping_el == -1){
                            newCost+=2;
                        }
                    }
                }
                int inner_distance = diff.GetDifference();
                newCost += inner_distance;
                // state->matrix[u_curr][v_curr] = newCost;
                ChangeCost(u_curr, v_curr, newCost, state->matrix, state->alpha, state->beta, true);
                // vertex.push_back(std::make_tuple(u_curr, v_curr, newCost));
            }
            int from_null = BranchEditDistanceFromNull(G2->GetBranch(v_curr));
            for (int v_nbr : G2->GetNeighbors(v_curr)) {
                if (state->inverse_mapping[v_nbr] != -1) {
                    from_null++;
                    // std::cout << "cur nbr " <<v_curr << " " << v_nbr<<  " " << from_null << "\n";
              }
            }
            for(int i = G1->GetNumVertices(); i < G2->GetNumVertices();i++){
              // std::cout << "cur nbr " <<v_curr << " " << i <<  " " << from_null << "\n";
              state->matrix[i][v_curr] = from_null;
              ChangeCost(i, v_curr, from_null, state->matrix, state->alpha, state->beta, true);
            }
            // ComputeBeta(state, v_curr);
        }
    // return vertex;
    }

#ifdef DD
std::pair<int, int>DHLowerBound(DHState *state){
    int u = 0;
    int v = 0;
    state->matrix.resize(G2->GetNumVertices(), std::vector<int>(G2->GetNumVertices(), 0));
    state->hungarian_assignment.resize(G2->GetNumVertices(), -1);
    state->alpha.resize(G2->GetNumVertices());
    int lb = 0, ub = 0;
    Timer hg_timer;
    Timer bd_timer;

    if(state->depth == -1){
        bd_timer.Start();
        ComputeBranchDistanceMatrixInitial(state);
        bd_timer.Stop();
        branchdistance_time += bd_timer.GetTime();

        N = G2->GetNumVertices();
        InitializeSolve_(state->alpha, state->beta);
        
        hg_timer.Start();
        Solve_(state->matrix, state->alpha, state->beta);
        
        exit(0);
        hg_timer.Stop();
        hungarian_time += hg_timer.GetTime();
        state->hungarian_assignment = assignment;
        std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
        std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
        std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
        std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            if(state->mapping[i] != -1) continue;
            hungarian_mapping[i] = state->hungarian_assignment[i];
            hungarian_inverse_mapping[state->hungarian_assignment[i]] = i;
        }
        ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
        lb = state->cost + ((total_cost + 1) / 2);
    }
    else{
        u = matching_order[state->depth];
        v = state->mapping[u];
        N = G2->GetNumVertices();
        state->matrix = static_cast<DHState *>(state->parent)->matrix;
        assignment = static_cast<DHState *>(state->parent)->hungarian_assignment; //copy to local assignment
        state->alpha = static_cast<DHState *>(state->parent)->alpha;
        state->beta = static_cast<DHState *>(state->parent)->beta;
        inverse_assignment.resize(G2->GetNumVertices(), -1);
        for(int i = 0; i < assignment.size();i++){ inverse_assignment[assignment[i]] = i;}

        bd_timer.Start();
        Match_(state);
        ComputeBranchDistanceMatrixDynamic(state);
        bd_timer.Stop();
        branchdistance_time += bd_timer.GetTime();
        
        std::vector<int> remain_left, remain_right;
        for(int i = 0 ; i < G1->GetNumVertices(); i++){if(state->mapping[i] == -1){remain_left.emplace_back(i);}}
        for(int i = G1->GetNumVertices(); i < G2->GetNumVertices();i++){remain_left.emplace_back(i);}
        for(int i = 0 ; i < G2->GetNumVertices(); i++){if(state->inverse_mapping[i] == -1){remain_right.emplace_back(i);}}

        int remaining = G2->GetNumVertices() - (state->depth + 1);
        std::vector<std::vector<int>> remain_matrix(remaining, (std::vector<int>(remaining, 0)));
        std::vector<int> remain_left_idx(N, -1);
        std::vector<int> remain_right_idx(N, -1);

        for(int v_idx = 0 ; v_idx < remain_right.size(); v_idx++){
          for(int u_idx = 0 ; u_idx < remain_left.size(); u_idx++){
            int u_curr = remain_left[u_idx];
            int v_curr = remain_right[v_idx];
            remain_matrix[u_idx][v_idx] = state->matrix[u_curr][v_curr];
          }
        }

        Hungarian hungarian(remain_matrix);
        hungarian.CopyAssignmentAlphaBeta(assignment, inverse_assignment, state->alpha, state->beta, remain_left, remain_right);

        // exit(0);


        hg_timer.Start();
        // SolvePartial_p(state, remain_left, remain_right);
        hungarian.SolvePartial();
        hg_timer.Stop();
        hungarian_time += hg_timer.GetTime();
        hungarian.CopyAssignmentAlphaBetaReverse(assignment, inverse_assignment, state->alpha, state->beta, remain_left, remain_right);


        state->hungarian_assignment = assignment;
        std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
        std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
        std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
        std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
        for(int i = 0 ; i < G1->GetNumVertices();i++){
            if(state->mapping[i] == -1){
                int u_ = i;
                int v_ = assignment[i];
                hungarian_mapping[u_] = v_;
                hungarian_inverse_mapping[v_] = u_;
            }
        }
        ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
        lb = state->cost + ((total_cost + 1) / 2);
    }
    return{lb, ub};
}

#endif
#ifdef CC
std::pair<int, int>DHLowerBound(DHState *state){
    int u = 0;
    int v = 0;
    state->matrix.resize(G2->GetNumVertices(), std::vector<int>(G2->GetNumVertices(), 0));
    state->hungarian_assignment.resize(G2->GetNumVertices(), -1);
    state->alpha.resize(G2->GetNumVertices());
    int lb = 0, ub = 0;
    
    Timer hg_timer;
    Timer bd_timer;

    if(state->depth == -1){
        bd_timer.Start();
        ComputeBranchDistanceMatrixInitial(state);
        bd_timer.Stop();
        branchdistance_time += bd_timer.GetTime();

        N = G2->GetNumVertices();
        InitializeSolve_(state->alpha, state->beta);
        
        hg_timer.Start();
        Solve_(state->matrix, state->alpha, state->beta);
        hg_timer.Stop();
        hungarian_time += hg_timer.GetTime();
        
        state->hungarian_assignment = assignment;
        std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
        std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
        std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
        std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            if(state->mapping[i] != -1) continue;
            hungarian_mapping[i] = state->hungarian_assignment[i];
            hungarian_inverse_mapping[state->hungarian_assignment[i]] = i;
        }
        ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
        lb = state->cost + ((total_cost + 1) / 2);
    }
    else{
        u = matching_order[state->depth];
        v = state->mapping[u];
        N = G2->GetNumVertices();
      
        state->matrix.resize(N, (std::vector<int>(N, 0)));    
        state->matrix = static_cast<DHState *>(state->parent)->matrix;
        assignment = static_cast<DHState *>(state->parent)->hungarian_assignment; //copy to local assignment
        state->alpha = static_cast<DHState *>(state->parent)->alpha;
        state->beta = static_cast<DHState *>(state->parent)->beta;
        inverse_assignment.resize(G2->GetNumVertices(), -1);

        for(int i = 0; i < assignment.size();i++){ inverse_assignment[assignment[i]] = i;}
        bd_timer.Start();
        Match_(state);

        ComputeBranchDistanceMatrixDynamic(state);
        bd_timer.Stop();
        branchdistance_time += bd_timer.GetTime();
        

        std::vector<int> remain_left, remain_right;
        for(int i = 0 ; i < G1->GetNumVertices(); i++){if(state->mapping[i] == -1){remain_left.emplace_back(i);}}
        for(int i = G1->GetNumVertices(); i < G2->GetNumVertices();i++){remain_left.emplace_back(i);}
        for(int i = 0 ; i < G2->GetNumVertices(); i++){if(state->inverse_mapping[i] == -1){remain_right.emplace_back(i);}}

        // std::cout << state->alpha << "\n";
        // std::cout << state->beta << "\n";

        hg_timer.Start();
        SolvePartial_p(state, remain_left, remain_right);
        hg_timer.Stop();
        hungarian_time += hg_timer.GetTime();
        state->hungarian_assignment = assignment;
        std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
        std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
        std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
        std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
        for(int i = 0 ; i < G1->GetNumVertices();i++){
            if(state->mapping[i] == -1){
                int u_ = i;
                int v_ = assignment[i];
                hungarian_mapping[u_] = v_;
                hungarian_inverse_mapping[v_] = u_;
            }
        }
        ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
        lb = state->cost + ((total_cost + 1) / 2);
    }
    return{lb, ub};
}

#endif


#ifdef AA
  std::pair<int, int> DHLowerBound(DHState* state) {
    // std::cerr << __LINE__ << "\n";
    // DifferenceVector diff;
    // diff.init(20);
    int remaining = G2->GetNumVertices() - (state->depth + 1);
    std::vector<int> rem_left, rem_right;
    std::vector<std::vector<int>> branch_distance_matrix(
        remaining, std::vector<int>(remaining, 0));

    Timer bd_timer;
    bd_timer.Start();
    ComputeBranchDistanceMatrix(state, branch_distance_matrix, rem_left,
                                rem_right);
    
    bd_timer.Stop();
    branchdistance_time += bd_timer.GetTime();
    Hungarian hungarian(branch_distance_matrix);

    Timer hg_timer;
    hg_timer.Start();
    hungarian_vertex_num += remaining;
    hungarian.Solve();
    hg_timer.Stop();
    hungarian_time += hg_timer.GetTime();


    if (DEBUG) hungarian.Print();

    auto& assignment = hungarian.GetAssignment();
    std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
    std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);

    std::memcpy(hungarian_mapping.data(), state->mapping,
                sizeof(int) * G1->GetNumVertices());
    std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping,
                sizeof(int) * G2->GetNumVertices());
    for (int i = 0; i < rem_left.size(); i++) {
      int u = rem_left[i];
      int v = rem_right[assignment[i]];
      hungarian_mapping[u] = v;
      hungarian_inverse_mapping[v] = u;

    }
    int u__ = matching_order[state->depth];
    int v__ = state->mapping[u__];

    int ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
    int lb = state->cost + ((hungarian.GetTotalCost() + 1) / 2);
    int u = matching_order[state->depth];
    int v = state->mapping[u];

    // if(dbg && state->depth == cnt - 1){
    //         std::cout << "depth : " << state->depth << "\n";
    //         std::cout << "u, v : "<<u__ << " " << v__ << "\n";
    //         std::cout <<"remain u : ";
    //         for(int i = 0 ; i < rem_left.size(); i++){
    //             std::cout << rem_left[i] << " ";
    //         }

            // PrintCostMatrix(branch_distance_matrix);
    //         std::cout << "\n";
    //         PrintStateMapping(state);
            // std::cout << "hg mapping : " <<  hungarian_mapping << "\n";
                // for(int i = 0 ;i < G1->GetNumVertices();i++){   std::cout << state->inverse_mapping[i] <<" ";}

    //         std::cout << "lb : " << lb << " ub : " << ub << "\n";
    //         std::cout << "hg totalcost : " << hungarian.GetTotalCost() <<"\n----------------------\n";
        // std::cout << "\n" << hungarian.GetTotalCost() << "\n";
    //     // std::cout << "\n";
    //     cnt++;
    // }
        // std::cout << "u v lb ub : " << u << " " << v << " " <<lb << " " << ub << "\n";
    // PrintCostMatrix(branch_distance_matrix);
    
      // std::cout << "u : " << u << " v : "<< v  << " cost : " << hungarian.GetTotalCost() <<  " ub : " << ub << " lb : " << lb  << "\n";
      // std::cout << hungarian_mapping << "\n--------------------------\n";
        // std::cout << "u v " << u << " " << v << "\n";
        // // std::cout << "lb, ub " << lb << " " << ub <<"\n";
        // std::cout << hungarian.GetTotalCost() << "\n";
        // std::cout << hungarian_mapping << "\n";
    return {lb, ub};
  }//DHLowerbound original version
#endif



double Gethgtime() const {return hungarian_time; }
double Getbdtime()const {return branchdistance_time; }
int64_t GetVertNum()const {return hungarian_vertex_num;}

int TotalCost(DHState *state){
    int total_cost = 0;
    for(int i = 0 ; i < G1->GetNumVertices();i++){
        if(state->mapping[i] == -1){
            total_cost += state->matrix[i][state->hungarian_assignment[i]];
        }
    }
    return total_cost;
}
void PrintCostMatrix(std::vector<std::vector<int>>& matrix){
    for(int i = 0 ; i < matrix.size();i++){
        for(int j = 0 ; j < matrix.size();j++){
            std::cout << matrix[i][j] << " "; 
        }
        std::cout << "\n";
    }
    std::cout << "-----------------------\n";
}
void PrintMatching(std::vector<int> &matching){
    for(int i = 0 ; i < G2->GetNumVertices();i++){
        std::cout << matching[i] << " ";
    }
    std::cout << "\n";
}
void PrintStateMapping(DHState * state){
    std::cout << "state mapping : ";
    for(int i = 0 ; i < G1->GetNumVertices();i++){
        std::cout << state->mapping[i] << " ";
    }
    std::cout << "\n";
}
void PrintCurrentMatching(DHState* state){
    int u = matching_order[state->depth];
    int v = state->mapping[u];
    std::cout << "u : " << u  << " v : " << v <<"\n";

}
};  // namespace GraphLib::GraphSimilarity
}