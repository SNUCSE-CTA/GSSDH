#pragma once
// #include "Base/Hungarian.h"
#include "Base/DynamicHungarian.h"
#include "GraphSimilarity/EditDistance.h"
#include "Base/Timer.h"
// #include "Base/DynamicHungarian.h"

namespace GraphLib::GraphSimilarity {
struct DHState : State {
  int next_mapping_order = -1;
  int ub = 1e9;
  std::vector<std::vector<int>> matrix;
  std::vector<int> hungarian_assignment;
};


class AStarDH : public GraphEditDistanceSolver {
  const bool DEBUG = false;
  const bool dbg = false;
  int cnt = 0;
  std::priority_queue<DHState*, std::vector<DHState*>, StateComparator> queue;
  std::map<int, int> num_hungarian;
  double hungarian_time = 0.0, branchdistance_time = 0.0;
//   std::vector<vector<int>> BranchDistanceMatrix;//initial matrix
 public:
 
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

  std::vector<std::tuple<int,int,int>> ComputeBranchDistanceMatrixDynamic(DHState *state){ //return u, v, cost
        DifferenceVector diff;
        diff.init(20);
        std::vector<std::tuple<int,int,int>> vertex;
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
                vertex.push_back(std::make_tuple(u_curr, v_curr, newCost));
            }
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
                vertex.push_back(std::make_tuple(u_curr, v_curr, newCost));
            }
        }
    return vertex;
  }

#ifdef BB
std::pair<int, int> DHLowerBound(DHState* state){
        state->matrix.resize(G2->GetNumVertices(),
                            std::vector<int>(G2->GetNumVertices(), 0));
        if(state->depth == -1){
            ComputeBranchDistanceMatrixInitial(state);
            Hungarian hungarian(state->matrix);
            hungarian.Solve();
            // if(dbg){
            //     hungarian.Print();
            //     exit(0);
            // }
            auto& assignment = hungarian.GetAssignment();
            
            std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
            std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
            for(int i = 0 ; i < G1->GetNumVertices(); i++){
                state->hungarian_assignment.push_back(assignment[i]);
                hungarian_mapping[i] = assignment[i];
                hungarian_inverse_mapping[assignment[i]] = i;
            }
            int ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
            int lb = state->cost + ((hungarian.GetTotalCost() + 1) / 2);

            // if(dbg){hungarian.PrintAssignment();}
        return {lb, ub};
        }
        else{
            int u = matching_order[state->depth];
            int v = state->mapping[u];

            // std::cout<< u << " " << v << " " <<state->depth << "where is segfault" << "\n";
            std::vector<std::tuple<int,int,int>> vertices; //u, v, cost
            state->matrix = static_cast<DHState*>(state->parent)->matrix;
            Timer bd_timer;
            bd_timer.Start();
            vertices = ComputeBranchDistanceMatrixDynamic(state);
            bd_timer.Stop();
            branchdistance_time += bd_timer.GetTime();

            Hungarian hungarian(state->matrix);
            state->hungarian_assignment = static_cast<DHState*>(state->parent)->hungarian_assignment;
            hungarian.CopyParent(static_cast<DHState*>(state->parent)->hungarian_assignment);

                        

            std::vector<bool> visited(G1->GetNumVertices(), 0);
            std::queue<int> vqueue;
            int u_prime = hungarian.Match(u, v, state->matrix);
            if(u_prime != u && u_prime != -1 && visited[u_prime] == false){
                vqueue.push(u_prime);
                visited[u_prime] = true;
            }


            for(int i = 0 ; i < vertices.size();i++){
                std::tuple<int, int, int> u_v_cost = vertices[i];
                int u_ = std::get<0>(vertices[i]);
                int v_ = std::get<1>(vertices[i]);
                int cost_ = std::get<2>(vertices[i]);
                u_prime = hungarian.ChangeCost(u_, v_, cost_);
                if(u_prime != -1 && visited[u_prime] == false){
                    vqueue.push(u_prime);
                    visited[u_prime] = true;
                }
            }

            Timer hg_timer;
            hg_timer.Start();
            hungarian.InitializeVariables();
            while(!vqueue.empty()){
                u_prime = vqueue.front();
                vqueue.pop();
                hungarian.SolvePartial(u_prime);
                // std::cout<< u_prime << " " << v << " " <<state->depth << "where is segfault" << "\n";
            }

            hg_timer.Stop();
            hungarian_time += hg_timer.GetTime();




            hungarian.TotalCost(state->mapping);//calculate Total Cost after partial solve
            state->hungarian_assignment = hungarian.GetAssignment();

            std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
            std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
            for(int i = 0 ; i < G1->GetNumVertices(); i++){
                hungarian_mapping[i] = state->hungarian_assignment[i];
                hungarian_inverse_mapping[state->hungarian_assignment[i]] = i;
                // std::cout << hungarian_mapping[i] << " " << hungarian_inverse_mapping[state->hungarian_assignment[i]] <<"\n";
            }
            int ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
            int lb = state->cost + ((hungarian.GetTotalCost() + 1) / 2);

        return {lb, ub};
    }
}
#endif

#ifdef AA
  std::pair<int, int> DHLowerBound(DHState* state) {
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

    int ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
    int lb = state->cost + ((hungarian.GetTotalCost() + 1) / 2);



    int u = matching_order[state->depth];
    int v = state->mapping[u];
    return {lb, ub};
  }//DHLowerbound original version
#endif



  double Gethgtime() {return hungarian_time; }
  double Getbdtime()const {return branchdistance_time; }

void PrintCostMatrix(std::vector<std::vector<int>>& matrix){
    for(int i = 0 ; i < G2->GetNumVertices();i++){
        for(int j = 0 ; j < G2->GetNumVertices();j++){
            std::cout << matrix[i][j] << " "; 
        }
        std::cout << "\n";
    }
}
void PrintCurrentMatching(DHState* state){
    int u = matching_order[state->depth];
    int v = state->mapping[u];
    std::cout << "u : " << u  << " v : " << v <<"\n";
}

};  // namespace GraphLib::GraphSimilarity
}