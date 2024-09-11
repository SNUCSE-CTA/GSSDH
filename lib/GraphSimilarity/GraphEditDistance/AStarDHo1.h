#pragma once
#include "GraphSimilarity/EditDistance.h"
#include "Base/Timer.h"

namespace GraphLib::GraphSimilarity{
struct DHoState : State{
    int next_mapping_order = -1;
    int ub = 1e9;
    std::vector<std::tuple<int, int, int>> dv; //save u, v, cost
    std::vector<int> hungarian_assignment;
    std::vector<int> alpha, beta;
};

class AStarDHo : public GraphEditDistanceSolver{
    std::priority_queue<DHoState*, std::vector<DHoState*>, StateComparator> queue;
    const int INF = 1e9;
    std::vector<std::vector<int>> global_matrix;
    std::vector<int> assignment, inverse_assignment;
    std::vector<int> left_visited, right_visited;
    std::vector<int> alpha, beta;
    int total_cost, N, theta;

    std::vector<int> restore_inverse_assignment; //size = G2->GetNumVertices()
    std::vector<int> restore_u_idx;
    std::vector<int> restore_v_idx;
    int cnt = 0;
    public: 

    void Initialize(DHoState* state){
        total_cost = INF;
        assignment = std::vector<int>(N, -1);
        inverse_assignment = std::vector<int>(N, -1);
        left_visited = std::vector<int>(N, 0);
        right_visited = std::vector<int>(N, 0);
        state->alpha = std::vector<int>(N, 0);
        state->beta = std::vector<int>(N, INF);
        theta = 0;
        for(int i = 0; i < N; i++){
            for(int j = 0 ; j < N; j++){
                state->beta[j] = std::min(state->beta[j], global_matrix[i][j]);
            }
        }
        alpha = state->alpha;
        beta = state->beta;
    }
    bool FindAungmentingPath(int i, std::vector<std::vector<int>> cost_matrix){
        left_visited[i] = true;
        for(int j = 0 ; j < N; j++){
            if(!right_visited[j] &&(alpha[i] + beta[j] == cost_matrix[i][j])){
                right_visited[j] = true;
                if(inverse_assignment[j] == -1 || FindAungmentingPath(inverse_assignment[j], cost_matrix)){
                    inverse_assignment[j] = i;
                    assignment[i] = j;
                    return true;
                }
            }
        }
        return false;
    }
    void RecalculatePotential(std::vector<std::vector<int>>& cost_matrix) {
        theta = INF;
        for(int i = 0 ; i < N; i++){
            if(left_visited[i]){
                for(int j = 0 ; j < N; j++){
                    if(!right_visited[j]){
                        theta = std::min(theta, cost_matrix[i][j] - alpha[i] - beta[j]);
                    }
                }
            }
        }
        for(int i = 0 ; i < N; i++){
            if(left_visited[i]){
                alpha[i] += theta;
            }
        }
        for(int j = 0 ; j < N ; j++){
            if(right_visited[j]){
                beta[j] -= theta;
            }
        }
    }
    void Solve(std::vector<std::vector<int>> cost_matrix){ //current matrix, current alpha beta, current assignment
        for(int i = 0 ; i < N; i++){
            if(assignment[i] == -1){
                if(cnt == 1) {
                    std::cout << "look " << i << "\n";
                }
                while(true){
                    std::fill(left_visited.begin(), left_visited.end(), 0);
                    std::fill(right_visited.begin(), right_visited.end(), 0);
                    if(FindAungmentingPath(i, cost_matrix)){
                        break;
                    }
                    RecalculatePotential(cost_matrix);
                }
            }
        }
        total_cost = 0;
        for(int i = 0 ; i < N; i++){
            total_cost += cost_matrix[i][assignment[i]];
            // std::cout << cost_matrix[i][assignment[i]] << " ";
        }
        // std::cout<< "\n";
    }

    void ComputeBranchDistanceMatrixInitial(DHoState* state){
    DifferenceVector diff;
    diff.init(20);
    global_matrix.resize(G2->GetNumVertices(), (std::vector<int>(G2->GetNumVertices(), 0)));
    int colSize = G2->GetNumVertices(); //right
    int rowSize = G1->GetNumVertices(); //left
    for(int v = 0 ; v < colSize; v++){
        auto &v_nbrs = G2->GetNeighbors(v);
        int u = 0;
        for(u = 0 ; u < rowSize; u++){
            auto &u_nbrs = G1->GetNeighbors(u);
            diff.reset();
            if(G1->GetVertexLabel(u) != G2->GetVertexLabel(v)){
                global_matrix[u][v] += 2;
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
                        global_matrix[u][v] += 2;
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
                        global_matrix[u][v] += 2;
                    }
                }
            }
            int inner_distance = diff.GetDifference();
            global_matrix[u][v] += inner_distance;
        }
        int from_null = BranchEditDistanceFromNull(G2->GetBranch(v));
        for(int v_nbr : G2->GetNeighbors(v)){
            if(state->inverse_mapping[v_nbr] != -1){
                from_null++;
            }
        }
        for(; u < colSize; u++){
            global_matrix[u][v] = from_null;
        }
    }
  }

    void ComputeBranchDistanceMatrixDynamic(DHoState *state){ //return u, v, cost
    // void ComputeBranchDistanceMatrixDynamic(DHState *state){ 
        DifferenceVector diff;
        diff.init(20);
        // std::vector<std::tuple<int,int,int>> edge_set;
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
                state->dv.push_back(std::make_tuple(u_curr, v_curr, newCost));
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
                state->dv.push_back(std::make_tuple(u_curr, v_curr, newCost));
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
            //   state->matrix[i][v_curr] = from_null;
                state->dv.push_back(std::make_tuple(i, v_curr, from_null));
            }
            // ComputeBeta(state, v_curr);
        }
        // return edge_set;
    }

    void MakeMatrix(DHoState* state,std::vector<std::vector<int>>& local_matrix_full, std::vector<std::vector<int>>& local_matrix, std::vector<int>& rem_left, std::vector<int>& rem_right){
        for(int v_idx = 0 ; v_idx < rem_right.size();v_idx ++){
            int v = rem_right[v_idx];
            int u_idx = 0;
            for(u_idx = 0 ; u_idx < rem_left.size();u_idx++){
                int u = rem_left[u_idx];
                local_matrix[u_idx][v_idx] = local_matrix_full[u][v];
            }
        }
    }

    void ChangeCost(DHoState* state, int u, int v, int newCost, std::vector<std::vector<int>>& local_matrix_full){
        int oldCost = local_matrix_full[u][v];
        local_matrix_full[u][v] = newCost;
        if(newCost > oldCost && state->hungarian_assignment[u] == v){
            state->hungarian_assignment[u] = -1;
            restore_inverse_assignment[v] = -1;
        }
    }

    void ChangeMatrix(DHoState* state, std::vector<std::vector<int>>& local_matrix_full){
        for(int i = 0 ; i < state->dv.size();i++){
            int u = std::get<0>(state->dv[i]);
            int v = std::get<1>(state->dv[i]);
            int cost = std::get<2>(state->dv[i]);
            local_matrix_full[u][v] = cost;
        }
    }
    void ChangeMatrixCurrent(DHoState* state, std::vector<std::vector<int>>& local_matrix_full){
        for(int i = 0; i <state->dv.size();i++){
            int u = std::get<0>(state->dv[i]);
            int v = std::get<1>(state->dv[i]);
            int cost = std::get<2>(state->dv[i]);
            ChangeCost(state, u, v, cost, local_matrix_full);
        }
    }
    
    void RestoreInverseAssignment(DHoState* state){
        restore_inverse_assignment.resize(G2->GetNumVertices(), -1);
        for(int i = 0; i < state->hungarian_assignment.size(); i++){
            int u = i;
            int v = state->hungarian_assignment[u];
            if(v != -1){
                restore_inverse_assignment[v] = u;
            }
        }
    }

    void Match(DHoState * state, std::vector<std::vector<int>>& local_matrix_full){
        int u = matching_order[state->depth];
        int v = state->mapping[u];
        for(int j = 0 ; j <G2->GetNumVertices();j++){
            if(j == v) continue;
            ChangeCost(state, u, j, INF, local_matrix_full);
        }
        for(int i = 0 ; i < G2->GetNumVertices();i++){
            if(i == u) continue;
            ChangeCost(state, i, v, INF, local_matrix_full);
        }
        restore_inverse_assignment[v] = u;
        state->hungarian_assignment[u] = v;
    }

    void CopyAssignment(DHoState* state, std::vector<int>&rem_left, std::vector<int>& rem_right){
        restore_u_idx.resize(G2->GetNumVertices(), -1);
        restore_v_idx.resize(G2->GetNumVertices(), -1);
        for(int i = 0 ; i < G1->GetNumVertices();i++){
            if(state->mapping[i] == -1){
                restore_u_idx[i] = rem_left.size();
                rem_left.emplace_back(i);
            }
        }
        for(int i = G1->GetNumVertices(); i < G2->GetNumVertices();i++){
            restore_u_idx[i] = rem_left.size();
            rem_left.emplace_back(i);
        }
        for(int i = 0 ; i < G2->GetNumVertices();i++){
            if(state->inverse_mapping[i] == -1){
                restore_v_idx[i] = rem_right.size();//restore_v_idx[v] = v_idx
                rem_right.emplace_back(i);
            }
        }
        assignment.resize(rem_left.size());
        fill(assignment.begin(), assignment.end(), -1);
        inverse_assignment.resize(rem_right.size());
        fill(inverse_assignment.begin(), inverse_assignment.end(), -1);
        
        for(int u_idx = 0; u_idx < rem_left.size(); u_idx++){
            int u = rem_left[u_idx];
            int v = state->hungarian_assignment[u];
            alpha[u_idx] = state->alpha[u];
            if(v != -1){
                int v_idx = restore_v_idx[v];
                assignment[u_idx] = v_idx;
            }
        }
        for(int v_idx = 0 ; v_idx < rem_right.size();v_idx++){
            int v = rem_right[v_idx];
            int u = restore_inverse_assignment[v];
            beta[v_idx] = state->beta[v];
            if(u != -1){
                int u_idx = restore_u_idx[u];
                inverse_assignment[v_idx] = u_idx;
            }
        }
    }

    void RestoreAssignment(DHoState* state, std::vector<int>& rem_left, std::vector<int>& rem_right){
        // std::cout << state->alpha << "\n";
        for(int u_idx = 0 ; u_idx < rem_left.size(); u_idx++){
            int u = rem_left[u_idx];
            // std::cout <<u << " " << u_idx << " " << rem_right[assignment[u_idx]] << "\n";
            state->alpha[u] = alpha[u_idx];
            state->hungarian_assignment[u] = rem_right[assignment[u_idx]];
        }
        for(int v_idx = 0 ; v_idx < rem_right.size(); v_idx++){
            int v = rem_right[v_idx];
            state->beta[v] = beta[v_idx];
        }
    }

    std::pair<int, int>DHoLowerBound(DHoState *state){
        int ub, lb;
        N = G2->GetNumVertices() - (state->depth + 1);
        if(state->depth == -1){
            ComputeBranchDistanceMatrixInitial(state);
            Initialize(state);
            Solve(global_matrix);
            state->hungarian_assignment = assignment;
            std::cout << assignment << "\n" << inverse_assignment <<"\n";
            
            state->alpha = alpha;
            state->beta = beta;
            std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
            std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
            std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
            std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
            for(int i = 0; i < G1->GetNumVertices(); i++){
                if(state->mapping[i] != -1) continue;
                hungarian_mapping[i] = state->hungarian_assignment[i];
                hungarian_inverse_mapping[state->hungarian_assignment[i]] = i;
            }
            ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
            lb = state->cost + ((total_cost + 1) / 2);
            // std::cout << assignment.size() << " " << assignment << "\n";
            // std::cout << inverse_assignment << "\n";
        }
        else{
            state->hungarian_assignment = static_cast<DHoState*>(state->parent)->hungarian_assignment;
            state->alpha = static_cast<DHoState*>(state->parent)->alpha;
            state->beta = static_cast<DHoState*>(state->parent)->beta;
            std::vector<int> rem_left, rem_right;
            std::vector<std::vector<int>> local_matrix_full = global_matrix;
            std::vector<std::vector<int>> local_matrix(N, (std::vector<int>(N, 0)));
            DHoState *p_state = static_cast<DHoState*>(state->parent);
            std::stack<DHoState*> stack;
            while(p_state != NULL){
                stack.push(p_state);
                p_state = static_cast<DHoState*>(p_state->parent);
            }
            while(!stack.empty()){
                DHoState * curr_state = stack.top();
                if(curr_state->depth != -1){
                    ChangeMatrix(curr_state, local_matrix_full);
                }
                stack.pop();
            }

            RestoreInverseAssignment(state);
            ComputeBranchDistanceMatrixDynamic(state);
            ChangeMatrixCurrent(state, local_matrix_full);//current matrix
            Match(state, local_matrix_full);
            CopyAssignment(state, rem_left, rem_right);
            MakeMatrix(state, local_matrix_full, local_matrix, rem_left, rem_right);
            // for(int i = 0 ; i < N;i++){
            //     std::cout << i << " " << inverse_assignment[i] << "\n";
            // }
            // std::cout << state->hungarian_assignment << "\n" << restore_inverse_assignment <<"\n";
            Solve(local_matrix);
            RestoreAssignment(state, rem_left, rem_right);
            // RestoreInverseAssignment(state);

            std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
            std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
            std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
            std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
            for(int u_idx = 0; u_idx < rem_left.size() - (G2->GetNumVertices() - G1->GetNumVertices()) ;u_idx++){
                int u = rem_left[u_idx];
                int v = rem_right[assignment[u_idx]];
                hungarian_mapping[u] = v;
                hungarian_inverse_mapping[v] = u;
            }
            
            ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
            lb = state->cost + ((total_cost + 1) / 2);
            std::cout << hungarian_mapping << "\n" << hungarian_inverse_mapping << "\n";
            exit(0);
        }
        return {lb, ub};
    }

    void ExtendState(DHoState* state){
        if(state->cost >= current_best) return;
        int depth = state->depth + 1;
        int u = matching_order[depth];
        for(int v = 0; v < G2->GetNumVertices(); v++){
            if(state->inverse_mapping[v] != -1) continue;
            DHoState* child_state = new DHoState(state);
            child_state->cost = GetChildEditCost(state, u, v);
            child_state->mapping[u] = v;
            child_state->inverse_mapping[v] = u;
            child_state->depth = depth;
            auto[lb, ub] = DHoLowerBound(child_state);
            child_state->lower_bound = lb;
            child_state->ub = ub;
            if(lb == ub) continue;
            if(depth == G1->GetNumVertices() - 1){
                current_best = std::min(current_best, child_state->lower_bound);
                memcpy(&current_best_mapping[0], child_state->mapping, sizeof(int) * NumG1Vertices);
                continue;
            }
            if(child_state->lower_bound >= current_best) continue;
            if(threshold > 0){
                if(child_state->lower_bound > threshold) continue;
            }
            queue.push(child_state);
        }
    }
    
    int GED(){
        PrepareGED();
        DHoState* initial_state = new DHoState(NULL);
        initial_state->cost = 0;
        initial_state->depth = -1;
        auto [lb, ub] = DHoLowerBound(initial_state);
        initial_state->lower_bound = lb;
        queue.push(initial_state);
        int64_t max_qsize = 1;
        while(!queue.empty()){
            DHoState* current_state = queue.top();
            num_nodes++;
            queue.pop();
            if(current_state->lower_bound >= current_best){
                queue = std::priority_queue<DHoState*, std::vector<DHoState*>, StateComparator>();
                break;
            }
            if(threshold >= 0){
                if(current_best < threshold){
                    queue = std::priority_queue<DHoState*, std::vector<DHoState*>, StateComparator>();
                    break;
                }
                if(current_state->lower_bound > threshold){
                    current_best = -1;
                    queue = std::priority_queue<DHoState*, std::vector<DHoState*>, StateComparator>();
                    break;
                }
            }
            ExtendState(current_state);
            max_qsize = std::max(max_qsize, (int64_t)queue.size());
        }
        if(threshold >= 0 and current_best > threshold){
            current_best = -1;
        }
        log.AddResult("MaxQueueSize", max_qsize, RESULT_INT64);
        log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
        log.AddResult("EditDistance", current_best, RESULT_INT);
        return current_best;
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


};


}//GraphLib::GraphSimilarity