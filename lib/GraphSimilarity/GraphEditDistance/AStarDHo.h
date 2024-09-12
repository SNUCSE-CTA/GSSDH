#pragma once
#include "GraphSimilarity/EditDistance.h"
#include "Base/Timer.h"

namespace GraphLib::GraphSimilarity{
struct DHoState : State{
    int next_mapping_order = -1;
    int ub = 1e9;
    std::vector<std::vector<int>> matrix;
    std::vector<int> hungarian_assignment;
    std::vector<int> hungarian_inverse_assignment;
    std::vector<int> alpha, beta;
};

class AStarDHo : public GraphEditDistanceSolver{
    std::priority_queue<DHoState*, std::vector<DHoState*>, StateComparator> queue;
    const int INF = 1e9;

    std::vector<int> assignment, inverse_assignment;
    std::vector<int> left_visited, right_visited;
    std::vector<int> alpha, beta;
    int total_cost, N, theta;
    int cnt = 0;

    std::vector<int> __left;
    std::vector<int> __right1, __right2;

	int acc = 0;

    int64_t functioncall = 0;
    
    double hgtime = 0.0;
    double bdtime = 0.0;


    std::vector<int> u_idxs, v_idxs;//u_idxs[u] = u_idx, v_idxs[v] = v_idx;

    public: 
    void Initialize(DHoState* state){
        total_cost = INF;
        assignment = std::vector<int>(N, -1);
        inverse_assignment = std::vector<int>(N, -1);
        left_visited = std::vector<int>(N, 0);
        right_visited = std::vector<int>(N, 0);
        u_idxs = std::vector<int>(N, -1);
        v_idxs = std::vector<int>(N, -1);
        state->alpha = std::vector<int>(N, 0);
        state->beta = std::vector<int>(N, INF);
        theta = 0;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                state->beta[j] = std::min(state->beta[j], state->matrix[i][j]);
            }
        }
        alpha = state->alpha;
        beta = state->beta;

		// for (int i = 0; i < N; ++i) {
		// 	for (int j = 0; j < N; ++j) {
		// 		std::cerr << state->matrix[i][j] << " ";
		// 	}
		// 	std::cerr << std::endl;
		// }

		// std::cerr << "a: ";
		// for (int i = 0; i < N; ++i) {
		// 	std::cerr << alpha[i] << " ";
		// }
		// std::cerr << std::endl;
		// std::cerr << "b: ";
		// for (int i = 0; i < N; ++i) {
		// 	std::cerr << beta[i] << " ";
		// }
		// std::cerr << std::endl;
		// std::cerr << std::endl;
    }

    bool FindAugmentingPath(int i, const std::vector<std::vector<int>>& cost_matrix) {
        __left.push_back(i);

        for (int j = 0; j < N; j++) {
            if (__builtin_expect(right_visited[j], true)) {
				continue;
			}

            if (__builtin_expect(alpha[i] + beta[j] != cost_matrix[i][j], true)) continue;
            right_visited[j] = true;
            __right1.push_back(j);

            if (inverse_assignment[j] == -1 || FindAugmentingPath(inverse_assignment[j], cost_matrix)) {
                inverse_assignment[j] = i;
                assignment[i] = j;
                return true;
            }
        }
		return false;
  	}

    void RecalculatePotential(const std::vector<std::vector<int>>& cost_matrix) {
        theta = INF;
        for (int i: __left) {
                for (int j = 0 ; j < N; j++) {
                    if (!right_visited[j]) {
                        theta = std::min(theta, cost_matrix[i][j] - alpha[i] - beta[j]);
                    }
                }
        }
        for (int i: __left) {
                alpha[i] += theta;
				acc += theta;
		}
        for (int j: __right1) {
                beta[j] -= theta;
				acc -= theta;
		}
    }

    void Solve(const std::vector<std::vector<int>>& cost_matrix, DHoState *state, std::vector<int>&rem_left, std::vector<int>&rem_right){

		acc = 0;
		for (int i = 0; i < N; ++i) {
			acc += alpha[i];
		}
		for (int j = 0; j < N; ++j) {
			acc += beta[j];
		}

		/* !! alpha-descending order */
		std::vector<std::pair<int, int>> mappingOrder;
		mappingOrder.reserve(N);
		for (int i = 0; i < N; ++i) {
			if (assignment[i] == -1) {
				mappingOrder.push_back(std::pair<int, int>{ alpha[i], i });
			}
		}
		sort(mappingOrder.rbegin(), mappingOrder.rend());
		for (const auto& [_, i]: mappingOrder) {
		/* simple order */
        // for (int i = 0; i < N; i++) {
        //     if (assignment[i] != -1) {
		// 		continue;
		// 	}
            functioncall++;

			while (true) {
            	int lb = state->cost + ((acc + 1) / 2);
				// if (false) {
				if (lb > threshold) {
					total_cost = 500;
					for (int i = 0; i < N; ++i) {
						assignment[i] = inverse_assignment[i] = i;
					}
					return;
				}

				std::fill(right_visited.begin(), right_visited.end(), 0);

				__left.clear();
				__right1.clear();

				const bool failed = FindAugmentingPath(i, cost_matrix);
				if (failed) {
					break;
				}

				RecalculatePotential(cost_matrix);
				// std::cerr << "alpha: ";
				// for (int i = 0; i < N; ++i) {
				// 	if (assignment[i] == -1) {
				// 		std::cerr << alpha[i] << " ";
				// 	}
				// }
				// std::cerr << std::endl;
				// std::cerr << " beta: ";
				// for (int j = 0; j < N; ++j) {
				// 	if ( inverse_assignment[j] != -1 ) {
				// 		std::cerr << beta[j] << " ";
				// 	}
				// }
				// std::cerr << std::endl;

			}
        }
        total_cost = 0;
        for (int i = 0; i < N; i++){
            total_cost += cost_matrix[i][assignment[i]];
            if (state->depth != -1) {
                int u = rem_left[i];
                int v = rem_right[i];
                int u_idx = u_idxs[u];
                int v_idx = v_idxs[v];
                state->alpha[u] = alpha[u_idx];
                state->beta[v] = beta[v_idx];
                state->hungarian_assignment[u] = rem_right[assignment[u_idx]];
                state->hungarian_inverse_assignment[v] = rem_left[inverse_assignment[v_idx]];
            }
        }
    }

    void ChangeCost(int i, int j, int newCost, DHoState* state){
        int oldCost = state->matrix[i][j];
        state->matrix[i][j] = newCost;
        if(newCost > oldCost && state->hungarian_assignment[i] == j){  
            state->hungarian_assignment[i] = -1;
            state->hungarian_inverse_assignment[j] = -1;
        }
    }

    void ComputeBranchDistanceMatrixInitial(DHoState *state){
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
    void ComputeBranchDistanceMatrixDynamic(DHoState *state){ 
        DifferenceVector diff;
        diff.init(20);
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
                ChangeCost(u_curr, v_curr, newCost, state);
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
                ChangeCost(u_curr, v_curr, newCost, state);
            }
            int from_null = BranchEditDistanceFromNull(G2->GetBranch(v_curr));
            for (int v_nbr : G2->GetNeighbors(v_curr)) {
                if (state->inverse_mapping[v_nbr] != -1) {
                    from_null++;
              }
            }
            for(int i = G1->GetNumVertices(); i < G2->GetNumVertices();i++){
            //   state->matrix[i][v_curr] = from_null;
              ChangeCost(i, v_curr, from_null, state);
            }
        }
    }
    void Match(DHoState* state){
        int u = matching_order[state->depth];
        int v = state->mapping[u];
        for(int j = 0 ; j < N; j++){
            if(j == v) continue;
            ChangeCost(u, j, INF, state);
        }
        for(int i = 0 ; i < N; i++){
            if(i == u) continue;
            ChangeCost(i, v, INF, state);
        }
        state->hungarian_inverse_assignment[v] = u;
        state->hungarian_assignment[u] = v;
    }
    void ComputeReducedMatrix(DHoState *state, std::vector<std::vector<int>> &local_matrix, std::vector<int>& rem_left, std::vector<int>& rem_right, int &remaining){
        left_visited.resize(remaining);
        right_visited.resize(remaining);
        alpha.resize(remaining);
        beta.resize(remaining);
        assignment.resize(remaining);
        inverse_assignment.resize(remaining);
        fill(alpha.begin(), alpha.end(), -1);
        fill(beta.begin(), beta.end(), -1);
        fill(assignment.begin(), assignment.end(), -1);
        fill(inverse_assignment.begin(), inverse_assignment.end(), -1);
        fill(u_idxs.begin(), u_idxs.end(), -1);
        fill(v_idxs.begin(), v_idxs.end(), -1);
        for(int u = 0; u < G1->GetNumVertices(); u++){
            if(state->mapping[u] == -1){
                u_idxs[u] = rem_left.size();
                rem_left.emplace_back(u);
                alpha[u_idxs[u]] = state->alpha[u];
            }
        }
        for(int u = G1->GetNumVertices(); u < G2->GetNumVertices();u++){
            u_idxs[u] = rem_left.size();
            rem_left.emplace_back(u);
            alpha[u_idxs[u]] = state->alpha[u];
        }
        for(int v = 0; v < G2->GetNumVertices(); v++){
            if(state->inverse_mapping[v] == -1){
                v_idxs[v] = rem_right.size();
                rem_right.emplace_back(v);
                beta[v_idxs[v]] = state->beta[v];;
            }
        }
        for(int i = 0 ; i < remaining; i++){
            if(state->hungarian_assignment[rem_left[i]] != -1){
                assignment[i] = v_idxs[state->hungarian_assignment[rem_left[i]]];//rem_left[i] == u, state->hungarian_assignment[rem_left[i]] = v,
            }
            if(state->hungarian_inverse_assignment[rem_right[i]] != -1){
                inverse_assignment[i] = u_idxs[state->hungarian_inverse_assignment[rem_right[i]]];
            }
        }
        // for(int i = 0; i < remaining; i++){
        //     int v = rem_right[i];
        //     int u = state->hungarian_inverse_assignment[v];
        //     if(u != -1){
        //         int u_idx = u_idxs[u];
        //         inverse_assignment[i] = u_idx; 
        //     }
        // }
        for(int v_idx = 0 ; v_idx < rem_right.size(); v_idx++){
            int v = rem_right[v_idx];
            for(int u_idx = 0;u_idx < rem_left.size();u_idx++){
                int u = rem_left[u_idx];
                local_matrix[u_idx][v_idx] = state->matrix[u][v];
            }
        }
                 // std::cerr << alpha << "\n" << beta << "\n" << assignment << "\n" << inverse_assignment << "\n";
    }
    void LocalToState(DHoState *state, std::vector<int>& rem_left, std::vector<int>& rem_right, int &remaining){
        for(int i = 0 ; i < remaining; i++){
            int u = rem_left[i];
            int v = rem_right[i];
            int u_idx = u_idxs[u];
            int v_idx = v_idxs[v];
            state->alpha[u] = alpha[u_idx];
            state->beta[v] = beta[v_idx];
            state->hungarian_assignment[u] = rem_right[assignment[u_idx]];
            state->hungarian_inverse_assignment[v] = rem_left[inverse_assignment[v_idx]];
        }
    }
    std::pair<int, int>DHoLowerBound(DHoState *state){
        int ub = 0, lb = 0;
        state->matrix.resize(G2->GetNumVertices(), std::vector<int>(G2->GetNumVertices(), 0));
        state->hungarian_assignment.resize(G2->GetNumVertices(), -1);
        state->alpha.resize(G2->GetNumVertices());
        state->beta.resize(G2->GetNumVertices());
        N = G2->GetNumVertices() - (state->depth + 1);
        // std::cout << N << "\n";
        Timer t;
        if(state->depth == -1){
            t.Start();
            ComputeBranchDistanceMatrixInitial(state);
            t.Stop();
            bdtime += t.GetTime();
            Initialize(state);
            std::vector<int> rem_left, rem_right;
            Timer b;
            b.Start();
            Solve(state->matrix, state, rem_left, rem_right);
            b.Stop();
            hgtime += b.GetTime();
            state->hungarian_assignment = assignment;
            state->hungarian_inverse_assignment = inverse_assignment;
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
        }
        else{
            state->matrix = static_cast<DHoState*>(state->parent)->matrix;
            state->alpha = static_cast<DHoState*>(state->parent)->alpha;
            state->beta = static_cast<DHoState*>(state->parent)->beta;
            state->hungarian_assignment = static_cast<DHoState*>(state->parent)->hungarian_assignment;
            state->hungarian_inverse_assignment = static_cast<DHoState*>(state->parent)->hungarian_inverse_assignment;
            t.Start();
            ComputeBranchDistanceMatrixDynamic(state);
            Match(state);
            t.Stop();
            bdtime += t.GetTime();
            int remaining = G2->GetNumVertices() - (state->depth + 1);
            std::vector<int> rem_left, rem_right;
            std::vector<std::vector<int>> local_matrix(remaining, std::vector<int>(remaining, 0));

			// std::cerr << state->alpha << std::endl;

            ComputeReducedMatrix(state, local_matrix, rem_left, rem_right, remaining);

			// std::cerr << alpha << std::endl;
			// std::cerr << std::endl;

            Timer a;
            a.Start();
            Solve(local_matrix, state, rem_left, rem_right);
            a.Stop();
            hgtime += a.GetTime();
            // LocalToState(state, rem_left, rem_right, remaining);

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
        PrepareGED(nullptr);
        DHoState* initial_state = new DHoState(NULL);
        hgtime = 0.0;
        bdtime = 0.0;
        cnt = 0;
        functioncall = 0;
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

double Gethgtime() const {return hgtime; }
double Getbdtime()const {return bdtime; }
int64_t GetCnt() const {return functioncall;}

};


}//GraphLib::GraphSimilarity
