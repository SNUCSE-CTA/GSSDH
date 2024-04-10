#pragma once
#include "Base.h"
#include "Timer.h"

/**
 * @brief Hungarian algorithm for Assignment problem (Minimum-weighted bipartite matching).
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
    Hungarian(std::vector<std::vector<int>> cost_matrix_) {
        this->N = cost_matrix_.size();
        InitializeSolver();
        cost_matrix = cost_matrix_;
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
            if (inverse_assignment[j] == -1 || FindAugmentingPath(inverse_assignment[j])) {
                inverse_assignment[j] = i;
                assignment[i] = j;
                return true;
            }
        }
        return false;
    }

    void RecalculatePotential() {
//        fprintf(stderr, "Invoked Hungarian::RecalculatePotential()\n");
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
//        fprintf(stderr, "Invoked Hungarian::Solve()\n");
//        Print();
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

    void Print() {
        fprintf(stderr, "Invoked Hungarian::Print() N = %d\n", N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(stderr, "%d ", cost_matrix[i][j]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "Assignment: ");
        for (int i = 0; i < N; i++) {
            fprintf(stderr, "%d(cost=%d)\n", assignment[i], cost_matrix[i][assignment[i]]);
        }
    }

    std::vector<int>& GetAssignment() {
        return assignment;
    }

    int GetTotalCost() {
        return total_cost;
    }
};