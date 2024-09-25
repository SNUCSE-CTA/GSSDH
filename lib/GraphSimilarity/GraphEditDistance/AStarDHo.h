#pragma once
#include "Base/Timer.h"
#include "GraphSimilarity/EditDistance.h"

namespace GraphLib::GraphSimilarity
{
struct DHoState : State
{
    int next_mapping_order = -1;
    int ub = 1e9;
    std::vector<std::vector<int>> matrix;
    std::vector<int> hungarian_assignment;
    std::vector<int> hungarian_inverse_assignment;
    std::vector<int> alpha, beta;
};

std::vector<DHoState *> states;

class AStarDHo : public GraphEditDistanceSolver
{
  public:
    struct StateIndexComparator
    {
        bool operator()(const int a, const int b) const
        {
            if (states[a]->lower_bound == (states[b])->lower_bound)
            {
                if (states[a]->depth == states[b]->depth)
                {
                    return states[a]->id < states[b]->id;
                }
                return states[a]->depth < states[b]->depth;
            }
            return states[a]->lower_bound > states[b]->lower_bound;
        }
    };

    std::priority_queue<int, std::vector<int>, StateIndexComparator> queue;
    std::stack<int> stk;
    const int INF = 1e9;

    std::vector<int> assignment, inverse_assignment;
    std::vector<int> left_visited, right_visited;
    std::vector<int> alpha, beta;
    int total_cost, N, theta;
    int cnt = 0;

    std::vector<int> __left;
    std::vector<int> __right1, __right2;
    std::vector<std::vector<int>> parikh1;
    std::vector<std::vector<int>> parikh2;

    int acc = 0;
    int64_t functioncall = 0;

    double hgtime = 0.0;
    double bdtime = 0.0;

    std::vector<int> u_idxs, v_idxs; // u_idxs[u] = u_idx, v_idxs[v] = v_idx;

  public:
    void Initialize(DHoState *state)
    {
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
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                state->beta[j] = std::min(state->beta[j], state->matrix[i][j]);
            }
        }
        alpha = state->alpha;
        beta = state->beta;

        using std::vector;
        parikh1.assign(G1->GetNumVertices(), vector<int>(std::max(G1->GetNumEdgeLabels(), G2->GetNumEdgeLabels()), 0));
        parikh2.assign(G2->GetNumVertices(), vector<int>(std::max(G1->GetNumEdgeLabels(), G2->GetNumEdgeLabels()), 0));
    }

    using ui = unsigned int;
    int Hungarian(char initialization, ui n, std::vector<std::vector<int>> &cost_matrix, DHoState *state,
                  std::vector<int> &rem_left, std::vector<int> &rem_right)
    {
        int *mx = assignment.data();
        int *my = inverse_assignment.data();
        int *lx = alpha.data();
        int *ly = beta.data();
        char *visX = new char[n];
        char *visY = new char[n];
        int *slack = new int[n];
        int *slackmy = new int[n];
        ui *prev = new ui[n];
        ui *queue = new ui[n];

        if (initialization)
        { // Initialization
            memset(mx, -1, sizeof(int) * n);
            memset(my, -1, sizeof(int) * n);
            memset(ly, 0, sizeof(int) * n);
            for (ui i = 0; i < n; i++)
            {
                lx[i] = INF;
                // ui *t_array = cost+i*n;
                for (ui j = 0; j < n; j++)
                    if (cost_matrix[i][j] < lx[i])
                        lx[i] = cost_matrix[i][j];
                for (ui j = 0; j < n; j++)
                    if (my[j] == -1 && cost_matrix[i][j] == lx[i])
                    {
                        mx[i] = j;
                        my[j] = i;
                        break;
                    }
            }
        }

        acc = 0;
        for (int i = 0; i < n; i++)
        {
            acc += alpha[i];
            // std::cout << acc << "\n";
        }
        for (int j = 0; j < n; j++)
        {
            acc += beta[j];
            // std::cout << acc << "\n";
        }

        for (int u = n - 1; u >= 0; u--)
            if (mx[u] == -1)
            { // Augmentation
                int lb = state->cost + ((acc + 1) / 2);
                if (lb > threshold)
                {
                    total_cost = 2 * threshold + 1000;
                    for (int i = 0; i < n; i++)
                    {
                        assignment[i] = inverse_assignment[i] = i;
                    }
                    return total_cost;
                }
                memset(visX, 0, sizeof(char) * n);
                memset(visY, 0, sizeof(char) * n);
                int q_n = 1;
                queue[0] = u;
                visX[u] = 1;
                // ui *t_array = cost + u*n;
                for (ui i = 0; i < n; i++)
                {
                    slack[i] = cost_matrix[u][i] - lx[u] - ly[i];
                    slackmy[i] = u;
                }
                int target = n, X;
                while (true)
                {
                    for (ui i = 0; i < q_n && target == n; i++)
                    {
                        ui v = queue[i];
                        // t_array = cost + v*n;
                        for (ui j = 0; j < n; j++)
                            if (!visY[j] && cost_matrix[v][j] == lx[v] + ly[j])
                            {
                                if (my[j] == -1)
                                {
                                    X = v;
                                    target = j;
                                    break;
                                }
                                visY[j] = 1;
                                X = my[j];
                                visX[X] = 1;
                                prev[X] = v;
                                queue[q_n++] = X;

                                // ui *tt_array = cost + X*n;
                                for (ui k = 0; k < n; k++)
                                    if (!visY[k] && cost_matrix[X][k] - lx[X] - ly[k] < slack[k])
                                    {
                                        slack[k] = cost_matrix[X][k] - lx[X] - ly[k];
                                        slackmy[k] = X;
                                    }
                            }
                    }
                    if (target != n)
                        break;

                    q_n = 0;
                    int delta = INF;
                    for (ui i = 0; i < n; i++)
                        if (!visY[i] && slack[i] < delta)
                            delta = slack[i];
                    for (ui i = 0; i < n; i++)
                    {
                        if (visX[i])
                        {
                            lx[i] += delta;
                            acc += delta;
                        }
                        if (visY[i])
                        {
                            ly[i] -= delta;
                            acc -= delta;
                        }
                        else
                            slack[i] -= delta;
                    }

                    for (ui i = 0; i < n; i++)
                        if (!visY[i] && slack[i] == 0)
                        {
                            if (my[i] == -1)
                            {
                                X = slackmy[i];
                                target = i;
                                break;
                            }
                            visY[i] = 1;
                            if (!visX[my[i]])
                            {
                                X = my[i];
                                visX[X] = 1;
                                prev[X] = slackmy[i];
                                queue[q_n++] = X;

                                // ui *tt_array = cost + X*n;
                                for (ui k = 0; k < n; k++)
                                    if (!visY[k] && cost_matrix[X][k] - lx[X] - ly[k] < slack[k])
                                    {
                                        slack[k] = cost_matrix[X][k] - lx[X] - ly[k];
                                        slackmy[k] = X;
                                    }
                            }
                        }
                }

                while (true)
                {
                    int ty = mx[X];
                    mx[X] = target;
                    my[target] = X;
                    if (X == u)
                        break;

                    X = prev[X];
                    target = ty;
                }
            }

        memset(visX, 0, sizeof(char) * n);
        memset(visY, 0, sizeof(char) * n);

        int res = 0;
        for (ui i = 0; i < n; i++)
        {
            res += cost_matrix[i][mx[i]];
            if (state->depth != -1)
            {
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
        delete[] visX;
        delete[] visY;
        delete[] slack;
        delete[] slackmy;
        delete[] prev;
        delete[] queue;
        return res;
    }

    // bool FindAugmentingPath(int i, std::vector<std::vector<int>>& cost_matrix){
    //     functioncall++;
    //     left_visited[i] = true;
    //     for(int j = 0 ; j < N; j++){
    //         if(!right_visited[j] &&(alpha[i] + beta[j] == cost_matrix[i][j])){
    //             right_visited[j] = true;
    //             if(inverse_assignment[j] == -1 ||
    //             FindAugmentingPath(inverse_assignment[j], cost_matrix)){
    //                 inverse_assignment[j] = i;
    //                 assignment[i] = j;
    //                 return true;
    //             }
    //         }
    //     }
    //     return false;
    // }
    bool FindAugmentingPath(int i, const std::vector<std::vector<int>> &cost_matrix)
    {
        // functioncall++;
        // left_visited[i] = true;
        __left.push_back(i);

        for (int j = 0; j < N; j++)
        {
            if (__builtin_expect(right_visited[j], true))
                continue;

            // for (int jdx = 0; jdx < (int)__right2.size(); ++jdx) {
            //     int j = __right2[jdx];

            if (__builtin_expect(alpha[i] + beta[j] != cost_matrix[i][j], true))
                continue;
            // functioncall++;
            right_visited[j] = true;
            __right1.push_back(j);
            //__right2[jdx] = __right2.back(); --jdx;
            //__right2.pop_back();

            if (inverse_assignment[j] == -1 || FindAugmentingPath(inverse_assignment[j], cost_matrix))
            {
                inverse_assignment[j] = i;
                assignment[i] = j;
                // functioncall++;
                return true;
            }
        }
        return false;
    }

    void RecalculatePotential(const std::vector<std::vector<int>> &cost_matrix)
    {
        // functioncall++;
        theta = INF;
        // for(int i = 0 ; i < N; i++){
        //     if(left_visited[i]){
        for (int i : __left)
        {
            for (int j = 0; j < N; j++)
            {
                if (!right_visited[j])
                {
                    // for (int j: __right2)
                    theta = std::min(theta, cost_matrix[i][j] - alpha[i] - beta[j]);
                    // functioncall++;
                }
            }
        }
        //    }
        //}
        // for(int i = 0 ; i < N; i++){
        //    if(left_visited[i]){
        for (int i : __left)
            alpha[i] += theta;
        //    }
        //}
        // for(int j = 0 ; j < N ; j++){
        //    if(right_visited[j]){
        for (int j : __right1)
            beta[j] -= theta;
        //    }
        //}
    }
    void Solve(std::vector<std::vector<int>> &cost_matrix, DHoState *state, std::vector<int> &rem_left,
               std::vector<int> &rem_right)
    {
        for (int i = N - 1; i >= 0; i--)
        {
            if (assignment[i] != -1)
                continue;
            functioncall++;
            while (true)
            {
                // std::fill(left_visited.begin(), left_visited.end(), 0);
                std::fill(right_visited.begin(), right_visited.end(), 0);

                __left.clear();

                __right1.clear();
                //__right2.resize(N);
                // for (int i = 0; i < N; ++i) __right2[i] =  i;

                // Timer t1;
                // t1.Start();
                bool flag = FindAugmentingPath(i, cost_matrix);
                // t1.Stop();
                // hgtime += t1.GetTime();
                if (flag)
                    break;
                // Timer t1;
                // t1.Start();

                RecalculatePotential(cost_matrix);
                //                     t1.Stop();
                // hgtime += t1.GetTime();
            }
        }
        total_cost = 0;
        for (int i = 0; i < N; i++)
        {
            total_cost += cost_matrix[i][assignment[i]];
            if (state->depth != -1)
            {
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
    void ChangeCost(int i, int j, int newCost, DHoState *state)
    {
        int oldCost = state->matrix[i][j];
        state->matrix[i][j] = newCost;
        if (newCost > oldCost && state->hungarian_assignment[i] == j)
        {
            state->hungarian_assignment[i] = -1;
            state->hungarian_inverse_assignment[j] = -1;
        }
    }
    void ComputeBranchDistanceMatrixInitial(DHoState *state)
    {
        DifferenceVector diff;
        diff.init(20);
        int N = state->matrix.size();
        int colSize = G2->GetNumVertices(); // right
        int rowSize = G1->GetNumVertices(); // left
        for (int v = 0; v < colSize; v++)
        {
            auto &v_nbrs = G2->GetNeighbors(v);
            int u = 0;
            for (u = 0; u < rowSize; u++)
            {
                auto &u_nbrs = G1->GetNeighbors(u);
                diff.reset();
                if (G1->GetVertexLabel(u) != G2->GetVertexLabel(v))
                {
                    state->matrix[u][v] += 2;
                }
                for (int l = 0; l < u_nbrs.size(); l++)
                {
                    int u_nbr = u_nbrs[l];
                    int u_el = G1->GetEdgeLabel(u, u_nbr);
                    if (state->mapping[u_nbr] == -1)
                    {
                        diff.update(u_el, 1);
                    }
                    else
                    {
                        int v_mapping_el = G2->GetEdgeLabel(v, state->mapping[u_nbr]);
                        if (v_mapping_el != u_el)
                        {
                            state->matrix[u][v] += 2;
                        }
                    }
                }
                for (int r = 0; r < v_nbrs.size(); r++)
                {
                    int v_nbr = v_nbrs[r];
                    int v_el = G2->GetEdgeLabel(v, v_nbr);
                    if (state->inverse_mapping[v_nbr] == -1)
                    {
                        diff.update(v_el, -1);
                    }
                    else
                    {
                        int u_mapping_el = G1->GetEdgeLabel(u, state->inverse_mapping[v_nbr]);
                        if (u_mapping_el == -1)
                        {
                            state->matrix[u][v] += 2;
                        }
                    }
                }
                int inner_distance = diff.GetDifference();
                state->matrix[u][v] += inner_distance;
            }
            int from_null = BranchEditDistanceFromNull(G2->GetBranch(v));
            for (int v_nbr : G2->GetNeighbors(v))
            {
                if (state->inverse_mapping[v_nbr] != -1)
                {
                    from_null++;
                }
            }
            for (; u < colSize; u++)
            {
                state->matrix[u][v] = from_null;
            }
        }
    }
    void ChangeAlphaBeta(DHoState *state, std::vector<bool> &row, std::vector<bool> &col)
    {
        for (int i = 0; i < row.size(); i++)
        {
            if (row[i] == true)
            {
                state->alpha[i] = INF;
                for (int j = 0; j < G2->GetNumVertices(); j++)
                {
                    state->alpha[i] = std::min(state->alpha[i], state->matrix[i][j] - state->beta[j]);
                }
            }
        }
        for (int j = 0; j < col.size(); j++)
        {
            if (col[j] == true)
            {
                state->beta[j] = INF;
                for (int i = 0; i < G2->GetNumVertices(); i++)
                {
                    state->beta[j] = std::min(state->beta[j], state->matrix[i][j] - state->alpha[i]);
                }
            }
        }
        for (int i = 0; i < G2->GetNumVertices(); i++)
        {
            int u = i;
            int v = state->hungarian_assignment[u];
            if (v != -1 && state->alpha[u] + state->beta[v] != state->matrix[u][v])
            {
                state->hungarian_assignment[u] = -1;
                state->hungarian_inverse_assignment[v] = -1;
            }
        }
    }

    void ComputeBranchDistanceMatrixDynamic(DHoState *state)
    {
        const int u = matching_order[state->depth];
        const int v = state->mapping[u];
        auto &u_nbrs = G1->GetNeighbors(u);
        auto &v_nbrs = G2->GetNeighbors(v);
        std::vector<bool> row(G1->GetNumVertices(), false);
        std::vector<bool> col(G2->GetNumVertices(), false);
        /* Compute Parikh vectors */
        using std::vector;
        const int n1 = G1->GetNumVertices();
        const int n2 = G2->GetNumVertices();
        const int np = std::max(G1->GetNumEdgeLabels(), G2->GetNumEdgeLabels());

        for (int u = 0; u < n1; ++u)
        {
            for (const int _u : G1->GetNeighbors(u))
                if (state->mapping[_u] == -1 || state->mapping[_u] == v)
                {
                    const int u_u = G1->GetEdgeLabel(u, _u);
                    ++parikh1[u][u_u];
                    ++parikh1[u][0];
                }
        }
        for (int v = 0; v < n2; ++v)
        {
            for (const int _v : G2->GetNeighbors(v))
                if (state->inverse_mapping[_v] == -1 || state->inverse_mapping[_v] == u)
                {
                    const int v_v = G2->GetEdgeLabel(v, _v);
                    ++parikh2[v][v_v];
                    ++parikh2[v][0];
                }
        }
        /* Update cost matrix */
        // CASE1 and CASE2: _u in Nbr(u)
        for (int i = 0; i < (int)u_nbrs.size(); ++i)
        {
            const int _u = u_nbrs[i];
            if (state->mapping[_u] != -1)
            {
                continue;
            }
            for (int _v = 0; _v < (int)G2->GetNumVertices(); ++_v)
            {
                if (state->inverse_mapping[_v] != -1)
                { // || G2->GetEdgeLabel(v, _v) != -1) {
                    continue;
                }

                const int l1 = G1->GetEdgeLabel(u, _u);
                const int l2 = G2->GetEdgeLabel(v, _v);

                if (l1 != l2)
                {
                    int delta = 0;
                    delta += 2;

                    delta -= ~(l1 | l2) || (parikh1[_u][0] > parikh2[_v][0]);
                    if (l1 != -1 && parikh1[_u][l1] <= parikh2[_v][l1])
                    {
                        delta += 1;
                    }
                    if (l2 != -1 && parikh1[_u][l2] >= parikh2[_v][l2])
                    {
                        delta += 1;
                    }
                    state->matrix[_u][_v] += delta;
                }
            }
            row[_u] = true;
        }

        // CASE3: _u not in Nbr(u) and _v in Nbr(v)
        for (int _u = 0; _u < (int)G1->GetNumVertices(); ++_u)
        {
            if (state->mapping[_u] != -1 || G1->GetEdgeLabel(u, _u) != -1)
            {
                continue;
            }
            for (int j = 0; j < (int)v_nbrs.size(); ++j)
            {
                const int _v = v_nbrs[j];
                if (state->inverse_mapping[_v] != -1)
                {
                    continue;
                }

                const int l1 = G1->GetEdgeLabel(u, _u);
                const int l2 = G2->GetEdgeLabel(v, _v);

                if (l1 != l2)
                {
                    int delta = 0;
                    delta += 2;

                    delta -= (parikh1[_u][0] < parikh2[_v][0]);
                    if (parikh1[_u][l2] >= parikh2[_v][l2])
                    {
                        delta += 1;
                    }
                    state->matrix[_u][_v] += delta;
                }
            }
        }

        for (int j = 0; j < (int)v_nbrs.size(); ++j)
        {
            const int _v = v_nbrs[j];
            if (state->inverse_mapping[_v] != -1)
            {
                continue;
            }
            for (int _u = G1->GetNumVertices(); _u < G2->GetNumVertices(); ++_u)
            {
                state->matrix[_u][_v] += 1;
            }
            col[_v] = true;
        }

        ChangeAlphaBeta(state, row, col);

        /* Release Parikh vectors */
        for (int u = 0; u < n1; ++u)
        {
            for (const int _u : G1->GetNeighbors(u))
                if (state->mapping[_u] == -1 || state->mapping[_u] == v)
                {
                    const int u_u = G1->GetEdgeLabel(u, _u);
                    parikh1[u][u_u] = 0;
                    parikh1[u][0] = 0;
                }
        }
        for (int v = 0; v < n2; ++v)
        {
            for (const int _v : G2->GetNeighbors(v))
                if (state->inverse_mapping[_v] == -1 || state->inverse_mapping[_v] == u)
                {
                    const int v_v = G2->GetEdgeLabel(v, _v);
                    parikh2[v][v_v] = 0;
                    parikh2[v][0] = 0;
                }
        }
    }

    void Match(DHoState *state)
    {
        int u = matching_order[state->depth];
        int v = state->mapping[u];
        for (int j = 0; j < N; j++)
        {
            if (j == v)
                continue;
            ChangeCost(u, j, INF, state);
        }
        for (int i = 0; i < N; i++)
        {
            if (i == u)
                continue;
            ChangeCost(i, v, INF, state);
        }
        state->hungarian_inverse_assignment[v] = u;
        state->hungarian_assignment[u] = v;
    }
    void ComputeReducedMatrix(DHoState *state, std::vector<std::vector<int>> &local_matrix, std::vector<int> &rem_left,
                              std::vector<int> &rem_right, int &remaining)
    {
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
        for (int u = 0; u < G1->GetNumVertices(); u++)
        {
            if (state->mapping[u] == -1)
            {
                u_idxs[u] = rem_left.size();
                rem_left.emplace_back(u);
                alpha[u_idxs[u]] = state->alpha[u];
            }
        }
        for (int u = G1->GetNumVertices(); u < G2->GetNumVertices(); u++)
        {
            u_idxs[u] = rem_left.size();
            rem_left.emplace_back(u);
            alpha[u_idxs[u]] = state->alpha[u];
        }
        for (int v = 0; v < G2->GetNumVertices(); v++)
        {
            if (state->inverse_mapping[v] == -1)
            {
                v_idxs[v] = rem_right.size();
                rem_right.emplace_back(v);
                beta[v_idxs[v]] = state->beta[v];
            }
        }
        for (int i = 0; i < remaining; i++)
        {
            if (state->hungarian_assignment[rem_left[i]] != -1)
            {
                assignment[i] = v_idxs[state->hungarian_assignment[rem_left[i]]];
            }
            if (state->hungarian_inverse_assignment[rem_right[i]] != -1)
            {
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
        for (int v_idx = 0; v_idx < rem_right.size(); v_idx++)
        {
            int v = rem_right[v_idx];
            for (int u_idx = 0; u_idx < rem_left.size(); u_idx++)
            {
                int u = rem_left[u_idx];
                local_matrix[u_idx][v_idx] = state->matrix[u][v];
            }
        }
        // std::cout << alpha << "\n" << beta << "\n" << assignment << "\n" <<
        // inverse_assignment << "\n";
    }
    void LocalToState(DHoState *state, std::vector<int> &rem_left, std::vector<int> &rem_right, int &remaining)
    {
        for (int i = 0; i < remaining; i++)
        {
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
    std::pair<int, int> DHoLowerBound(DHoState *state)
    {
        int ub = 0, lb = 0;
        // state->matrix.resize(G2->GetNumVertices(), std::vector<int>(G2->GetNumVertices(), 0));
        // state->hungarian_assignment.resize(G2->GetNumVertices(), -1);
        // state->alpha.resize(G2->GetNumVertices());
        // state->beta.resize(G2->GetNumVertices());
        N = G2->GetNumVertices() - (state->depth + 1);
        // std::cout << N << "\n";
        // Timer t;
        if (state->depth == -1)
        {   
            Timer ttt;
            ttt.Start();
            state->matrix.resize(G2->GetNumVertices(), std::vector<int>(G2->GetNumVertices(), 0));
            // t.Start();
            ComputeBranchDistanceMatrixInitial(state);
            // t.Stop();
            // bdtime += t.GetTime();
            Initialize(state);
            std::vector<int> rem_left, rem_right;
            // Timer b;
            // b.Start();
            total_cost = Hungarian(1, N, state->matrix, state, rem_left, rem_right);
            // b.Stop();
            // hgtime += b.GetTime();

            state->hungarian_assignment = assignment;
            state->hungarian_inverse_assignment = inverse_assignment;
            state->alpha = alpha;
            state->beta = beta;

            std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
            std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
            std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
            std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
            for (int i = 0; i < G1->GetNumVertices(); i++)
            {
                if (state->mapping[i] != -1)
                    continue;
                hungarian_mapping[i] = state->hungarian_assignment[i];
                hungarian_inverse_mapping[state->hungarian_assignment[i]] = i;
            }
            // t.Start();
            ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
            lb = state->cost + ((total_cost + 1) / 2);
            ttt.Stop();bdtime += ttt.GetTime();
        }
        else
        {   
            Timer ttt;
            ttt.Start();
            state->matrix = static_cast<DHoState *>(state->parent)->matrix;
            state->alpha = static_cast<DHoState *>(state->parent)->alpha;
            state->beta = static_cast<DHoState *>(state->parent)->beta;
            state->hungarian_assignment = static_cast<DHoState *>(state->parent)->hungarian_assignment;
            state->hungarian_inverse_assignment = static_cast<DHoState *>(state->parent)->hungarian_inverse_assignment;
            ttt.Stop();
            bdtime += ttt.GetTime();
            // t.Start();
            ComputeBranchDistanceMatrixDynamic(state);
            // Timer ttt;
            // ttt.Start();
            Match(state);
            // ttt.Stop();
            // bdtime += ttt.GetTime();
            // t.Stop();
            // bdtime += t.GetTime();

            int remaining = G2->GetNumVertices() - (state->depth + 1);
            std::vector<int> rem_left, rem_right;
            std::vector<std::vector<int>> local_matrix(remaining, std::vector<int>(remaining, 0));
            ComputeReducedMatrix(state, local_matrix, rem_left, rem_right, remaining);


            Timer a;
            a.Start();
            total_cost = Hungarian(0, N, local_matrix, state, rem_left, rem_right);
            a.Stop();
            hgtime += a.GetTime();
  
            std::vector<int> hungarian_mapping(G1->GetNumVertices(), -1);
            std::vector<int> hungarian_inverse_mapping(G2->GetNumVertices(), -1);
            std::memcpy(hungarian_mapping.data(), state->mapping, sizeof(int) * G1->GetNumVertices());
            std::memcpy(hungarian_inverse_mapping.data(), state->inverse_mapping, sizeof(int) * G2->GetNumVertices());
            for (int u_idx = 0; u_idx < rem_left.size() - (G2->GetNumVertices() - G1->GetNumVertices()); u_idx++)
            {
                int u = rem_left[u_idx];
                int v = rem_right[assignment[u_idx]];
                hungarian_mapping[u] = v;
                hungarian_inverse_mapping[v] = u;
            }
            ub = ComputeDistance(hungarian_mapping, hungarian_inverse_mapping);
            lb = state->cost + ((total_cost + 1) / 2);
        }
        if(lb == 3){
          std::cout << matching_order[state->depth] << " " << state->mapping[matching_order[state->depth]] << " " << lb << " " << ub << "\n";
          std::cout << state->hungarian_assignment << "\n";
        }
        return {lb, ub};
    }

    void ExtendState(DHoState *state)
    {
        if (state->cost >= current_best)
            return;
        int depth = state->depth + 1;
        int u = matching_order[depth];
        for (int v = 0; v < G2->GetNumVertices(); v++)
        {
            if (state->inverse_mapping[v] != -1)
                continue;
            DHoState *child_state = new DHoState(state);
            child_state->cost = GetChildEditCost(state, u, v);
            child_state->mapping[u] = v;
            child_state->inverse_mapping[v] = u;
            child_state->depth = depth;
            auto [lb, ub] = DHoLowerBound(child_state);
            child_state->lower_bound = lb;
            child_state->ub = ub;
            if (lb == ub)
            {
                delete child_state;
                continue;
            }
            if (depth == G1->GetNumVertices() - 1)
            {
                current_best = std::min(current_best, child_state->lower_bound);
                memcpy(&current_best_mapping[0], child_state->mapping, sizeof(int) * NumG1Vertices);
                continue;
            }
            if (child_state->lower_bound >= current_best)
            {
                delete child_state;
                continue;
            }
            if (threshold > 0)
            {
                if (child_state->lower_bound > threshold)
                { 
                    delete child_state;
                    continue;
                }
            }
            states.emplace_back(child_state);
            // num_nodes++;
            stk.push(states.size() -1);
            // queue.push(states.size() - 1);
        }
    }

    int GED()
    {
        PrepareGED(nullptr);
        DHoState *initial_state = new DHoState(NULL);
        hgtime = 0.0;
        bdtime = 0.0;
        cnt = 0;
        functioncall = 0;
        initial_state->cost = 0;
        initial_state->depth = -1;
        auto [lb, ub] = DHoLowerBound(initial_state);
        initial_state->lower_bound = lb;
        // num_nodes++;
        states.emplace_back(initial_state);
        queue.push(states.size() - 1);
        stk.push(states.size() -1);
        int64_t max_qsize = 1;
        // while (!queue.empty())
        while(!stk.empty())
        {
            // DHoState *current_state = states[queue.top()];
              DHoState *current_state = states[stk.top()];
            num_nodes++;
            // queue.pop();
            stk.pop();
            if (current_state->lower_bound >= current_best)
            {     
                // queue = std::priority_queue<int, std::vector<int>, StateIndexComparator>();
                stk = std::stack<int>();
                break;
            }
            if (threshold >= 0)
            {
                if (current_best < threshold)
                {
                    // queue = std::priority_queue<int, std::vector<int>, StateIndexComparator>();
                    stk = std::stack<int>();
                    break;
                }
                if (current_state->lower_bound > threshold)
                {
                    current_best = -1;
                    // queue = std::priority_queue<int, std::vector<int>, StateIndexComparator>();
                    stk = std::stack<int>();
                    break;
                }
            }
            ExtendState(current_state);
            max_qsize = std::max(max_qsize, (int64_t)queue.size());
        }
        // fprintf(stderr, "States_size = %d\n", states.size());
        for (int i = 0; i < states.size(); i++)
        {
            delete states[i];
        }
        states.clear();
        if (threshold >= 0 and current_best > threshold)
        {
            current_best = -1;
        }
        log.AddResult("MaxQueueSize", max_qsize, RESULT_INT64);
        log.AddResult("AStarNodes", num_nodes, RESULT_INT64);
        log.AddResult("EditDistance", current_best, RESULT_INT);
        return current_best;
    }

    void PrintCostMatrix(std::vector<std::vector<int>> &matrix)
    {
        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix.size(); j++)
            {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "-----------------------\n";
    }

    double Gethgtime() const
    {
        return hgtime;
    }
    double Getbdtime() const
    {
        return bdtime;
    }
    int64_t GetCnt() const
    {
        return functioncall;
    }
};
} // namespace GraphLib::GraphSimilarity
