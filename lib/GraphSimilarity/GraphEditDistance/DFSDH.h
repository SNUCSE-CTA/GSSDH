#pragma once
#include "Base/Timer.h"
#include "GraphSimilarity/EditDistance.h"

namespace GraphLib::GraphSimilarity
{
class DFSDH : public GraphEditDistanceSolver
{
  public:
    std::vector<std::vector<int>> matrix;
    std::vector<int> mapping; // real mapping
    std::vector<int> inverse_mapping;
    std::vector<std::vector<int>> assignment; // for hungarian, assignment[depth][u]
    std::vector<std::vector<int>> inverse_assignment;
    std::vector<std::vector<int>> parikh1;
    std::vector<std::vector<int>> parikh2;

    std::vector<char> left_visited;
    std::vector<char> right_visited;
    std::vector<std::vector<int>> alpha;
    std::vector<std::vector<int>> beta;

    std::vector<bool> row;
    std::vector<bool> col;

    char *visX;
    char *visY;
    int *slack;
    int *slackmy;
    unsigned int *prev;
    unsigned int *queue;
    bool flag = false;

    const int INF = 1e3;
    const int INF2 = 1e5;
    int N = 0;
    int acc = 0;
    int total_cost = 0;
    // int lb = 0, ub = 0;
    int depth = -1;
    int tau = -1;

    int cost = 0; // mapping cost

    bool init = false;

    double hgtime = 0.0;
    double bdtime = 0.0;
    int64_t cnt = 0;

    void Initialize()
    {
        matrix = std::vector<std::vector<int>>(G2->GetNumVertices(), std::vector<int>(G2->GetNumVertices(), 0));
        mapping = std::vector<int>(G2->GetNumVertices(), -1);
        inverse_mapping = std::vector<int>(G2->GetNumVertices(), -1);
        assignment = std::vector<std::vector<int>>(G1->GetNumVertices() + 1,
                                                   std::vector<int>(G2->GetNumVertices(), -1)); // hungarian
        inverse_assignment =
            std::vector<std::vector<int>>(G2->GetNumVertices() + 1, std::vector<int>(G2->GetNumVertices(), -1));
        alpha = std::vector<std::vector<int>>(G1->GetNumVertices() + 1, std::vector<int>(G2->GetNumVertices(), 0));
        beta = std::vector<std::vector<int>>(G1->GetNumVertices() + 1, std::vector<int>(G2->GetNumVertices(), 0));
        row = std::vector<bool>(G1->GetNumVertices(), false);
        col = std::vector<bool>(G2->GetNumVertices(), false);
        parikh1.assign(G1->GetNumVertices(),
                       std::vector<int>(std::max(G1->GetNumEdgeLabels(), G2->GetNumEdgeLabels()), 0));
        parikh2.assign(G2->GetNumVertices(),
                       std::vector<int>(std::max(G1->GetNumEdgeLabels(), G2->GetNumEdgeLabels()), 0));
        N = G2->GetNumVertices();
        visX = new char[N];
        visY = new char[N];
        slack = new int[N];
        slackmy = new int[N];
        prev = new unsigned int[N];
        queue = new unsigned int[N];
        depth = -1;
        acc = 0;
        total_cost = 0;
        // lb = 0;
        // ub = 0;
        depth = -1;
        cost = 0; // mapping cost
        flag = false;
        cnt = 0;
        init = false;
    }
    using ui = unsigned int;

    void IterativeDeepningSearch()
    {
        for (int i = 0; i <= threshold; i++)
        {
            tau = i;
            init = false;
            DFS(-1, -1);
            if (flag)
            {
                break;
            }
        }
    }

    int Hungarian(char initialization)
    {
        int *mx = assignment[depth].data();
        int *my = inverse_assignment[depth].data();
        int *lx = alpha[depth].data();
        int *ly = beta[depth].data();
        ui n = N;
        if (initialization)
        { // Initialization
            memset(mx, -1, sizeof(int) * n);
            memset(my, -1, sizeof(int) * n);
            memset(ly, 0, sizeof(int) * n);
            for (ui i = 0; i < n; i++)
            {
                lx[i] = INF;
                for (ui j = 0; j < n; j++)
                    if (matrix[i][j] < lx[i])
                        lx[i] = matrix[i][j];
                for (ui j = 0; j < n; j++)
                    if (my[j] == -1 && matrix[i][j] == lx[i])
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
            if (mapping[i] == -1)
            {
                acc += lx[i];
            }
        }
        for (int j = 0; j < n; j++)
        {
            if (inverse_mapping[j] == -1)
            {
                acc += ly[j];
            }
        }
        int lb = cost + ((acc + 1) / 2);
        if (lb > tau)
        {
            return acc;
        }

        for (int u = n - 1; u >= 0; u--)
            if (mx[u] == -1)
            { // Augmentation
                memset(visX, 0, sizeof(char) * n);
                memset(visY, 0, sizeof(char) * n);
                int q_n = 1;
                queue[0] = u;
                visX[u] = 1;
                for (ui i = 0; i < n; i++)
                {
                    slack[i] = matrix[u][i] - lx[u] - ly[i];
                    slackmy[i] = u;
                }
                int target = n, X;
                while (true)
                {
                    for (ui i = 0; i < q_n && target == n; i++)
                    {
                        ui v = queue[i];
                        for (ui j = 0; j < n; j++)
                            if (!visY[j] && matrix[v][j] == lx[v] + ly[j])
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

                                for (ui k = 0; k < n; k++)
                                    if (!visY[k] && matrix[X][k] - lx[X] - ly[k] < slack[k])
                                    {
                                        slack[k] = matrix[X][k] - lx[X] - ly[k];
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
                    lb = cost + ((acc + 1) / 2);
                    if (lb > tau)
                    {
                        return acc;
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
                                    if (!visY[k] && matrix[X][k] - lx[X] - ly[k] < slack[k])
                                    {
                                        slack[k] = matrix[X][k] - lx[X] - ly[k];
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
        return acc;
    }

    void ChangeAlphaBeta(std::vector<bool> &row, std::vector<bool> &col)
    {
        for (int i = 0; i < row.size(); i++)
        {
            if (row[i] == true)
            {
                alpha[depth][i] = INF;
                for (int j = 0; j < G2->GetNumVertices(); j++)
                {
                    alpha[depth][i] = std::min(alpha[depth][i], matrix[i][j] - beta[depth][j]);
                }
            }
        }
        for (int j = 0; j < col.size(); j++)
        {
            if (col[j] == true)
            {
                beta[depth][j] = INF;
                for (int i = 0; i < G2->GetNumVertices(); i++)
                {
                    beta[depth][j] = std::min(beta[depth][j], matrix[i][j] - alpha[depth][i]);
                }
            }
        }
        for (int i = 0; i < G2->GetNumVertices(); i++)
        {
            int u = i;
            int v = assignment[depth][u];
            if (v != -1 && alpha[depth][u] + beta[depth][v] != matrix[u][v])
            {
                assignment[depth][u] = -1;
                inverse_assignment[depth][v] = -1;
            }
        }
    }

    void ComputeBranchDistanceMatrixInitial()
    {
        DifferenceVector diff;
        diff.init(20);
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
                    matrix[u][v] += 2;
                }
                for (int l = 0; l < u_nbrs.size(); l++)
                {
                    int u_nbr = u_nbrs[l];
                    int u_el = G1->GetEdgeLabel(u, u_nbr);
                    diff.update(u_el, 1);
                }
                for (int r = 0; r < v_nbrs.size(); r++)
                {
                    int v_nbr = v_nbrs[r];
                    int v_el = G2->GetEdgeLabel(v, v_nbr);
                    diff.update(v_el, -1);
                }
                int inner_distance = diff.GetDifference();
                matrix[u][v] += inner_distance;
            }
            int from_null = BranchEditDistanceFromNull(G2->GetBranch(v));
            for (; u < colSize; u++)
            {
                matrix[u][v] = from_null;
            }
        }
    }
    void ComputeBranchDistanceMatrixDynamic(const int u, const int v, bool update)
    {
        auto &u_nbrs = G1->GetNeighbors(u);
        auto &v_nbrs = G2->GetNeighbors(v);
        fill(row.begin(), row.end(), 0);
        fill(col.begin(), col.end(), 0);
        /* Update cost matrix */
        // CASE1 and CASE2: _u in Nbr(u)
        for (int i = 0; i < (int)u_nbrs.size(); ++i)
        {
            const int _u = u_nbrs[i];
            if (mapping[_u] != -1)
            {
                continue;
            }
            for (int _v = 0; _v < (int)G2->GetNumVertices(); ++_v)
            {
                if (inverse_mapping[_v] != -1)
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
                    if (update)
                    {
                        matrix[_u][_v] += delta;
                    }
                    else
                    {
                        matrix[_u][_v] -= delta;
                    }
                }
            }
            row[_u] = true;
        }
        // CASE3: _u not in Nbr(u) and _v in Nbr(v)
        for (int _u = 0; _u < (int)G1->GetNumVertices(); ++_u)
        {
            if (mapping[_u] != -1 || G1->GetEdgeLabel(u, _u) != -1)
            {
                continue;
            }
            for (int j = 0; j < (int)v_nbrs.size(); ++j)
            {
                const int _v = v_nbrs[j];
                if (inverse_mapping[_v] != -1)
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
                    if (update)
                    {
                        matrix[_u][_v] += delta;
                    }
                    else
                    {
                        matrix[_u][_v] -= delta;
                    }
                }
            }
        }

        for (int j = 0; j < (int)v_nbrs.size(); ++j)
        {
            const int _v = v_nbrs[j];
            if (inverse_mapping[_v] != -1)
            {
                continue;
            }
            for (int _u = G1->GetNumVertices(); _u < G2->GetNumVertices(); ++_u)
            {
                if (update)
                {
                    matrix[_u][_v] += 1;
                }
                else
                {
                    matrix[_u][_v] -= 1;
                }
            }
            col[_v] = true;
        }
        if (update)
        {
            ChangeAlphaBeta(row, col);
        }
    }

    void ComputeParikhVector()
    {
        const int n1 = G1->GetNumVertices();
        const int n2 = G2->GetNumVertices();
        const int np = std::max(G1->GetNumEdgeLabels(), G2->GetNumEdgeLabels());
        for (int u = 0; u < n1; ++u)
        {
            for (const int _u : G1->GetNeighbors(u))
                if (mapping[_u] == -1)
                {
                    const int u_u = G1->GetEdgeLabel(u, _u);
                    ++parikh1[u][u_u];
                    ++parikh1[u][0];
                }
        }
        for (int v = 0; v < n2; ++v)
        {
            for (const int _v : G2->GetNeighbors(v))
                if (inverse_mapping[_v] == -1)
                {
                    const int v_v = G2->GetEdgeLabel(v, _v);
                    ++parikh2[v][v_v];
                    ++parikh2[v][0];
                }
        }
    }
    void UpdateParikhVector(const int u, const int v)
    {
        auto &u_nbrs = G1->GetNeighbors(u);
        auto &v_nbrs = G2->GetNeighbors(v);
        for (auto _u : u_nbrs)
        {
            if (mapping[_u] == -1)
            {
                const int u_u = G1->GetEdgeLabel(u, _u);
                parikh1[_u][u_u]--;
                parikh1[_u][0]--;
            }
        }
        for (auto _v : v_nbrs)
        {
            if (inverse_mapping[_v] == -1)
            {
                const int v_v = G2->GetEdgeLabel(v, _v);
                parikh2[_v][v_v]--;
                parikh2[_v][0]--;
            }
        }
    }
    void RestoreParikhVector(const int u, const int v)
    {
        auto &u_nbrs = G1->GetNeighbors(u);
        auto &v_nbrs = G2->GetNeighbors(v);
        for (auto _u : u_nbrs)
        {
            if (mapping[_u] == -1)
            {
                const int u_u = G1->GetEdgeLabel(u, _u);
                parikh1[_u][u_u]++;
                parikh1[_u][0]++;
            }
        }
        for (auto _v : v_nbrs)
        {
            if (inverse_mapping[_v] == -1)
            {
                const int v_v = G2->GetEdgeLabel(v, _v);
                parikh2[_v][v_v]++;
                parikh2[_v][0]++;
            }
        }
    }

    int ChildEditCost(int u, int v)
    {
        // int child_cost = cost;
        int child_cost = 0;
        int u_label = G1->GetVertexLabel(u), v_label = G2->GetVertexLabel(v);
        if (u_label != v_label)
        {
            child_cost++;
        }
        int num_u_edges = 0;
        for (int u_nbr : G1->GetNeighbors(u))
        {
            if (mapping[u_nbr] != -1)
            {
                num_u_edges++;
            }
        }
        int ec = num_u_edges;
        for (int vprime : G2->GetNeighbors(v))
        {
            int uprime = inverse_mapping[vprime];
            if (uprime == -1)
                continue;
            ec++;
            int l1 = G1->GetEdgeLabel(u, uprime);
            int l2 = G2->GetEdgeLabel(v, vprime);
            if (l1 == -1)
                continue;
            if (l1 == l2)
                ec -= 2;
            else
                ec--;
        }
        child_cost += ec;
        return child_cost;
    }

    std::pair<int, int> LowerBound()
    {
        int ub = tau + 1, lb = 0;
        if (!init)
        {
            // if(tau == 0) {total_cost = Hungarian(1);}
            // else{
            // total_cost = Hungarian(0);
            // }
            total_cost = Hungarian(1);
            init = true;
        }
        else
        {
            total_cost = Hungarian(0);
        }
        lb = cost + ((total_cost + 1) / 2);
        // std::cout<<depth<<std::endl;
        // std::cout<<assignment[depth]<<std::endl;
        // std::cout<<inverse_assignment[depth]<<std::endl;
        if (lb <= tau)
        {
            ub = ComputeDistance(assignment[depth], inverse_assignment[depth]);
        }
        // std::cout<<lb<<' '<<ub<<std::endl;
        return {lb, ub};
    }

    void ChangeCostINF(int u, int v)
    {
        matrix[u][v] += INF;
        if (assignment[depth][u] == v)
        {
            assignment[depth][u] = -1;
            inverse_assignment[depth][v] = -1;
        }
    }
    void Match(int u, int v)
    {
        for (int j = 0; j < N; j++)
        {
            if (j != v)
            {
                ChangeCostINF(u, j);
            }
        }
        for (int i = 0; i < N; i++)
        {
            if (i != u)
            {
                ChangeCostINF(i, v);
            }
        }
    }
    void RemoveMatch(int u, int v)
    { // have to fix
        for (int j = 0; j < N; j++)
        {
            if (j != v)
            {
                matrix[u][j] -= INF;
            }
        }
        for (int i = 0; i < N; i++)
        {
            if (i != u)
            {
                matrix[i][v] -= INF;
            }
        }
    }

    // stop dfs
    void DFS(int u, int v)
    {
        cnt++;
        depth++;
        int chcost = 0;
        if (depth == G1->GetNumVertices())
        {
            depth--;
            return;
        }
        // std::cout << u << " " << v <<  " " << flag <<'\n';
        if (depth != 0)
        {
            assignment[depth] = assignment[depth - 1]; // copy from parent
            inverse_assignment[depth] = inverse_assignment[depth - 1];
            alpha[depth] = alpha[depth - 1];
            beta[depth] = beta[depth - 1];
            mapping[u] = v;
            inverse_mapping[v] = u;
            Match(u, v);
            chcost = ChildEditCost(u, v);
            cost += chcost;
            ComputeBranchDistanceMatrixDynamic(u, v, 1);
            UpdateParikhVector(u, v);
        }
        auto [lb, ub] = LowerBound();
        // std::cout <<u << " " << v << " " << lb << " " <<ub << "\n";
        int uprime = matching_order[depth];
        for (int i = 0; i < G2->GetNumVertices() - depth; i++)
        {
            if (i != 0)
            {
                std::tie(lb, ub) = LowerBound();
            }
            // std::cout<<"?? "<<' '<<lb<<' '<<ub<<std::endl;

            int vprime = assignment[depth][uprime];
            // if(uprime == 10){
            // std::cout<< uprime<<' '<<vprime<<' '<<lb<<' '<<ub<<std::endl;
            // std::cout<< assignment[depth] << "\n";
            // PrintMatrix();
            // }
            if (lb > tau)
            {
                break;
            }
            if (ub <= tau)
            {
                flag = true;
                return;
            }
            DFS(uprime, vprime);
            // std::cout<< "ret "<< uprime<<' '<<vprime<<' '<<lb<<' '<<ub<<std::endl;

            if (flag)
                return;
            matrix[uprime][vprime] += INF2;
            assignment[depth][uprime] = -1;
            inverse_assignment[depth][vprime] = -1;
        }
        for (int i = 0; i < G2->GetNumVertices(); i++)
        {
            if (matrix[uprime][i] >= INF2)
            {
                matrix[uprime][i] -= INF2;
            }
        }

        if (depth != 0)
        {
            RestoreParikhVector(u, v);
            ComputeBranchDistanceMatrixDynamic(u, v, 0);
            mapping[u] = -1;
            inverse_mapping[v] = -1;
            RemoveMatch(u, v);
            cost -= chcost;
        }
        depth--;
    }
    int GED()
    {
        PrepareGED(nullptr);
        Initialize();
        ComputeBranchDistanceMatrixInitial();
        ComputeParikhVector();
        // std::cout << parikh1 << "\n";
        // PrintMatrix();
        // exit(0);
        tau = threshold;
        DFS(-1, -1);
        // IterativeDeepningSearch();
        if (flag)
        {
            return 0;
        }
        else
            return -1;
    }

    void PrintMatrix()
    {
        for (int i = 0; i < G2->GetNumVertices(); i++)
        {
            for (int j = 0; j < G2->GetNumVertices(); j++)
            {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "-------------------------------------------------\n";
    };
    int64_t GetCnt() const
    {
        return cnt;
    }
    double Gethgtime() const
    {
        return hgtime;
    }
    double Getbdtime() const
    {
        return bdtime;
    }
    int64_t Getcnt() const
    {
        return cnt;
    }
};

} // namespace GraphLib::GraphSimilarity