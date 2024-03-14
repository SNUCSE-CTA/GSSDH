#include <iostream>
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "DataStructure/LabeledGraphDatabase.h"
#include "GraphSimilarity/EditDistance.h"
#include "DataStructure/Graph.h"
using namespace std;
using namespace GraphLib;

Timer timer;
std::string data_root = "../data/GraphSimilaritySearch/";
int tau = 4;
std::vector<ResultLogger> logs;
std::vector<std::string> log_entries = {
        "Ans", "NumCandidates", "TotalFilteringTime", "TotalAStarTime", "TotalSearchSpace"
};

// 24576, 1
int main(int argc, char *argv[]) {
    std::string dataset = "AIDS";
    std::string dataset_path = data_root + "AIDS.txt";
    std::string query_path = data_root + "AIDS_query100.txt";

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd':
                    dataset_path = argv[i + 1];
                break;
                case 'q':
                    query_path = argv[i + 1];
                break;
                case 't':
                    tau = atoi(argv[i + 1]);
                break;
                default:break;
            }
        }
    }

    LabeledGraphDatabase DB;
    DB.LoadGraphDatabase(dataset_path, 1);
    DB.LoadGraphDatabase(query_path, 0);

    auto data = DB.GetData();
    auto queries = DB.GetQueries();

    EditDistanceSolver ED(DB);
    int total_num_candidates = 0, ans = 0;
    for (int i = 0; i < queries.size(); i++) {
        for (int j = 0; j < data.size(); j++) {
            bool verify = ED.Initialize(queries[i], data[j], tau);
            if (ED.GetCurrentBestGED() <= tau) {
                total_num_candidates++;
                ans++; continue;
            }
            if (verify) {
                total_num_candidates++;
                int ged = ED.AStar();
                if (ged != -1) {
                    ans++;
                }
            }
            ResultLogger log = ED.GetLog();
            logs.push_back(log);
        }
    }
    ResultLogger aggregated_log;
    aggregated_log.AddResult("TotalFilteringTime", Total(logs, "FilteringTime"), RESULT_DOUBLE_FIXED);
    aggregated_log.AddResult("TotalAStarTime", Total(logs, "AStarTime"), RESULT_DOUBLE_FIXED);
    aggregated_log.AddResult("NumCandidates", total_num_candidates, RESULT_INT);
    aggregated_log.AddResult("Ans", ans, RESULT_INT);
    aggregated_log.AddResult("TotalSearchSpace", (int64_t)Total(logs, "AStarNodes") , RESULT_INT64);
    aggregated_log.AddResult("TotalMaxQueueSize", (int64_t)Total(logs, "MaxQueueSize"), RESULT_INT64);
    aggregated_log.PrintResults();
    return 0;
}