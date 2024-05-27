#include <iostream>
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "GraphSimilarity/GraphSimilaritySearch.h"
#include "GraphSimilarity/EditDistance.h"
#include "DataStructure/Graph.h"
using namespace std;
using namespace GraphLib::GraphSimilarity;

Timer timer;
std::string data_root = "../data/GraphSimilaritySearch/";
int tau = 4;
std::vector<ResultLogger> logs;
std::vector<std::string> log_entries = {
        "Ans", "NumCandidates", "TotalFilteringTime", "TotalAStarTime", "TotalSearchSpace"
};

int main(int argc, char *argv[]) {
    std::string dataset = "AIDS";
    std::string dataset_path = data_root + "AIDS.database";
    std::string query_path = data_root + "AIDS.100query";

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

    GraphSimilaritySearch DB;
    DB.LoadGraphDatabase(dataset_path, -1);
//    auto data = DB.GetData();
    DB.BuildBranches();
    DB.LoadGraphDatabase(query_path, tau);
    DB.GetLog().PrintResults();
    return 0;
}