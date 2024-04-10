#include <iostream>
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "GraphSimilarity/GraphSimilaritySearch.h"
#include "GraphSimilarity/EditDistance.h"
#include "DataStructure/Graph.h"
using namespace std;
using namespace GraphLib;

Timer timer;
std::string data_root = "../data/GraphSimilaritySearch/", dataset = "aids100/";
std::string dataset_path, query_path;

std::vector<std::string> log_entries = {
    "EditDistance", "AStarNodes", "MaxQueueSize", "InitializeTime", "AStarTime"
};
int main(int argc, char *argv[]) {
    dataset_path = data_root+dataset+"8929.txt";
    query_path = data_root+dataset+"16881.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd':
                    dataset_path = argv[i + 1];
                break;
                case 'q':
                    query_path = argv[i + 1];
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
    ED.GCombined = new LabeledGraphDatabaseEntry();
    ED.Initialize(data[0], queries[0], -1);
    ED.AStar();
    auto log = ED.GetLog();
    log.PrintResults(log_entries);
    return 0;
}