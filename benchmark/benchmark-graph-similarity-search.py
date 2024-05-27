from utility import *
import pandas as pd 
import ast 

def run(args):
    binary = "GraphSimilaritySearch"
    data, query, tau = args
    command = f"../build/{binary} -d {data} -q {query} -t {tau} > ./results/run_{tau}.out"
    if platform.system() == "Darwin":
        cmd = f"(/usr/bin/time -l {command}) 2>&1"
    else:
        cmd = f"(/usr/bin/time -v {command}) 2>&1"
    r = {"tau": tau}
    with Popen(cmd, shell=True, stdout=PIPE, stderr=sys.stderr, universal_newlines=True,bufsize=1) as process:
        for _line in iter(process.stdout.readline, ''):
            line = _line.strip()
            if "[Result]" in line:
                _, key, _, value = line.strip().split()
                r[key] = to_number(value)
            if "maximum resident set size" in line:
                maxrss = to_number(line.split()[0])
            if "Maximum resident set size" in line:
                maxrss = to_number(line.split()[-1])
    process.wait()
    r["MaxRSS(MB)"] = round(maxrss / ((2**20) if platform.system() == "Darwin" else (2**10)), 2)
    print(r)
    return r

exp_ident = ""
if len(sys.argv) >= 2:
    exp_ident = "-"+sys.argv[1]

import itertools, multiprocessing
if __name__ == '__main__':
    os.makedirs("../build", exist_ok=True)
    os.system(f"cd ../build && cmake .. && make GraphSimilaritySearch")
    os.makedirs("results/GraphSimilaritySearch", exist_ok=True)
    gitver = get_current_repository_ver()[-8:]
    currentdate = get_experiment_date().replace("-", "")
    data_path = "../data/GraphSimilaritySearch/AIDS.database"
    query_path = "../data/GraphSimilaritySearch/AIDS.100query"
    workloads = [(data_path, query_path, tau) for tau in [7, 6, 5, 4, 3, 2, 1, 0]]
    P = multiprocessing.Pool(processes=3)
    results = []
    for result in P.imap_unordered(run, workloads):
        results.append(result)
    agg_result = (pd.DataFrame(results)).sort_values(by=['tau'])
    agg_result.to_csv(f"results/GraphSimilaritySearch/{currentdate}-{gitver}{exp_ident}.csv", index=False)
    agg_result.to_csv(f"results/GraphSimilaritySearch/last-results.csv", index=False)