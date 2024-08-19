from utility import *
import pandas as pd
import ast

AIDS_AStarLSa_result = pd.read_csv("./reference/GED-AIDS-AStar-LSa.csv")


def run(args):
    q_l, q_r = args
    r = {"q_left": q_l, "q_right": q_r}
    binary = "GraphEditDistance"
    q_l = f"../data/GraphSimilaritySearch/aids100/{q_l}.txt"
    q_r = f"../data/GraphSimilaritySearch/aids100/{q_r}.txt"
    command = f"../build/{binary} -d {q_l} -q {q_r}"
    if platform.system() == "Darwin":
        cmd = f"(/usr/bin/time -l {command}) 2>&1"
    else:
        cmd = f"(/usr/bin/time -v {command}) 2>&1"
    with Popen(
        cmd,
        shell=True,
        stdout=PIPE,
        stderr=sys.stderr,
        universal_newlines=True,
        bufsize=1,
    ) as process:
        for _line in iter(process.stdout.readline, ""):
            line = _line.strip()
            if "[Result]" in line:
                _, key, _, value = line.strip().split()
                r[key] = to_number(value)
            if "maximum resident set size" in line:
                maxrss = to_number(line.split()[0])
            if "Maximum resident set size" in line:
                maxrss = to_number(line.split()[-1])
    r["MaxRSS(MB)"] = round(
        maxrss / ((2**20) if platform.system() == "Darwin" else (2**10)), 2
    )
    process.wait()
    print(r)
    return r


exp_ident = ""
if len(sys.argv) >= 2:
    exp_ident = "-" + sys.argv[1]

import itertools, multiprocessing

if __name__ == "__main__":
    os.makedirs("../build", exist_ok=True)
    os.system(f"cd ../build && cmake .. && make GraphEditDistance")
    os.makedirs("results/GraphEditDistance", exist_ok=True)
    gitver = get_current_repository_ver()[-8:]
    currentdate = get_experiment_date().replace("-", "")
    AIDS_AStarLSa_result.sort_values(by=["edit_distance", "query"], inplace=True)
    queries = AIDS_AStarLSa_result["query"].apply(ast.literal_eval).tolist()[:100]
    P = multiprocessing.Pool(processes=2)
    results = []
    for result in P.imap_unordered(run, queries):
        results.append(result)
    agg_result = pd.DataFrame(results)
    agg_result.to_csv(
        f"results/GraphEditDistance/{currentdate}-{gitver}{exp_ident}.csv", index=False
    )
