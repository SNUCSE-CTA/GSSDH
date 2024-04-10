#pragma once
#include "DataStructure/LabeledGraph.h"
#include "WeisfeilerLehman.h"
#include "Base/Logger.h"
#include "Base/Timer.h"
#include "State.h"

/**
 * @brief For two multisets A and B with integer elements in range [0, sz),
 *        maintain the set-difference-size between them efficiently.
 *        Used for graph similarity search.
 */
struct DifferenceVector {
    std::vector<int> value;
    int pos = 0, neg = 0;
    DifferenceVector(int sz = 0) {value.resize(sz, 0);}

    void init(int sz) {value.resize(sz, 0); pos = neg = 0;}
    void reset() { std::fill(value.begin(), value.end(), 0); pos = neg = 0;}

    unsigned long size() const {return value.size();}

    int& operator[](int index) {
        return value[index];
    }

    // val should be +1 or -1
    void update(int idx, int val) {
        if (val > 1) {
            for (int i = 0; i < val; i++) update(idx, 1);
            return;
        }
        if (val < -1) {
            for (int i = 0; i < -val; i++) update(idx, -1);
            return;
        }
        if (value[idx] == 0) {
            (val > 0 ? pos : neg)++;
        }
        else if (value[idx] < 0) {
            neg += (val > 0 ? -1 : 1);
        }
        else if (value[idx] > 0) {
            pos += (val > 0 ? 1 : -1);
        }
        value[idx] += val;
    }

    int GetDifference() {
        return std::max(pos, neg);
    }
};
