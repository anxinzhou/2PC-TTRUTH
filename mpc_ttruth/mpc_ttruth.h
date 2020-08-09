//
// Created by anxin on 8/8/20.
//

#ifndef TTRUTH_MPC_TTRUTH_H
#define TTRUTH_MPC_TTRUTH_H

#include "common/mpc_util.h"

namespace MPC {
    extern const uint CLUSTER_NUM;

    vector<vector<uint64_t>>
    latent_truth_discovery(vector<vector<uint64_t>> &obs, uint iter, ABYParty *pt, e_role role);

    vector<vector<uint64_t>> sphere_kmeans(vector<vector<uint64_t>> &points, uint iter, ABYParty *pt, e_role role);

    vector<vector<int>>
    ttruth(vector<vector<vector<vector<uint64_t>>>> &all_kvec, vector<vector<vector<uint64_t>>> &answers, uint topK, ABYParty *pt,
           e_role role);

}

#endif //TTRUTH_MPC_TTRUTH_H
