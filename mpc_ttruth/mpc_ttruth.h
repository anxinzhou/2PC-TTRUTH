//
// Created by anxin on 8/8/20.
//

#ifndef TTRUTH_MPC_TTRUTH_H
#define TTRUTH_MPC_TTRUTH_H

#include "common/mpc_util.h"
extern const uint CLUSTER_NUM = 10;

namespace MPC {
    vector<vector<uint64_t>> latent_truth_discovery(vector<vector<uint64_t>> &obs, uint iter);
}


#endif //TTRUTH_MPC_TTRUTH_H
