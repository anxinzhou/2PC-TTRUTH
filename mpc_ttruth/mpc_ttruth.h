//
// Created by anxin on 8/8/20.
//

#ifndef TTRUTH_MPC_TTRUTH_H
#define TTRUTH_MPC_TTRUTH_H

#include "common/mpc_util.h"
#include <boost/serialization/vector.hpp>
//#include <boost/serialization/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <filesystem>

namespace MPC {
    extern const uint CLUSTER_NUM;
    vector<vector<uint64_t>> load_vector_vector(const string &path);
    void cache_vector_vector(vector<vector<uint64_t>> &cluster_index,const string &path);
    vector<vector<uint64_t>>
    latent_truth_discovery(vector<vector<uint64_t>> &obs, uint iter, ABYParty *pt, e_role role);

    vector<vector<uint64_t>> sphere_kmeans(vector<vector<uint64_t>> &points, uint iter, ABYParty *pt, e_role role);

    vector<vector<int>>
    ttruth(vector<vector<vector<vector<uint64_t>>>> &all_kvec, vector<vector<vector<uint64_t>>> &answers, uint topK, ABYParty *pt,
           e_role role);

}

#endif //TTRUTH_MPC_TTRUTH_H
