//
// Created by anxin on 8/8/20.
//

#ifndef TTRUTH_TTRUTH_H
#define TTRUTH_TTRUTH_H

#include <vector>
#include <cmath>
using namespace std;

extern const int CLUSTER_NUM;

vector<vector<int>> latent_truth_discovery(vector<vector<vector<int>>> &all_obs,  uint iter);

#endif //TTRUTH_TTRUTH_H
