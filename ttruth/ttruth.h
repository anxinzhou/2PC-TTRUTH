//
// Created by anxin on 8/8/20.
//

#ifndef TTRUTH_TTRUTH_H
#define TTRUTH_TTRUTH_H

#include <vector>
#include <cmath>
#include<bits/stdc++.h>
using namespace std;

vector<vector<int>> latent_truth_discovery(vector<vector<vector<int>>> &all_obs,  uint iter);
vector<vector<int>> sphere_kmeans(vector<vector<double>>points, uint iter);
// return top-k score
vector<vector<int>> ttruth(vector<vector<vector<vector<double>>>> &all_kvec, int topK);

#endif //TTRUTH_TTRUTH_H
