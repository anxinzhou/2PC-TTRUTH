//
// Created by anxin on 8/8/20.
//

#include "ttruth.h"

const uint ALPHA[4] = {90, 10, 50, 50};

const int CLUSTER_NUM = 10;

vector<vector<int>> latent_truth_discovery(vector<vector<vector<int>>> &all_obs, int iter) {
    uint question_num = all_obs.size();
    uint user_num = all_obs[0].size();
    vector<vector<int>> tls(question_num, vector<int>(user_num, 0));
    // random initialization truth label
    for (int i = 0; i < question_num; i++) {
        for (int j = 0; j < user_num; j++) {
            float p = double(rand()) / RAND_MAX;
            tls[i][j] = p > 0.5 ? 1 : 0;
        }
    }

    // initialize posterior counts
    vector<vector<int>> pos_counts(user_num, vector<int>(4, 0));
    for (int j = 0; j < user_num; j++) {
        // n,0,0
        for (int i = 0; i < question_num; i++) {
            auto &tl = tls[i];
            auto &ob = all_obs[i][j];

            for (int k = 0; k < tl.size(); k++) {
                int t = tl[k];
                int o = ob[k];
                pos_counts[j][t * 2 + o] += 1;
            }
        }
    }

    // iteratively update truth label and prior count
    for (int i = 0; i < question_num; i++) {
        auto &tl = tls[i];
        for (int k = 0; k < tl.size(); k++) {
            double p_t = 0;
            double p_negt = 0;
            int t = tl[k];
            for (int j = 0; j < user_num; j++) {
                auto &ob = all_obs[i][j];
                int o = ob[k];
                p_t += log2(pos_counts[j][t * 2 + o] - 1 + ALPHA[t * 2 + o]);
                p_t -= log2(pos_counts[j][t * 2] + pos_counts[j][t * 2 + 1] - 1 + ALPHA[t * 2] + ALPHA[t * 2 + 1]);
                p_negt += log2(pos_counts[j][(1 - t) * 2 + o] + ALPHA[(1 - t) * 2 + o]);
                p_negt -= log2(pos_counts[j][(1 - t) * 2] + pos_counts[j][(1 - t) * 2 + 1] + ALPHA[(1 - t) * 2]
                               + ALPHA[(1 - t) * 2 + 1]);
            }

            // update statistics
            double threshold = 1.0 / (1 + exp(p_t - p_negt));
            double p = double(rand()) / RAND_MAX;
            if (p > threshold) {
                for (int j = 0; j < user_num; j++) {
                    auto &ob = all_obs[i][j];
                    int o = ob[k];
                    pos_counts[j][t * 2 + o] -= 1;
                    pos_counts[j][(1 - t) * 2 + o] += 1;
                }
                tl[k] = 1 - t;
            }
        }
    }
    return tls;
}