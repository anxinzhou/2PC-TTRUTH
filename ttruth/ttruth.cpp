//
// Created by anxin on 8/8/20.
//

#include "ttruth.h"

const uint ALPHA[4] = {90, 10, 50, 50};

const int CLUSTER_NUM = 10;
const uint SKM_ITER = 10;
const uint LTM_ITER = 10;

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

template <class T>
T inner_product(vector<T> &a, vector<T> &b) {
    int dim = a.size();
    double res = 0;
    for (int i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}

double distance(vector<double> &a, vector<double> &b) {
    return 2 - inner_product(a, b);
}

vector<vector<double>> cluster_init(vector<vector<double>> &points) {
    int cluster_num = CLUSTER_NUM;
    vector<vector<double>> cluster_centers(cluster_num);
    int dim = points[0].size();
    // pick first cluster center
    int r = rand() % cluster_num;
    cluster_centers[0] = points[r];

    vector<double> D(points.size(), INT_MAX);
    //pick other clusters
    for (int i = 1; i < cluster_num; i++) {
        //update min distance
        double total_dis = 0;
        for (int j = 0; j < points.size(); j++) {
            D[j] = min(D[j], distance(cluster_centers[i - 1], points[j]));
            total_dis += D[j];
        }

        // sample next cluster
        double p = double(rand()) / RAND_MAX;
        p = p * total_dis;

        vector<double> bar_D(points.size(), INT_MAX);
        for (int j = 0; j < points.size(); j++) {
            bar_D[j] = j == 0 ? D[0] : bar_D[j - 1];
        }

        int cluster_index = 0;
        for (int j = 0; j < points.size(); j++) {
            if (p < bar_D[j]) {
                cluster_index = j;
                break;
            }
        }
        cluster_centers[i] = points[cluster_index];
    }
    return cluster_centers;
}

vector<vector<int>> sphere_kmeans(vector<vector<double>> points, uint iter) {
    int cluster_num = CLUSTER_NUM;
    int dim = points[0].size();
    vector<vector<int>> cluster_index(points.size(), vector<int>(cluster_num, 0));
    // initialization
    auto cluster_centers = cluster_init(points);
    // iteration
    for (int t = 0; t < iter; t++) {
        // assign cluster index
        for (int j = 0; j < points.size(); j++) {
            double min_dis = -INT_MAX;
            int min_index = -1;
            for (int i = 0; i < cluster_num; i++) {
                double dis = distance(points[j], cluster_centers[i]);
                if (dis < min_dis) {
                    min_dis = dis;
                    min_index = i;
                }
            }
            cluster_index[j][min_index] = 1;
        }

        // update new cluster center
        vector<vector<double>> new_cluster_centers(cluster_num, vector<double>(dim, 0));
        for (int j = 0; j < points.size(); j++) {
            int l = -1;
            for (int i = 0; i < cluster_num; i++) {
                if (cluster_index[j][i] == 1) {
                    l = i;
                    break;
                }
            }
            for (int k = 0; k < dim; k++) {
                new_cluster_centers[l][k] += points[j][k];
            }
        }

        for (int i = 0; i < cluster_num; i++) {
            double val = inner_product(cluster_centers[i], cluster_centers[i]);
            if (val == 0) val = INT_MAX;
            val = 1 / (sqrt(val));
            for (int j = 0; j < dim; j++) {
                cluster_centers[i][j] *= val;
            }
        }
    }
}

vector<vector<int>> ttruth(vector<vector<vector<double>>> &all_kvec, int topK) {
    int question_num = all_kvec.size();
    int user_num = all_kvec[0].size();
    int cluster_number = CLUSTER_NUM;

    vector<vector<vector<int>>> all_obs(question_num,
                                        vector<vector<int>>(user_num,
                                                            vector<int>(cluster_number, 0)));

    //cluster keywords
    vector<vector<vector<vector<int>>>> all_cls(question_num,
                                                vector<vector<vector<int>>>(user_num));


    for (int i = 0; i < question_num; i++) {
        auto &obs = all_obs[i];
        auto &kvec = all_kvec[i];
        auto &cls = all_cls[i];
        vector<uint> key_num(user_num);
        vector<vector<double>> points;
        for (int j = 0; j < user_num; j++) {
            key_num[j] = kvec[j].size();
            points.insert(points.end(), kvec[j].begin(), kvec[j].end());
        }
        auto cluster_index = sphere_kmeans(points, SKM_ITER);

        int c = 0;
        int u = 0;
        for (int j = 0; j < cluster_index.size(); j++) {
            while (c + key_num[u] < j + 1) {
                c += key_num[u];
                u = u + 1;
            }
            for (int k = 0; k < cluster_number; k++) {
                obs[u][k] += cluster_index[j][k];
            }
            cls[u].push_back(std::move(cluster_index[j]));
        }

        // normalize observation
        for (int j = 0; j < user_num; j++) {
            auto &ob = obs[j];
            for (int k = 0; k < cluster_number; k++) {
                ob[k] = ob[k]>0? 1:0;
            }
        }
    }

    //latent truth discovery
    auto tls = latent_truth_discovery(all_obs, LTM_ITER);

    // calculate score of each user for each question
    vector<vector<int>>topk_index(question_num,vector<int>(topK,0));
    vector<vector<int>> score(question_num, vector<int>(user_num, 0));
    for (int i = 0; i < question_num; i++) {
        vector<pair<int,int>>score_and_index(user_num);
        auto &cls = all_cls[i];
        auto &tl = tls[i];

        for (int j = 0; j < user_num; j++) {
            auto &cl = cls[j];
            uint64_t s = 0;
            for (int k = 0; k < cl.size(); k++) {
                s += inner_product(tl, cl[k]);
            }
            score[i][j] = s;
            score_and_index.emplace_back(s,j);
        }

        sort(score_and_index.begin(), score_and_index.end(),greater<pair<int,int>>());
        for(int j=0; j<topK; j++) {
            topk_index[i][j] = score_and_index[j].second;
        }
    }

    return topk_index;
}