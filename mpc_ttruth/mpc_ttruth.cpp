//
// Created by anxin on 8/8/20.
//

#include "mpc_ttruth.h"

namespace MPC {
    vector<vector<uint64_t>> latent_truth_discovery(vector<vector<vector<uint64_t>>> &all_obs, uint iter,ABYParty *pt, e_role role) {
        // truth label
        uint question_num = all_obs.size();
        uint user_num = all_obs[0].size();
        vector<vector<uint64_t>> tls(question_num, vector<uint64_t>(user_num,0));
        // random initialization truth label
        for(int i=0; i<question_num; i++) {
            for(int j=0;j<user_num;j++) {
                uint r = random(pt,role);
                uint threshold = role==SERVER? RAND_MAX/2:0;
                uint t = compare(r,threshold,pt,role);
                tls[i][j] = t;
            }
        }

        // initialize posterior counts
        vector<vector<uint64_t>> pos_counts(user_num,vector<uint64_t>(4,0));
        for(int j=0; j<user_num; j++) {
            // n,0,0
            for(int i=0; i<question_num; i++) {
                auto &tl = tls[i];
                // n, 0,0
                auto &ob = all_obs[i][j];
                vector<uint64_t>one_minus_tl(tl);
                for(int k=0; k<tl.size();k++) {
                    one_minus_tl[k] = -one_minus_tl[k];
                    if(role == SERVER) {
                        one_minus_tl[k] +=1;
                    }
                }

                vector<uint64_t>one_minus_ob(ob);
                for(int k=0;k<ob.size();k++) {
                    one_minus_ob[k] = -one_minus_ob[k];
                    if(role==SERVER) {
                        one_minus_ob[k] +=1;
                    }
                }

                pos_counts[j][0]+=inner_product(one_minus_tl, one_minus_ob, pt, role);
                // n,0,1
                pos_counts[j][1]+=inner_product(one_minus_tl, ob, pt, role);
                // n,1,0
                pos_counts[j][2]+=inner_product(tl, one_minus_ob, pt, role);
                // n,1,1
                pos_counts[j][3]+=inner_product(tl, ob, pt, role);
            }
        }

        // iteratively update truth label and prior count
        for(int i=0; i<question_num; i++) {
            auto& tl = tls[i];
            for(int k=0;k<tl.size();k++) {
                for (int j = 0; j < user_num; j++) {
                    // load paramers

                }
            }
        }
    }
}