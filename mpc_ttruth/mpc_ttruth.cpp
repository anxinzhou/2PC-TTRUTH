//
// Created by anxin on 8/8/20.
//

#include "mpc_ttruth.h"

const uint ALPHA[4] = {90, 10, 50, 50};
const uint CLUSTER_NUM = 10;

namespace MPC {
    vector<vector<uint64_t>>
    latent_truth_discovery(vector<vector<vector<uint64_t>>> &all_obs, uint iter,
                           ABYParty *pt, e_role role) {
        // truth label
        uint question_num = all_obs.size();
        uint user_num = all_obs[0].size();
        vector<vector<uint64_t>> tls(question_num, vector<uint64_t>(user_num, 0));
        // random initialization truth label
        for (int i = 0; i < question_num; i++) {
            for (int j = 0; j < user_num; j++) {
                uint r = random(pt, role);
                uint t = share_gt_const(r, RAND_MAX / 2, pt, role);
                tls[i][j] = t;
            }
        }

        // initialize posterior counts
        vector<vector<uint64_t>> pos_counts(user_num, vector<uint64_t>(4, 0));
        for (int j = 0; j < user_num; j++) {
            // n,0,0
            for (int i = 0; i < question_num; i++) {
                auto &tl = tls[i];
                auto &ob = all_obs[i][j];
                vector<uint64_t> one_minus_tl(tl);
                for (int k = 0; k < tl.size(); k++) {
                    one_minus_tl[k] = -one_minus_tl[k];
                    if (role == SERVER) {
                        one_minus_tl[k] += 1;
                    }
                }

                vector<uint64_t> one_minus_ob(ob);
                for (int k = 0; k < ob.size(); k++) {
                    one_minus_ob[k] = -one_minus_ob[k];
                    if (role == SERVER) {
                        one_minus_ob[k] += 1;
                    }
                }
                // n, 0,0
                pos_counts[j][0] += inner_product(one_minus_tl, one_minus_ob, pt, role);
                // n,0,1
                pos_counts[j][1] += inner_product(one_minus_tl, ob, pt, role);
                // n,1,0
                pos_counts[j][2] += inner_product(tl, one_minus_ob, pt, role);
                // n,1,1
                pos_counts[j][3] += inner_product(tl, ob, pt, role);
            }
        }

        // iteratively update truth label and prior count
        for (int i = 0; i < question_num; i++) {
            auto &tl = tls[i];
            for (int k = 0; k < tl.size(); k++) {
                uint p_t = 0;
                uint p_negt = 0;
                auto tl_eq_1 = eq(tl[k], role == SERVER ? 1 : 0, pt, role);
                uint64_t tl_eq_0 = role == SERVER ? -tl_eq_1 : -tl_eq_1 + 1;
                vector<vector<uint64_t>> tmp_prior(user_num,vector<uint64_t>(2,0));
                for (int j = 0; j < user_num; j++) {
                    // load paramers
                    auto &ob = all_obs[i][j];
                    auto ob_eq_1 = eq(ob[k], role == SERVER ? 1 : 0, pt, role);
                    uint64_t ob_eq_0 = role == SERVER ? -ob_eq_1 : -ob_eq_1 + 1;

                    uint64_t n_u_t_0 = product(tl_eq_1, pos_counts[j][2], pt, role) +
                                       product(tl_eq_0, pos_counts[j][0], pt, role);
                    uint64_t n_u_t_1 = product(tl_eq_1, pos_counts[j][3], pt, role) +
                                       product(tl_eq_0, pos_counts[j][1], pt, role);
                    uint64_t n_u_negt_0 = product(tl_eq_0, pos_counts[j][2], pt, role) +
                                          product(tl_eq_1, pos_counts[j][0], pt, role);
                    uint64_t n_u_negt_1 = product(tl_eq_0, pos_counts[j][3], pt, role) +
                                          product(tl_eq_1, pos_counts[j][1], pt, role);

                    uint64_t n_u_t_o = product(ob_eq_1, n_u_t_1, pt, role) +
                                       product(ob_eq_0, n_u_t_0, pt, role);
                    uint64_t n_u_negt_o = product(ob_eq_0, n_u_t_1, pt, role) +
                                          product(ob_eq_1, n_u_t_0, pt, role);

                    uint64_t alpha_t_0 = product(tl_eq_1, role == SERVER ? ALPHA[2] : 0, pt, role) +
                                         product(tl_eq_0, role == SERVER ? ALPHA[0] : 0, pt, role);
                    uint64_t alpha_t_1 = product(tl_eq_1, role == SERVER ? ALPHA[3] : 0, pt, role) +
                                         product(tl_eq_0, role == SERVER ? ALPHA[1] : 0, pt, role);
                    uint64_t alpha_negt_0 = product(tl_eq_0, role == SERVER ? ALPHA[2] : 0, pt, role) +
                                            product(tl_eq_1, role == SERVER ? ALPHA[0] : 0, pt, role);
                    uint64_t alpha_negt_1 = product(tl_eq_0, role == SERVER ? ALPHA[3] : 0, pt, role) +
                                            product(tl_eq_1, role == SERVER ? ALPHA[1] : 0, pt, role);
                    uint64_t alpha_t_o = product(ob_eq_1, alpha_t_1, pt, role) +
                                         product(ob_eq_0, alpha_t_0, pt, role);
                    uint64_t alpha_negt_o = product(ob_eq_0, alpha_t_1, pt, role) +
                                            product(ob_eq_1, alpha_t_0, pt, role);
                    // update p_t
                    uint tmp1 = n_u_t_o + alpha_t_o;
                    if(role==SERVER) {
                        tmp1-=1;
                    }
                    tmp1 = log(tmp1 << FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                    uint tmp2 = n_u_t_1 + n_u_t_0 + alpha_t_0 + alpha_t_1;
                    if(role==SERVER) {
                        tmp2-=1;
                    }
                    tmp2 = log(tmp2<<FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                    p_t += tmp1 - tmp2;
                    // update p_(1-t)
                    uint tmp3 = n_u_negt_o + alpha_negt_o;
                    tmp3 = log(tmp3<<FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                    uint tmp4 = n_u_negt_0 + n_u_negt_1 + alpha_negt_0 + alpha_negt_1;
                    tmp4 = log(tmp4<<FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                    p_negt = tmp3 - tmp4;

                    tmp_prior[j][0] = n_u_t_o;
                    tmp_prior[j][1] = n_u_negt_o;
                }
                uint64_t threshold_p = sigmoid(p_negt-p_t,FLOAT_SCALE_FACTOR,FLOAT_SCALE_FACTOR,pt,role);
                uint64_t r = random(pt, role);
                threshold_p *= RAND_MAX;
                r = r<<FLOAT_SCALE_FACTOR;
                uint64_t flip = gt(threshold_p, r, pt, role);
                // update statistics
                for(int j=0;j<user_num;j++) {
                    auto &ob = all_obs[i][j];
                    auto ob_eq_1 = eq(ob[k], role == SERVER ? 1 : 0, pt, role);
                    uint64_t ob_eq_0 = role == SERVER ? -ob_eq_1 : -ob_eq_1 + 1;

                    tmp_prior[j][0] -= flip;
                    tmp_prior[j][1] += flip;
                    uint64_t  n_u_t_o = tmp_prior[j][0];
                    uint64_t  n_u_negt_o = tmp_prior[j][1];
                    pos_counts[j][0] = product(ob_eq_0,
                                               product(tl_eq_0, n_u_t_o, pt, role) +
                                               product(tl_eq_1, n_u_negt_o,pt,role),pt,role);

                    pos_counts[j][1] = product(ob_eq_1,
                                               product(tl_eq_0, n_u_t_o, pt, role) +
                                               product(tl_eq_1, n_u_negt_o,pt,role),pt,role);
                    pos_counts[j][2] = product(ob_eq_0,
                                               product(tl_eq_1, n_u_t_o, pt, role) +
                                               product(tl_eq_0, n_u_negt_o,pt,role),pt,role);
                    pos_counts[j][0] = product(ob_eq_1,
                                               product(tl_eq_1, n_u_t_o, pt, role) +
                                               product(tl_eq_0, n_u_negt_o,pt,role),pt,role);
                }
                // update truth label
                tl[k] = product(flip, role==SERVER? 1-tl[k]:-tl[k], pt, role) +
                        product(role==SERVER? 1-flip:-flip, tl[k],pt,role);
            }
            return tls;
        }
    }

    // only for float number
    uint64_t distance(vector<uint64_t>&a, vector<uint64_t>&b, ABYParty *pt, e_role role) {
        uint64_t gamma = role==SERVER?2:0;
        uint64_t  tmp = inner_product(a,b,pt,role);
        tmp = right_shift_const(tmp,FLOAT_SCALE_FACTOR,pt,role);
        return gamma - tmp;
    }

    //only for float number
    // sphere kmeans++ for init
    vector<vector<uint64_t>> cluster_init(vector<vector<uint64_t>> &points, ABYParty *pt, e_role role) {
        vector<vector<uint64_t>> cluster_centers(CLUSTER_NUM);
        uint k = CLUSTER_NUM;
        uint dim = points[0].size();
        // pick first cluster
        uint r = random(pt, role);
        r = right_shift(k*r,UINT_LEN,pt,role);
        for(int i=0;i<points.size();i++) {
            uint64_t cmp = share_eq_const(r,i,pt,role);
            vector<uint64_t> tmp(dim, cmp);
            tmp = product(tmp,points[i],pt,role);
            for(int j=0; j<dim;j++) {
                cluster_centers[0][j] += tmp[j];
            }
        }
        // right shift
        cluster_centers[0] = right_shift_const(cluster_centers[0],FLOAT_SCALE_FACTOR,pt,role);

        vector<uint64_t> D(points.size(),UINT_MAX);
        //pick other clusters

        for(int i=1;i<k;i++) {
            for(int j=0;j<points.size();j++) {
                D[j] = min(D[j],
                           distance(points[j],cluster_centers[i-1],pt,role),pt,role);
            }
            // pick new cluster
            vector<uint64_t> bar_D(points.size());
            for(int j=0;j<points.size();j++) {
                bar_D[j] = j==0? D[0]: bar_D[j-1] + D[j];
            }
            uint64_t J = bar_D[points.size()-1];

            uint r = random(pt, role);
            uint64_t  p = product(r, J, pt, role);
            p = right_shift_const(p, UINT_LEN, pt, role);


        }
    }

    vector<uint64_t> sphere_kmeans(vector<vector<uint64_t>> &point, uint iter, ABYParty *pt, e_role role) {

    }
}