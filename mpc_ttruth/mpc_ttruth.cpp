//
// Created by anxin on 8/8/20.
//

#include <ABY/extern/ENCRYPTO_utils/src/ENCRYPTO_utils/crypto/crypto.h>
#include "mpc_ttruth.h"

namespace MPC {
    const uint ALPHA[4] = {80, 20, 40, 60};
    const uint CLUSTER_NUM = 10;
    const uint SKM_ITER = 1;
    const uint LTM_ITER = 1;
    const uint ANSWER_LEN = 128;
    const uint64_t NEG_PAD_VAL = uint64_t(1)<<32;

    void print_share(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        auto tmp = open_share(a,pt,role);
        for(int i=0; i<a.size();i++) {
            cout<<tmp[i]<<" ";
        }
        cout<<endl;
    }

    void print_share(uint64_t a, ABYParty *pt, e_role role) {
        auto tmp = open_share(a,pt,role);
        cout<<tmp<<endl;
    }

    void print_scaled_share(uint64_t a, ABYParty *pt, e_role role) {
        auto tmp = open_share(a,pt,role);
        double val;
        if(tmp>INT_MAX) val = 0-double(-tmp);
        else val = tmp;
        cout<<val/(1<<FLOAT_SCALE_FACTOR)<<endl;
    }

    void print_scaled_share(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        auto tmp = open_share(a,pt,role);
        for(int i=0; i<a.size();i++) {
            double val;
            if(tmp[i]>INT_MAX) val = 0-double(-tmp[i]);
            else val = tmp[i];
            cout<<val/(1<<FLOAT_SCALE_FACTOR)<<" ";
        }
        cout<<endl;
    }

    void print_distance(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        auto tmp = open_share(a,pt,role);
        for(int i=0; i<a.size();i++) {
            double val;
            if(tmp[i]>uint64_t(INT_MAX)<<FLOAT_SCALE_FACTOR) val = 0-double(-tmp[i]);
            else val = tmp[i];

            cout<<val/(uint64_t(1)<<(FLOAT_SCALE_FACTOR)*2)<<" ";
        }
        cout<<endl;
    }

    void print_distance(uint64_t a, ABYParty *pt, e_role role) {
        auto tmp = open_share(a,pt,role);

        double val = tmp>>FLOAT_SCALE_FACTOR;
        cout<<val/(1<<FLOAT_SCALE_FACTOR)<<endl;
    }

    vector<vector<uint64_t>>
    latent_truth_discovery(vector<vector<vector<uint64_t>>> &all_obs, uint iter,
                           ABYParty *pt, e_role role) {
        // truth label
        uint question_num = all_obs.size();
        uint user_num = all_obs[0].size();
        uint cluster_num = CLUSTER_NUM;
        vector<vector<uint64_t>> tls(question_num, vector<uint64_t>(cluster_num, 0));
        // random initialization truth label
        cout<<"initial truth label"<<endl;
        auto start = clock();
        for (int i = 0; i < question_num; i++) {
            for (int j = 0; j < cluster_num; j++) {
                uint64_t r = random(pt, role);
//                cout<<"------------------------"<<endl;
//                print_share(r,pt,role);
//                cout<< (1<<RANDOMNESS_BIT) / 2 <<endl;
                uint64_t t = share_gt_const(r, (1<<RANDOMNESS_BIT) / 2, pt, role);
                tls[i][j] = t;
            }
//            print_share(tls[i],pt,role);
//            exit(-1);
        }
        auto end = clock();
//        cout<<"Init truth label time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

        // initialize posterior counts
        start = clock();
        vector<vector<uint64_t>> pos_counts(user_num, vector<uint64_t>(4, 0));
        for (int j = 0; j < user_num; j++) {
            // n,0,0
            for (int i = 0; i < question_num; i++) {
                auto &tl = tls[i];
                auto &ob = all_obs[i][j];
                vector<uint64_t> one_minus_tl(tl);
                for (int k = 0; k < cluster_num; k++) {
                    one_minus_tl[k] = -one_minus_tl[k];
                    if (role == SERVER) {
                        one_minus_tl[k] += 1;
                    }
                }

                vector<uint64_t> one_minus_ob(ob);
                for (int k = 0; k < cluster_num; k++) {
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

        end = clock();
//        cout<<"Init posterior time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

//        cout<<"initial posterior counts"<<endl;
//        for(int j=0; j<user_num; j++) {
//            print_share(pos_counts[j],pt,role);
//        }

        // iteratively update truth label and prior count
        cout<<"iteratively update truth label and prior count"<<endl;
        start = clock();
        for (int t = 0; t < iter; t++) {
            for (int i = 0; i < question_num; i++) {

                string pos_count_path = string("../pos_count_cache/")
                        +"iter"+to_string(t)+"_q"+to_string(i)+"_role"+to_string(role);
                string truth_label_path = string("../truth_label_cache/")
                                          +"iter"+to_string(t)+"_q"+to_string(i)+"_role"+to_string(role);

                if(false and filesystem::exists(pos_count_path) && filesystem::exists(truth_label_path)) {
                    cout<<"iter"<<t<<" question "<<i<<"-th"<<" truth statistics exists"<<endl;
                    pos_counts = load_vector_vector(pos_count_path);
                    tls = load_vector_vector(truth_label_path);
                } else {
                    cout<<"update truth label for "<<i<<"-th question"<<endl;
//                cout<<pt->GetTotalGates()<<endl;
//                cout<<pt->GetTotalDepth()<<endl;
                    auto &tl = tls[i];
                    for (int k = 0; k < tl.size(); k++) {
                        uint64_t p_t = 0;
                        uint64_t p_negt = 0;
                        auto tl_eq_1 = eq(tl[k], role == SERVER ? 1 : 0, pt, role);
                        uint64_t tl_eq_0 = role == SERVER ? -tl_eq_1 : -tl_eq_1 + 1;
                        vector<vector<uint64_t>> tmp_prior(user_num, vector<uint64_t>(2, 0));

                        auto start = clock();
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
                            uint64_t n_u_negt_o = product(ob_eq_1, n_u_negt_1, pt, role) +
                                                  product(ob_eq_0, n_u_negt_0, pt, role);

                            //                        cout<<"-----------------------------"<<endl;
//                        print_share(tl[k],pt,role);
//                        print_share(tl_eq_1,pt,role);
//                        print_share(ob[k], pt, role);
//                        print_share(ob_eq_1,pt,role);
//                        print_share(pos_counts[j],pt,role);
//                        print_share(n_u_t_o, pt, role);
//                        print_share(n_u_negt_o,pt,role);
//                        print_share(n_u_t_0,pt,role);
//                        print_share(n_u_t_1,pt,role);
//                        print_share(n_u_negt_0,pt,role);
//                        print_share(n_u_negt_1,pt,role);
//                        exit(-1);

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
                            uint64_t alpha_negt_o = product(ob_eq_1, alpha_negt_1, pt, role) +
                                                    product(ob_eq_0, alpha_negt_0, pt, role);
                            // update p_t
                            uint64_t tmp1 = n_u_t_o + alpha_t_o;
                            if (role == SERVER) {
                                tmp1 -= 1;
                            }
//                        print_share(tmp1,pt,role);


                            tmp1 = left_shift_const(tmp1, FLOAT_SCALE_FACTOR, pt, role);

                            auto start = clock();
                            tmp1 = log(tmp1, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                            auto end = clock();
//                        cout<<"Log time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//                        print_scaled_share(tmp1,pt,role);
//                        exit(-1);
                            uint64_t tmp2 = n_u_t_1 + n_u_t_0 + alpha_t_0 + alpha_t_1;
                            if (role == SERVER) {
                                tmp2 -= 1;
                            }
                            tmp2= left_shift_const(tmp2, FLOAT_SCALE_FACTOR, pt, role);
                            tmp2 = log(tmp2, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                            p_t += tmp1 - tmp2;
                            // update p_(1-t)
                            uint64_t tmp3 = n_u_negt_o + alpha_negt_o;
                            tmp3 = left_shift_const(tmp3, FLOAT_SCALE_FACTOR, pt, role);
                            tmp3 = log(tmp3, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                            uint64_t tmp4 = n_u_negt_0 + n_u_negt_1 + alpha_negt_0 + alpha_negt_1;
                            tmp4 = left_shift_const(tmp4, FLOAT_SCALE_FACTOR, pt, role);
                            tmp4 = log(tmp4, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                            p_negt += tmp3 - tmp4;

//                        print_scaled_share(p_t,pt,role);
//                        print_scaled_share(p_negt,pt,role);
//                        exit(-1);

                            tmp_prior[j][0] = n_u_t_o;
                            tmp_prior[j][1] = n_u_negt_o;
                        }

                        end = clock();
//                    cout<<"Calculate prob time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

                        start = clock();

                        uint64_t threshold_p = sigmoid(p_negt - p_t, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
                        end = clock();
//                    cout<<"Sigmoid time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//                    cout<<"-----------------------------"<<endl;
//                    print_scaled_share(p_t,pt,role);
//                    print_scaled_share(p_negt,pt,role);
//                    print_scaled_share(threshold_p,pt,role);
//                    exit(-1);

                        uint64_t r = random(pt, role);

                        r = left_shift_const(r, FLOAT_SCALE_FACTOR, pt, role);
//                    print_scaled_share(threshold_p,pt,role);
                        threshold_p = left_shift_const(threshold_p, RANDOMNESS_BIT, pt, role);
                        uint64_t flip = gt(threshold_p, r, pt, role);
//                    cout<<"------------------"<<endl;
//                    print_share(r,pt,role);
//                    print_share(threshold_p,pt,role);
//                    cout<<"flip"<<endl;
//                    print_share(flip, pt, role);
//                    print_share(r,pt,role);
//                    print_share(threshold_p,pt,role);
//                    cout<<"---------------------"<<endl;
                        // update statistics
                        for (int j = 0; j < user_num; j++) {
                            auto &ob = all_obs[i][j];
                            auto ob_eq_1 = eq(ob[k], role == SERVER ? 1 : 0, pt, role);
                            uint64_t ob_eq_0 = role == SERVER ? -ob_eq_1 : -ob_eq_1 + 1;

                            tmp_prior[j][0] -= flip;
                            tmp_prior[j][1] += flip;
                            uint64_t n_u_t_o = tmp_prior[j][0];
                            uint64_t n_u_negt_o = tmp_prior[j][1];

//                        cout<<"---------------------"<<endl;
//                        print_share(tl[k],pt,role);
//                        print_share(ob[k],pt,role);
//                        print_share(tl_eq_0, pt, role);
//                        print_share(tl_eq_1, pt, role);
//                        print_share(ob_eq_0, pt, role);
//                        print_share(ob_eq_1,pt,role);
//                        print_share(pos_counts[j],pt,role);
//                        print_share(n_u_t_o,pt,role);
//                        print_share(n_u_negt_o,pt,role);
//                        exit(-1);

                            pos_counts[j][0] = product(ob_eq_0,
                                                       product(tl_eq_0, n_u_t_o, pt, role) +
                                                       product(tl_eq_1, n_u_negt_o, pt, role), pt, role)
                                               + product(ob_eq_1, pos_counts[j][0],pt,role);

                            pos_counts[j][1] = product(ob_eq_1,
                                                       product(tl_eq_0, n_u_t_o, pt, role) +
                                                       product(tl_eq_1, n_u_negt_o, pt, role), pt, role)
                                               + product(ob_eq_0, pos_counts[j][1],pt,role);

                            pos_counts[j][2] = product(ob_eq_0,
                                                       product(tl_eq_1, n_u_t_o, pt, role) +
                                                       product(tl_eq_0, n_u_negt_o, pt, role), pt, role)
                                               + product(ob_eq_1, pos_counts[j][2],pt,role);


                            pos_counts[j][3] = product(ob_eq_1,
                                                       product(tl_eq_1, n_u_t_o, pt, role) +
                                                       product(tl_eq_0, n_u_negt_o, pt, role), pt, role)
                                               + product(ob_eq_0, pos_counts[j][3],pt,role);
//                        exit(-1);
                        }
                        // update truth label
                        tl[k] = product(flip, role == SERVER ? 1 - tl[k] : -tl[k], pt, role) +
                                product(role == SERVER ? 1 - flip : -flip, tl[k], pt, role);
                    }
//                    cache_vector_vector(pos_counts, pos_count_path);
//                    cache_vector_vector(tls, truth_label_path);
                }

//                for(int j=0;j<pos_counts.size();j++) {
//                    if(i==1)
//                    print_share(pos_counts[j],pt,role);
//                }

            }

//            cout<<" truth label"<<endl;
//            print_share(tls[0],pt,role);

        }
        end = clock();
//        cout<<"Update truth label time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

//        cout<<"after latent truth posterior counts"<<endl;
//        for(int j=0; j<user_num; j++) {
//            print_share(pos_counts[j],pt,role);
//        }
        cout<<"final truth label"<<endl;
        for(int i=0; i<question_num; i++){
            print_share(tls[i],pt,role);
        }
        return tls;
    }

    // only for float number
    uint64_t distance(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty *pt, e_role role) {
        uint64_t gamma = role == SERVER ? 2<<FLOAT_SCALE_FACTOR : 0;
        auto start = clock();
        uint64_t tmp = inner_product(a, b, pt, role);
        auto end = clock();
        cout<<"Inner product Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        start = clock();
        tmp = right_shift_const(tmp, FLOAT_SCALE_FACTOR, pt, role);
        end = clock();
        cout<<"Right shift time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        return gamma - tmp;
    }

    // only for float number
    uint64_t non_right_shift_distance(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty *pt, e_role role) {
        uint64_t gamma = role == SERVER ? uint64_t (2)<<(FLOAT_SCALE_FACTOR+FLOAT_SCALE_FACTOR) : 0;
//        auto start = clock();

        uint64_t tmp = inner_product(a, b, pt, role);

//        print_scaled_share(a,pt,role);
//        cout<<"-----------------------------------"<<endl;

//        auto end = clock();
//        cout<<"Inner product Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        return gamma - tmp;
    }

    //only for float number
    // sphere kmeans++ for init
    vector<vector<uint64_t>> cluster_init(vector<vector<uint64_t>> &points, ABYParty *pt, e_role role) {
        uint cluster_num = CLUSTER_NUM;
        uint dim = points[0].size();
        vector<vector<uint64_t>> cluster_centers(CLUSTER_NUM,vector<uint64_t>(dim,0));
        // pick first cluster
        cout<<"number of point "<<points.size()<<endl;
        uint64_t r = random(pt, role);
        r = right_shift_const(r*points.size(), RANDOMNESS_BIT, pt, role);
//        print_share(r,pt,role);
//        auto start = clock();
        vector<uint64_t> vec_r(points.size(), r);
        vector<uint64_t> index(points.size());
        for (int i = 0; i < points.size(); i++) index[i] = role == SERVER ? i : 0;
        vector<uint64_t> cmp_vec = eq(vec_r, index, pt, role);

//        print_share(cmp_vec,pt,role);

        for (int i = 0; i < points.size(); i++) {
            uint64_t cmp = cmp_vec[i];
            vector<uint64_t> tmp(dim, cmp);
//            auto start = clock();
            tmp = product(tmp, points[i], pt, role);
//            auto end = clock();
//            cout<<"Vector product time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
            for (int j = 0; j < dim; j++) {
                cluster_centers[0][j] += tmp[j];
            }
        }
//        print_scaled_share(cluster_centers[0],pt,role);

//        auto end = clock();
//        cout<<"Pick first cluster Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

        vector<uint64_t> D(points.size(), uint64_t(10)<<FLOAT_SCALE_FACTOR<<FLOAT_SCALE_FACTOR);
        //pick other clusters
        for (int i = 1; i < cluster_num; i++) {
//            start = clock();
            // calculate new min distance
            vector<uint64_t> new_dis(points.size());
            for(int j=0; j<points.size(); j++) {
                new_dis[j] = non_right_shift_distance(points[j], cluster_centers[i - 1], pt, role);
            }
//            print_distance(new_dis,pt,role);
//            cout<<"-------------------"<<endl;
//            exit(-1);

            for (int j = 0; j < points.size(); j++) {
//                auto start = clock();
                D[j] = min(D[j],
                           non_right_shift_distance(points[j], cluster_centers[i - 1], pt, role), pt, role);
//                auto end = clock();
//                cout<<"Min Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
            }
//            print_distance(D,pt,role);
//            exit(-1);
            // pick new cluster
            vector<uint64_t> bar_D(points.size());
            for (int j = 0; j < points.size(); j++) {
                bar_D[j] = j == 0 ? D[0] : bar_D[j - 1] + D[j];
            }
            uint64_t J = bar_D[points.size() - 1];

            r = random(pt, role);
//            print_share(r,pt,role);
            uint64_t p = J;
            p = right_shift_const(p, RANDOMNESS_BIT, pt, role);
            p = product(r, p, pt, role);

            // batch comparision
            vector<uint64_t>p_vec(points.size(),p);
            cmp_vec = gt(bar_D, p_vec, pt,role);
//            print_distance(p_vec,pt,role);
//            print_distance(bar_D,pt,role);
//            print_share(cmp_vec,pt,role);
//            cout<<"-------------------------------------"<<endl;
            for (int j = 0; j < points.size(); j++) {
                if (j == points.size() - 1) {
                    for (int k = 0; k < dim; k++) {
                        cluster_centers[i][k] += points[j][k];
                    }
                    continue;
                }

//                uint64_t cmp = gt(bar_D[j], p, pt, role);
                uint64_t cmp = cmp_vec[j];
                vector<uint64_t> tmp(dim, cmp);
                auto diff = minus(points[j], points[j + 1]);
                tmp = product(tmp, diff, pt, role);
                for (int k = 0; k < dim; k++) {
                    cluster_centers[i][k] += tmp[k];
                }
            }
//            end = clock();
//            cout<<"Select a new cluster Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        }
        //test
//        for(int i=0; i<cluster_num; i++) {
//            print_scaled_share(cluster_centers[i],pt,role);
//            cout<<"------------------"<<endl;
////            exit(-1);
//        }
        return cluster_centers;
    }

    vector<vector<uint64_t>> sphere_kmeans(vector<vector<uint64_t>> &points, uint iter, ABYParty *pt, e_role role) {
        cout<<"convergence points "<<points.size()<<endl;
        cout<<"-----------------------------"<<endl;
//        for(int i=0; i<points.size();i++) {
//            print_scaled_share(points[i],pt,role);
//        }
//        exit(-1);
        cout<<"-----------------------------"<<endl;
        int cluster_num = CLUSTER_NUM;
        int dim = points[0].size();
        vector<vector<uint64_t>> cluster_index(points.size(), vector<uint64_t>(cluster_num, 0));
        // initialization
//        cout<<"start cluster init"<<endl;
//        cout<<"point size "<<
        auto start = clock();
        auto cluster_centers = cluster_init(points, pt, role);
        auto end = clock();
//        cout<<"Cluster init time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//        cout<<"cluster init finish!"<<endl;
        // iteration
        cout<<"start clustering"<<endl;
        for (int t = 0; t < iter; t++) {
            // assign cluster index
//            auto start = clock();

            for (int j = 0; j < points.size(); j++) {
                vector<uint64_t> dis(cluster_num);
                for (int i = 0; i < cluster_num; i++) {
                    dis[i] = non_right_shift_distance(points[j], cluster_centers[i], pt, role);
                }
                auto start = clock();
                cluster_index[j] = argmin_vector(dis, pt, role);
                auto end = clock();
//                cout<<"argmin time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

//                print_scaled_share(points[j],pt,role);
//                cout<<"--------------------------------"<<endl;
//                print_distance(dis,pt,role);
//                cout<<"--------------------------------"<<endl;
//                print_share(cluster_index[j],pt,role);
//                exit(-1);
            }
//            auto end = clock();
//            cout<<"assign cluster index Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;


            // update new cluster center
//            start = clock();
            vector<vector<uint64_t>> new_cluster_centers(cluster_num, vector<uint64_t>(dim, 0));
            for (int i = 0; i < cluster_num; i++) {
                for (int j = 0; j < points.size(); j++) {
                    auto l = cluster_index[j][i];
                    vector<uint64_t> tmp(dim, l);
//                    auto start = clock();
                    tmp = product(tmp, points[j], pt, role);
//                    auto end = clock();
//                    cout<<"product Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
                    for (int k = 0; k < dim; k++) {
                        new_cluster_centers[i][k] += tmp[k];
                    }
                }
            }

//            end = clock();
//            cout<<"update cluster center time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
            //normalize cluster center
            start = clock();
            for (int i = 0; i < cluster_num; i++) {
//                print_scaled_share(cluster_centers[i],pt,role);
                uint64_t val = inner_product(new_cluster_centers[i], new_cluster_centers[i], pt, role);
                val = right_shift_const(val, FLOAT_SCALE_FACTOR, pt, role);

//                auto start = clock();
//                cout<<"val"<<endl;
//                print_scaled_share(val,pt,role);
                val = rep_square_root(val, FLOAT_SCALE_FACTOR, FLOAT_SCALE_FACTOR, pt, role);
//                print_scaled_share(val,pt,role);
//                auto end = clock();
//                cout<<"req square root Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
                vector<uint64_t> tmp(dim, val);
                new_cluster_centers[i] = product(tmp, new_cluster_centers[i], pt, role);

//                print_distance(cluster_centers[i],pt,role);
//                start = clock();

                for(int j=0;j<dim;j++) {
                    new_cluster_centers[i][j] += role==SERVER? NEG_PAD_VAL << FLOAT_SCALE_FACTOR:0;
                    new_cluster_centers[i][j] = right_shift_const(new_cluster_centers[i][j],FLOAT_SCALE_FACTOR,pt,role);
                    new_cluster_centers[i][j] -= role==SERVER? NEG_PAD_VAL:0;
                }

//                print_scaled_share(cluster_centers[i],pt,role);
//                end = clock();
//                cout<<"batch right shift time "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//                exit(-1);
            }

//            print_scaled_share(cluster_centers[5],pt,role);
//            cout<<"---------------"<<endl;
            // test convergence
            double diff = 0;
            for(int i=0; i<cluster_num; i++) {
                auto newc = open_share(new_cluster_centers[i],pt,role);
                auto oldc = open_share(cluster_centers[i],pt,role);
//                if(t==SKM_ITER-1) {
//                    cout<<"---------------"<<endl;
//                    print_scaled_share(newc,pt,role);
//                    print_scaled_share(oldc,pt,role);
//                }
                for(int j=0; j<dim;j++) {
                    double x;
                    double y;
                    if(newc[j]>INT_MAX) x = 0-double(-newc[j]);
                    else x = newc[j];
                    if(oldc[j]>INT_MAX) y = 0-double(-oldc[j]);
                    else y = oldc[j];
                    x = x/(1<<FLOAT_SCALE_FACTOR);
                    y = y/(1<<FLOAT_SCALE_FACTOR);
//                    cout<<x<<" "<<y<<endl;
                    diff+= (x-y)*(x-y);
                }
            }
            cout<<t<<"-th round "<<"convergence: "<<diff<<endl;

//            cout<<"one iteration finish"<<endl;
            swap(cluster_centers, new_cluster_centers);
//            print_scaled_share(cluster_centers[5],pt,role);
//            cout<<"---------------"<<endl;

            end = clock();

//            cout<<"------------------------"<<endl;
//            for(int i=0;i<cluster_num;i++) {
//                auto val = inner_product(cluster_centers[i],cluster_centers[i],pt,role);
//                print_distance(val,pt,role);
//            }
//            cout<<"------------------------"<<endl;

            vector<uint64_t>statistics(cluster_num,0);
            for(int i=0;i<points.size();i++) {
                for(int j=0;j<cluster_num;j++) {
                    statistics[j]+=cluster_index[i][j];
                }
            }
            statistics = open_share(statistics,pt,role);
            for(int j=0;j<cluster_num;j++) {
                cout<<statistics[j]<<" ";
            }
            cout<<endl;
//            cout<<"Normalize time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        }

//        for(int i=0; i<cluster_num;i++) {
//            print_scaled_share(cluster_centers[i],pt,role);
//            cout<<"---------------------"<<endl;
//        }


        return cluster_index;
    }

    vector<vector<uint64_t>> load_vector_vector(const string &path) {
        vector<vector<uint64_t>> cluster_index;
        ifstream index_file(path);
        if (!index_file.is_open()) {
            cout << "can not open index file " << path<< endl;
            cout << strerror(errno) << endl;
            exit(0);
        }
        boost::archive::binary_iarchive ia(index_file);
        ia>>cluster_index;
        index_file.close();
        return cluster_index;
    }

    void cache_vector_vector(vector<vector<uint64_t>> &cluster_index,const string &path) {
        ofstream index_file(path);
        if (!index_file.is_open()) {
            cout << "can not open index file " << path<< endl;
            cout << strerror(errno) << endl;
            exit(0);
        }
        boost::archive::binary_oarchive oa(index_file);
        oa << cluster_index;
        index_file.close();
    }

    vector<vector<int>>
    ttruth(vector<vector<vector<vector<uint64_t>>>> &all_kvec, vector<vector<vector<uint64_t>>> &answers, uint topK, ABYParty *pt,
           e_role role) {
        uint question_num = all_kvec.size();
        uint user_num = all_kvec[0].size();
        uint cluster_number = CLUSTER_NUM;
        vector<vector<vector<uint64_t>>> all_obs(question_num,
                                                 vector<vector<uint64_t>>(user_num,
                                                                          vector<uint64_t>(cluster_number,
                                                                                           0)));
        // clustering keywords
        vector<vector<vector<vector<uint64_t>>>> all_cls(question_num,
                                                         vector<vector<vector<uint64_t>>>(user_num));

        auto start_all = clock();
        auto start = clock();
        for (int i = 0; i < question_num; i++) {
            auto &obs = all_obs[i];
            auto &kvec = all_kvec[i];
            auto &cls = all_cls[i];
            vector<uint> key_num(user_num);
            vector<vector<uint64_t>> points;
            for (int j = 0; j < user_num; j++) {
                key_num[j] = kvec[j].size();
                points.insert(points.end(), kvec[j].begin(), kvec[j].end());
            }

            vector<vector<uint64_t>> cluster_index;
            string index_path = string("../index_cache/q"+to_string(i)+"_role"+to_string(role));
//            cout<< filesystem::exists(index_path) <<endl;
            if(false and filesystem::exists(index_path)) {
                cout<<"question "<<i<<"-th"<<" cluster index exists"<<endl;
                cluster_index = load_vector_vector(index_path);

                vector<uint64_t>statistics(CLUSTER_NUM,0);
                for(int i=0;i<points.size();i++) {
                    for(int j=0;j<CLUSTER_NUM;j++) {
                        statistics[j]+=cluster_index[i][j];
                    }
                }
                statistics = open_share(statistics,pt,role);
                for(int j=0;j<CLUSTER_NUM;j++) {
                    cout<<statistics[j]<<" ";
                }
                cout<<endl;

            } else {
                cout<<"compute question "<<i<<"-th"<<" cluster index"<<endl;
                auto start=clock();
                cluster_index = sphere_kmeans(points, SKM_ITER, pt, role);
//                cache_vector_vector(cluster_index, index_path);
                auto end=clock();
//                cout<<"worker num "<<user_num<<" kmeans time "<<double(end-start)/CLOCKS_PER_SEC<<endl;
            }
            uint c = 0;
            uint u = 0;
            for (int j = 0; j < points.size(); j++) {
                while (c + key_num[u] < j + 1) {
                    c += key_num[u];
                    u = u + 1;
                }
                for (int k = 0; k < cluster_number; k++) {
                    obs[u][k] += cluster_index[j][k];
                }
                cls[u].push_back(std::move(cluster_index[j]));
            }
        }
        auto end = clock();
        cout<<"Cluster time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

        start = clock();
        for(int i=0;i<question_num;i++) {
            // normalize observation
            auto &obs = all_obs[i];
            for (int j = 0; j < user_num; j++) {
                auto &ob = obs[j];
                vector<uint64_t> tmp(cluster_number, 0);
                ob = gt(ob, tmp, pt, role);
//                cout<<"observation"<<endl;
//                print_share(ob,pt,role);
//                cout<<"----------------------------"<<endl;
            }
        }
        end = clock();
//        cout<<"Observation update time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        // latent truth discovery

        start = clock();
        cout<<"start latent truth"<<endl;
        auto tls = latent_truth_discovery(all_obs, LTM_ITER, pt, role);
        end = clock();
        cout<<"Latent truth time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;

        // calculate score of each user for each question
        start = clock();
        vector<vector<uint64_t>> score(question_num, vector<uint64_t>(user_num, 0));
        for (int i = 0; i < question_num; i++) {
            auto &cls = all_cls[i];
            auto &tl = tls[i];
            auto &ans = answers[i];
            for (int j = 0; j < user_num; j++) {
                auto &cl = cls[j];
                uint64_t s = 0;
                for (int k = 0; k < cl.size(); k++) {
                    s += inner_product(tl, cl[k], pt, role);
//                    print_share(tl,pt,role);
//                    print_share(cl[k], pt,role);
//                    print_share(s,pt,role);
                }
                score[i][j] = s;
//                print_share(s,pt,role);
//                cout<<"------------------------------"<<endl;
            }
            // answer selection
            auto best_user = argmax_vector(score[i], pt, role);
            vector<uint64_t> as(ANSWER_LEN, 0);
            for (int j = 0; j < user_num; j++) {
                vector<uint64_t> tmp(ANSWER_LEN, best_user[j]);
                tmp = product(ans[j], tmp, pt, role);
                for (int k = 0; k < as.size(); k++) {
                    as[k] += tmp[k];
                }
            }
        }

        auto end_all = clock();
        cout<<"TextTruth time: "<<(double)(end_all - start_all) / CLOCKS_PER_SEC<<"S"<<endl;

        // for test only return topk index for benchmark
        vector<vector<int>>topk_index(question_num,vector<int>(topK,0));
        vector<int> best_index(question_num);
        for (int i = 0; i < question_num; i++) {
            auto &ss = score[i];
            ss = open_share(ss,pt,role);
            vector<pair<int,int>>score_and_index(user_num);
            for(int j=0;j<user_num;j++) {
                score_and_index.emplace_back(ss[j],j);
            }

            sort(score_and_index.begin(), score_and_index.end(),greater<>());
            for(int j=0; j<topK; j++) {
                topk_index[i][j] = score_and_index[j].second;
            }
        }
        return topk_index;
    }
}