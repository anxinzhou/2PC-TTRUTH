/**
 \file 		innerproduct_test.cpp
 \author	sreeram.sadasivam@cased.de
 \copyright	ABY - A Framework for Efficient Mixed-protocol Secure Two-party Computation
			Copyright (C) 2019 Engineering Cryptographic Protocols Group, TU Darmstadt
			This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU Lesser General Public License as published
            by the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.
            ABY is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
            GNU Lesser General Public License for more details.
            You should have received a copy of the GNU Lesser General Public License
            along with this program. If not, see <http://www.gnu.org/licenses/>.
 \brief		Inner Product Test class implementation.
 */

//Utility libs
#include <ENCRYPTO_utils/crypto/crypto.h>
#include <ENCRYPTO_utils/parse_options.h>
//ABY Party class
#include "abycore/aby/abyparty.h"

#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"
#include "common/mpc_util.h"

#include <algorithm>
#include <random>
#include <filesystem>
#include <iterator>

#include "ttruth/ttruth.h"
#include "mpc_ttruth/mpc_ttruth.h"

const int QUESTION_NUM = 87;

int32_t read_test_options(int32_t *argcp, char ***argvp, e_role *role,
                          uint32_t *bitlen, uint32_t *numbers, uint32_t *secparam, std::string *address,
                          uint16_t *port, int32_t *test_op) {

    uint32_t int_role = 0, int_port = 0;

    parsing_ctx options[] =
            {{(void *) &int_role, T_NUM, "r", "Role: 0/1",                                          true,  false},
             {(void *) numbers,   T_NUM, "n", "Number of elements for inner product, default: 128", false, false},
             {(void *) bitlen,    T_NUM, "b", "Bit-length, default 16",                             false, false},
             {(void *) secparam,  T_NUM, "s", "Symmetric Security Bits, default: 128",              false, false},
             {(void *) address,   T_STR, "a", "IP-address, default: localhost",                     false, false},
             {(void *) &int_port, T_NUM, "p", "Port, default: 7766",                                false, false},
             {(void *) test_op,   T_NUM, "t", "Single test (leave out for all operations), default: off",
                                                                                                    false, false}};

    if (!parse_options(argcp, argvp, options,
                       sizeof(options) / sizeof(parsing_ctx))) {
        print_usage(*argvp[0], options, sizeof(options) / sizeof(parsing_ctx));
        std::cout << "Exiting" << std::endl;
        exit(0);
    }

    assert(int_role < 2);
    *role = (e_role) int_role;

    if (int_port != 0) {
        assert(int_port < 1 << (sizeof(uint16_t) * 8));
        *port = (uint16_t) int_port;
    }

    return 1;
}

vector<string> string_split(const string &s) {
    vector<string> tokens;
    istringstream iss(s);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
    return tokens;
}

vector<vector<vector<string>>> load_keyword() {
    const string kv_bath_dir = "../answer_grading/keywords";
    vector<vector<vector<string>>> all_keywords(QUESTION_NUM);
    // read answer
    filesystem::directory_iterator end_itr;

    for (filesystem::directory_iterator itr(kv_bath_dir); itr != end_itr; ++itr) {
        // If it's not a directory, list it. If you want to list directories too, just remove this check.
        if (filesystem::is_regular_file(itr->path())) {
            // assign current file name to current_file and echo it out to the console.
            auto current_file = itr->path();
            int question_id = atof(current_file.filename().c_str());
            ifstream keywords_file(current_file.string());
            vector<vector<string>> user_keywords;
            if (!keywords_file.is_open()) {
                cout<< "can not open answer file"<<endl;
                cout << strerror(errno) << endl;
                exit(0);
            }
            string raw_keywords;
            while (getline(keywords_file, raw_keywords)) {
                auto keywords = string_split(raw_keywords);
                user_keywords.push_back(std::move(keywords));
            }
            all_keywords[question_id] = (std::move(user_keywords));
        }
    }

    return all_keywords;
}


vector<vector<double>> load_score() {
    const string score_bath_dir = "../answer_grading/scores";
    vector<vector<double>> score(QUESTION_NUM);

    filesystem::directory_iterator end_itr;
    for (filesystem::directory_iterator itr(score_bath_dir); itr != end_itr; ++itr) {
        if (filesystem::is_regular_file(itr->path())) {
            auto current_file = itr->path();
            int question_id = atof(current_file.filename().c_str());
            ifstream score_file(current_file.string());
            vector<double> sub_scores;
            if (!score_file.is_open()) {
                cout << "can not open score file" << endl;
                cout << strerror(errno) << endl;
                exit(0);
            }
            double s;
            while (score_file >> s) {
                sub_scores.push_back(s);
            }
            score[question_id] = std::move(sub_scores);
        }
    }
    return score;
}

unordered_map<string, vector<double>> load_pretrained_vec() {
    const string model_path = "../answer_grading/pretrained_vec";
    unordered_map<string, vector<double>> pretrained_vec;

    ifstream word_embedding_file(model_path);
    if (!word_embedding_file.is_open()) {
        cout << "can not open word embedding file" << endl;
        cout << strerror(errno) << endl;
        exit(0);
    }
    string line;
    while (getline(word_embedding_file, line)) {
        istringstream iss(line);
        string word;
        iss >> word;
        vector<double> vec;
        double v;
        while (iss >> v) {
            vec.push_back(v);
        }
        pretrained_vec.emplace(word,vec);

    }
    word_embedding_file.close();
//    cout<<pretrained_vec.size()<<endl;
    return pretrained_vec;
}

void test_ttruth() {
    auto pretrained_vec = load_pretrained_vec();
    auto all_keywords = load_keyword();
    auto score = load_score();

    if(score.size()!=all_keywords.size()) {
        cout<<score.size()<<" "<<all_keywords.size()<<endl;
        cout<<"wrong size at line"<<__LINE__<<" "<<__FILE__<<endl;
        exit(-1);
    }
    int question_num = all_keywords.size();
    int user_num = all_keywords[0].size();  // max 24

    vector<vector<vector<vector<double>>>> all_kvec(question_num,
                                                    vector<vector<vector<double>>>(user_num));

    for(int i=0; i<question_num; i++) {
        for(int j=0; j<user_num; j++) {
            auto&keywords = all_keywords[i][j];
            vector<vector<double>> key_vecs;
            for (auto &w:keywords) {
                key_vecs.push_back(pretrained_vec.at(w));
            }
            all_kvec[i][j] = std::move(key_vecs);
        }
    }

    int topK=24;
    auto topk_index = ttruth(all_kvec, topK);
    vector<double>total_avg(topK,0);
    double all_socre = 0;
    int all_count = question_num * topK;
    for(int i=0;i<question_num;i++) {
        vector<double>avg(topK,0);
        cout<<"question "<<i+1<<" ";
        for(int j=0;j<topK;j++) {
            int index = topk_index[i][j];
            cout<<score[i][index] << " ";
            avg[j] += score[i][index];
            if(j!=0) {
                avg[j] += avg[j-1];
            }
            all_socre+=score[i][index];
        }
        cout<<endl;
        for(int j=0;j<topK;j++) {
            avg[j] /= (j+1);
            total_avg[j] += avg[j];
        }
    }
    for(int i=0; i<topK;i++) {
        total_avg[i]/=question_num;
        cout<< total_avg[i]<<" ";
    }
    cout<<endl;
    cout<<"avg: " << all_socre/all_count;
    cout<<endl;
    vector<double>total_avg2(topK,0);
    for(int i=0;i<question_num;i++) {
        vector<double>avg(topK,0);
        for(int j=0;j<topK;j++) {
            int index = topk_index[i][j];
//            cout<<score[i][index] << " ";
            total_avg2[j] += score[i][index];
        }
    }
//    cout<<endl;
    for(int i=0; i<topK;i++) {
        total_avg2[i]/=question_num;
        cout<< total_avg2[i]<<" ";
    }

}

void testMPCTextTruth(ABYParty *pt, e_role role) {
    auto pretrained_vec = load_pretrained_vec();
    auto all_keywords = load_keyword();
    auto score = load_score();

    if(score.size()!=all_keywords.size()) {
        cout<<score.size()<<" "<<all_keywords.size()<<endl;
        cout<<"wrong size at line"<<__LINE__<<" "<<__FILE__<<endl;
        exit(-1);
    }
    int question_num = all_keywords.size();
    int user_num = all_keywords[0].size();  // max 24

    vector<vector<vector<vector<uint64_t>>>> all_kvec(question_num,
                                                    vector<vector<vector<uint64_t>>>(user_num));

    // convert keywords vector to sharing vector
    for(int i=0; i<question_num; i++) {
        for(int j=0; j<user_num; j++) {
            auto&keywords = all_keywords[i][j];
            vector<vector<uint64_t>> key_vecs;
            for (auto &w:keywords) {
                vector<double>tmp(pretrained_vec.at(w));
                vector<uint64_t>vec(tmp.size(),0);
                for(int k=0; k<tmp.size(); k++) {
                    double t = tmp[k] * FLOAT_SCALE_FACTOR;
                    if (t<0) {
                        vec[k] = uint64_t(-t);
                    } else {
                        vec[k] = t;
                    }
                    auto s = MPC::Share(vec[k], role);
                    vec[k] = s.val();
                }
                key_vecs.push_back(std::move(vec));
            }
            all_kvec[i][j] = std::move(key_vecs);
        }
    }

    vector<vector<vector<uint64_t>>> answers(question_num,
                                             vector<vector<uint64_t>>(user_num,
                                                     vector<uint64_t>(150,0)));

    int topK=24;
    auto topk_index = MPC::ttruth(all_kvec, answers, topK, pt, role);
    vector<double>total_avg(topK,0);
    double all_socre = 0;
    int all_count = question_num * topK;
    for(int i=0;i<question_num;i++) {
        vector<double>avg(topK,0);
        cout<<"question "<<i+1<<" ";
        for(int j=0;j<topK;j++) {
            int index = topk_index[i][j];
            cout<<score[i][index] << " ";
            avg[j] += score[i][index];
            if(j!=0) {
                avg[j] += avg[j-1];
            }
            all_socre+=score[i][index];
        }
        cout<<endl;
        for(int j=0;j<topK;j++) {
            avg[j] /= (j+1);
            total_avg[j] += avg[j];
        }
    }
    for(int i=0; i<topK;i++) {
        total_avg[i]/=question_num;
        cout<< total_avg[i]<<" ";
    }
    cout<<endl;
    cout<<"avg: " << all_socre/all_count;
    cout<<endl;
    vector<double>total_avg2(topK,0);
    for(int i=0;i<question_num;i++) {
        vector<double>avg(topK,0);
        for(int j=0;j<topK;j++) {
            int index = topk_index[i][j];
//            cout<<score[i][index] << " ";
            total_avg2[j] += score[i][index];
        }
    }
//    cout<<endl;
    for(int i=0; i<topK;i++) {
        total_avg2[i]/=question_num;
        cout<< total_avg2[i]<<" ";
    }
}

int main(int argc, char **argv) {

    e_role role;
    uint32_t bitlen = UINT64_LEN, numbers = 128, secparam = 128, nthreads = 1;
    uint16_t port = 7766;
    std::string address = "127.0.0.1";
    int32_t test_op = -1;
    e_mt_gen_alg mt_alg = MT_OT;

    read_test_options(&argc, &argv, &role, &bitlen, &numbers, &secparam, &address, &port, &test_op);

    seclvl seclvl = get_sec_lvl(secparam);

    // call inner product routine. set size with cmd-parameter -n <size>
    ABYParty *pt = MPC::init_party(role, address, port, seclvl, UINT64_LEN, nthreads, mt_alg);

//    test_ttruth();
testMPCTextTruth(pt,role);
    delete pt;
    return 0;
}

