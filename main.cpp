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

int32_t read_test_options(int32_t* argcp, char*** argvp, e_role* role,
                          uint32_t* bitlen, uint32_t* numbers, uint32_t* secparam, std::string* address,
                          uint16_t* port, int32_t* test_op) {

    uint32_t int_role = 0, int_port = 0;

    parsing_ctx options[] =
            { { (void*) &int_role, T_NUM, "r", "Role: 0/1", true, false },
              { (void*) numbers, T_NUM, "n",	"Number of elements for inner product, default: 128", false, false },
              {	(void*) bitlen, T_NUM, "b", "Bit-length, default 16", false, false },
              { (void*) secparam, T_NUM, "s", "Symmetric Security Bits, default: 128", false, false },
              {	(void*) address, T_STR, "a", "IP-address, default: localhost", false, false },
              {	(void*) &int_port, T_NUM, "p", "Port, default: 7766", false, false },
              { (void*) test_op, T_NUM, "t", "Single test (leave out for all operations), default: off",
                      false, false } };

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



vector<vector<vector<string>>> load_keyword() {
    const string kv_bath_dir = "../answer_grading/keywords";
    const vector<string> qdir{"q1","q2","q3","q4","q5","q6","q7","q8","q9","q10","q11","q12"};
}



vector<vector<double>> load_score() {
    const string score_bath_dir = "../answer_grading/scores";
    const vector<string> qdir{"q1","q2","q3","q4","q5","q6","q7","q8","q9","q10","q11","q12"};
    vector<vector<double>> score;
    for(int i=0; i<qdir.size(); i++) {
        string dirpath = score_bath_dir + "/" + qdir[i];
        filesystem::directory_iterator end_itr;
        for(filesystem::directory_iterator itr(dirpath);itr!=end_itr;++itr) {
            if(filesystem::is_regular_file(itr->path())) {
                auto current_file = itr->path();
                int question_id = atof(current_file.filename().c_str());
                ifstream score_file(current_file.string());
                vector<double> sub_scores;
                if (!score_file.is_open()) {
                    cout<< "can not open score file"<<endl;
                    cout << strerror(errno) << endl;
                    exit(0);
                }
                float score;
                while (score_file >> score) {
                    sub_scores.push_back(score);
                }
                scores[question_id] = sub_scores;
            }
        }
    }

    int question_num = questions.size();
    scores = vector<vector<float>> (question_num,vector<float>());
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr(score_file_path); itr != end_itr; ++itr) {
// If it's not a directory, list it. If you want to list directories too, just remove this check.
        if (boost::filesystem::is_regular_file(itr->path())) {
// assign current file name to current_file and echo it out to the console.
            auto current_file = itr->path();
            int question_id = atof(current_file.filename().c_str());
            ifstream score_file(current_file.string());
            vector<float> sub_scores;
            if (!score_file.is_open()) {
                cout<< "can not open score file"<<endl;
                cout << strerror(errno) << endl;
                exit(0);
            }
            float score;
            while (score_file >> score) {
                sub_scores.push_back(score);
            }
            scores[question_id] = sub_scores;
        }
    }


}

unordered_map<string,vector<double>> load_pretrained_vec() {
    const string model_path = "../answer_grading/pretrained_vec";

}


void AnswerGradingData::load_dataset(const string &dir_path) {
    const string question_file_path = (boost::filesystem::path(dir_path) / "questions/questions").string();
    const string answers_file_path = (boost::filesystem::path(dir_path) / "answers").string();
    const string score_file_path = (boost::filesystem::path(dir_path) / "scores").string();
    load_question(question_file_path);
    load_answers(answers_file_path);
    load_scores(score_file_path);
}

void AnswerGradingData::load_question(const string &question_file_path) {
    ifstream question_file(question_file_path);
    if (!question_file.is_open()) {
        cout<<"can not open question file"<<endl;
        cout << strerror(errno) << endl;
        exit(0);
    }

    string line;
    while (getline(question_file, line)) {
        questions.emplace_back(std::move(line),key_factor_number);
    }
    question_file.close();
}

void AnswerGradingData::load_answers(const string &answer_file_path) {
    // read answer
    int question_num = questions.size();
    answers = vector<vector<Answer>> (question_num,vector<Answer>());
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr(answer_file_path); itr != end_itr; ++itr) {
        // If it's not a directory, list it. If you want to list directories too, just remove this check.
        if (boost::filesystem::is_regular_file(itr->path())) {
            // assign current file name to current_file and echo it out to the console.
            auto current_file = itr->path();
            int question_id = atof(current_file.filename().c_str());
            ifstream answers_file(current_file.string());
            vector<Answer> sub_answers;
            if (!answers_file.is_open()) {
                cout<< "can not open answer file"<<endl;
                cout << strerror(errno) << endl;
                exit(0);
            }
            string answer;
            int count = 0;
            while (getline(answers_file, answer)) {
                sub_answers.emplace_back(count, question_id, std::move(answer));
                count+=1;
            }

            answers[question_id] = sub_answers;
        }
    }
}

void AnswerGradingData::load_scores(const string &score_file_path) {
    int question_num = questions.size();
    scores = vector<vector<float>> (question_num,vector<float>());
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr(score_file_path); itr != end_itr; ++itr) {
        // If it's not a directory, list it. If you want to list directories too, just remove this check.
        if (boost::filesystem::is_regular_file(itr->path())) {
            // assign current file name to current_file and echo it out to the console.
            auto current_file = itr->path();
            int question_id = atof(current_file.filename().c_str());
            ifstream score_file(current_file.string());
            vector<float> sub_scores;
            if (!score_file.is_open()) {
                cout<< "can not open score file"<<endl;
                cout << strerror(errno) << endl;
                exit(0);
            }
            float score;
            while (score_file >> score) {
                sub_scores.push_back(score);
            }
            scores[question_id] = sub_scores;
        }
    }
}

int main(int argc, char** argv) {

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

//    test_argmax4(pt, role);
//    test_inner_product(pt, role);
//    test_right_shift(pt, role);
//test_log(pt,role);
test_rep_square_root(pt,role);
//test_eq(pt,role);
//    test_right_shift2(pt,role);
//test_sigmoid(pt,role);

    delete pt;
    return 0;
}

