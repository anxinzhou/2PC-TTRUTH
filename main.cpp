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

void test_inner_product(ABYParty *pt, e_role role) {
    int len = 100;
    srand(1);
    vector<uint64_t> a(len);
    vector<uint64_t> b(len);
    uint64_t v_sum = 0;
    for(int i=0;i<len;i++) {
        a[i] = rand()%100;
        b[i] = rand()%100;
        v_sum += a[i] * b[i];
    }

    MPC::VecShare a_share (a,role);
    MPC::VecShare b_share (b,role);

    std::vector<Sharing*>& sharings = pt->GetSharings();
    ArithmeticCircuit* circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();


    auto s_out = MPC::build_inner_product_circuit(a_share.val(),b_share.val(),pt, role);
    s_out = circ->PutOUTGate(s_out, ALL);

    pt->ExecCircuit();

    auto output = s_out->get_clear_value<uint64_t>();
    pt->Reset();
//    assert(output == v_sum);
//
    std::cout << "\nCircuit Result: " << output;
    std::cout << "\nVerification Result: " << v_sum << std::endl;
}

void test_argmax(ABYParty *pt, e_role role) {
    int len = 20;
    srand(1);
    vector<uint64_t> a(len);

    for(int i=0;i<len;i++) {
        a[i] = i+1;
    }
    shuffle(a.begin(), a.end(),std::default_random_engine(1));
    int v_max = 0;
    int v_index = 0;
    for(int i=0;i<len;i++) {
        if(v_max<a[i]) {
            v_max = a[i];
            v_index = i;
        }
    }


    MPC::VecShare a_share (a,role);
    std::vector<Sharing*>& sharings = pt->GetSharings();
    BooleanCircuit* circ = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
//    ArithmeticCircuit *arc = (ArithmeticCircuit *) sharings[S_ARITH]->GetCircuitBuildRoutine();
//
//    len = 50;
//    vector<uint64_t> t1(len);
//    vector<uint64_t> t2(len);
//
//    auto x = circ->PutSIMDINGate(len,t1.data(),UINT64_LEN, SERVER);
//    auto y =  circ->PutSIMDINGate(len,t2.data(),UINT64_LEN,CLIENT);
//    auto output = circ->PutGTGate(x,y);
//    output = circ->PutOUTGate(output, ALL);

//    pt->ExecCircuit();

    auto s_out = MPC::build_argmax_circuit(a_share.val(), pt, role);
    auto s_out2 = MPC::build_argmax_circuit(a_share.val(), pt, role);

    s_out = circ->PutOUTGate(s_out, ALL);
    s_out2 = circ->PutOUTGate(s_out2, ALL);

    pt->ExecCircuit();

    auto output = s_out->get_clear_value<uint64_t>();

    pt->Reset();

    std::cout << "\nCircuit Result: " << output;
    std::cout << "\nVerification Result: " << v_index << std::endl;
    cout<<a[v_index]<<endl;
}

void test_argmax2(ABYParty *pt, e_role role) {
    int len = 10;
    srand(1);
    vector<uint64_t> a(len);

    for(int i=0;i<len;i++) {
        a[i] = i+1;
    }
    shuffle(a.begin(), a.end(),std::default_random_engine(1));
    int v_max = 0;
    int v_index = 0;
    for(int i=0;i<len;i++) {
        if(v_max<a[i]) {
            v_max = a[i];
            v_index = i;
        }
    }


    MPC::VecShare a_share (a,role);
    std::vector<Sharing*>& sharings = pt->GetSharings();
    BooleanCircuit* circ = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

    auto s_out = MPC::build_argmax_circuit2(a_share.val(), pt, role);

    s_out = circ->PutOUTGate(s_out, ALL);

    pt->ExecCircuit();

    auto output = s_out->get_clear_value<uint64_t>();

    pt->Reset();

    std::cout << "\nCircuit Result: " << output;
    std::cout << "\nVerification Result: " << v_index << std::endl;
    cout<<a[v_index]<<endl;
}

void test_argmax3(ABYParty *pt, e_role role) {
    int len = 10;
    srand(1);
    vector<uint64_t> a(len);

    for(int i=0;i<len;i++) {
        a[i] = i+1;
    }
    shuffle(a.begin(), a.end(),std::default_random_engine(1));
    int v_max = 0;
    int v_index = 0;
    for(int i=0;i<len;i++) {
        if(v_max<a[i]) {
            v_max = a[i];
            v_index = i;
        }
    }


    MPC::VecShare a_share (a,role);
    std::vector<Sharing*>& sharings = pt->GetSharings();
    auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();

    auto s_out = MPC::argmax(a_share.val(), pt, role);


    std::cout << "\nCircuit Result: " << s_out;
    std::cout << "\nVerification Result: " << v_index << std::endl;
    cout<<a[v_index]<<endl;
}

void test_argmax4(ABYParty *pt, e_role role) {
    int len = 10;
    srand(1);
    vector<uint64_t> a(len);

    for(int i=0;i<len;i++) {
        a[i] = i+1;
    }
    shuffle(a.begin(), a.end(),std::default_random_engine(1));
    int v_max = 0;
    int v_index = 0;
    for(int i=0;i<len;i++) {
        if(v_max<a[i]) {
            v_max = a[i];
            v_index = i;
        }
    }


    MPC::VecShare a_share (a,role);
    std::vector<Sharing*>& sharings = pt->GetSharings();
    auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();

    auto s_out = MPC::argmax_vector(a_share.val(), pt, role);

    for(uint i=0;i<len;i++) {
        if(s_out[i] == 1) {
            std::cout << "\nCircuit Result: " << i<<endl;

        }
    }
    std::cout << "\nVerification Result: " << v_index << std::endl;
    cout<<a[v_index]<<endl;
}

void test_max(ABYParty *pt, e_role role) {
    int len = 2;
    srand(1);
    vector<uint64_t> a(len);

    for(int i=0;i<len;i++) {
        a[i] = i+1;
    }
    shuffle(a.begin(), a.end(),std::default_random_engine(1));
    int v_max = 0;
    int v_index = 0;
    for(int i=0;i<len;i++) {
        if(v_max<a[i]) {
            v_max = a[i];
            v_index = i;
        }
    }


    MPC::VecShare a_share (a,role);
    std::vector<Sharing*>& sharings = pt->GetSharings();
    BooleanCircuit* circ = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();

    auto s_out = MPC::build_max_circuit(a_share.val(), pt, role);

    s_out = circ->PutOUTGate(s_out, ALL);

    pt->ExecCircuit();

    auto output = s_out->get_clear_value<uint64_t>();

    pt->Reset();

    std::cout << "\nCircuit Result: " << output;
    std::cout << "\nVerification Result: " << v_max << std::endl;
    cout<<a[v_index]<<endl;
}

void test_min(ABYParty *pt, e_role role) {
    uint64_t a =5;
    uint64_t b= 10;
    auto output = MPC::min(a,b,pt,role);
    cout<<output<<endl;
}

void test_max2N(ABYParty *pt, e_role role) {
    uint64_t a=1024;
    cout<< MPC::max2N(a,pt,role) <<endl;
}

void test_right_shift(ABYParty *pt, e_role role) {
    srand(1);
    uint64_t a = 2048;
    uint64_t tmp = rand();
    if (role == SERVER) {
        a = tmp;
    } else {
        a = a-tmp;
    }
    auto res = MPC::right_shift(a,3,pt,role);
//    cout<<a<<endl;
    cout<<res<<endl;
}

void test_log(ABYParty *pt, e_role role) {
    clock_t start, end;
    start = clock();
    for(int i=0; i<1;i++) {
        uint64_t a = 21<<FLOAT_SCALE_FACTOR;
        srand(2);
        uint64_t x = rand();
        a = role==SERVER? a-x:x;

        uint64_t res = MPC::log(a,FLOAT_SCALE_FACTOR,FLOAT_SCALE_FACTOR,pt,role);
            cout<<res<<endl;
    }
    end = clock();
    cout<<"Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//    cout<<res<<endl;
}

void test_rep_square_root(ABYParty *pt, e_role role) {
    clock_t start, end;
    start = clock();
    for(int i=0; i<1;i++) {
        uint64_t a = (0)<<20;
        srand(2);
        uint64_t x = rand();
        a = role==SERVER? a-x:x;
        auto res = MPC::rep_square_root(a, 20,20, pt, role);
        cout<<res<<endl;
    }
    end = clock();
    cout<<"Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//    cout<<res<<endl;
}

void test_sigmoid(ABYParty *pt, e_role role) {
    clock_t start, end;
    start = clock();
    for(int i=0; i<1;i++) {
        uint64_t a = (3)<<16;
        srand(2);
        uint64_t x = rand();
        a = role==SERVER? a-x:x;
        auto res = MPC::sigmoid(a, 16,16, pt, role);
        cout<<res<<endl;
    }
    end = clock();
    cout<<"Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
//    cout<<res<<endl;
}

void test_eq(ABYParty *pt, e_role role) {
    auto res = MPC::eq(3,2,pt,role);
    cout<<res<<endl;
}

void test_right_shift2(ABYParty *pt, e_role role) {
    uint64_t a = role==SERVER? 4:5;
    auto res = MPC::right_shift(a,1,pt,role);
    cout<<res<<endl;
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

