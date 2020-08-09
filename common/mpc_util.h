/**
 \file 		innerproduct.h
 \author 	sreeram.sadasivam@cased.de
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
 \brief		Implementation of the Inner Product using ABY Framework.
 */

#ifndef __INNERPRODUCT_H_
#define __INNERPRODUCT_H_

#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/aby/abyparty.h"
#include <math.h>
#include <vector>
#include <cassert>
#include <iostream>

extern const uint UINT64_LEN;
extern const uint UINT_LEN;
extern const uint FLOAT_SCALE_FACTOR;

using namespace std;

namespace MPC {
    class Share {
    public:
        uint64_t a;
        uint64_t b;
        e_role role;

        Share() {}

        Share(uint64_t x, e_role r) {
            role = r;
            srand(1);
            a = rand();
            b = x - a;
        }

        uint64_t &val() {
            return role == SERVER ? a : b;
        }
    };

    class VecShare {
    public:
        vector<uint64_t> a;
        vector<uint64_t> b;
        e_role role;

        VecShare() {}

        VecShare(vector<uint64_t> &x, e_role r) {
            role = r;
            srand(1);
            for (auto v:x) {
                uint64_t v1 = rand();
                uint64_t v2 = v - v1;
                a.push_back(v1);
                b.push_back(v2);
            }
        }

        vector<uint64_t> &val() {
            return role == SERVER ? a : b;
        }
    };

    class MulTripletVec {
    public:
        VecShare a;
        VecShare b;
        VecShare c;

//    e_role role;
        explicit MulTripletVec(int len, e_role r) {
//        role = r;
            vector<uint64_t> x;
            vector<uint64_t> y;
            vector<uint64_t> z;
            srand(2);
            for (uint64_t i = 0; i < len; i++) {
                uint64_t p = rand();
                uint64_t q = rand();
                x.push_back(p);
                y.push_back(q);
                z.push_back(p * q);
            }
            a = VecShare(x, r);
            b = VecShare(y, r);
            c = VecShare(z, r);
        }

        vector<uint64_t> &getA() {
            return a.val();
        }

        vector<uint64_t> &getB() {
            return b.val();
        }

        vector<uint64_t> &getC() {
            return c.val();
        }
    };

    class MulTriplet {
    public:
        Share a;
        Share b;
        Share c;

//    e_role role;
        explicit MulTriplet(e_role r) {
//        role = r;
            srand(2);
            uint64_t p = rand();
            uint64_t q = rand();
            a = Share(p, r);
            b = Share(q, r);
            c = Share(p*q, r);
        }

        uint64_t &getA() {
            return a.val();
        }

        uint64_t &getB() {
            return b.val();
        }

        uint64_t &getC() {
            return c.val();
        }
    };

    ABYParty *init_party(e_role role, const std::string &address, uint16_t port, seclvl seclvl, uint32_t bitlen,
                         uint32_t nthreads, e_mt_gen_alg mt_alg);

    uint64_t share_gt_const(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    uint64_t gt(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    uint64_t eq(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    vector<uint64_t> eq(vector<uint64_t>&a, vector<uint64_t>&b, ABYParty *pt, e_role role);
    uint64_t share_eq_const(uint64_t a, uint64_t b, ABYParty *pt, e_role role);

    vector<uint64_t> argmax_vector(vector<uint64_t>&a, ABYParty *pt, e_role role);
    vector<uint64_t> argmin_vector(vector<uint64_t>&a, ABYParty *pt, e_role role);
    uint64_t min(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    uint64_t argmax(vector<uint64_t>&a, ABYParty *pt, e_role role);
    uint64_t argmax_test(vector<uint64_t>&a, ABYParty *pt, e_role role);

    uint64_t product(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    vector<uint64_t> product(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty*pt, e_role role);
    uint64_t inner_product(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty *pt, e_role role);
    vector<uint64_t> minus(vector<uint64_t> &a, vector<uint64_t> &b);

    uint64_t rep_square_root(uint64_t a, uint64_t scale_factor,uint64_t already_scaled_factor, ABYParty *pt, e_role role);
    uint64_t log(uint64_t a, uint64_t scale_factor, uint64_t already_scaled_factor,ABYParty *pt,e_role role);
    uint64_t sigmoid(uint64_t a, uint64_t scale_factor,uint64_t already_scaled_factor,ABYParty *pt, e_role role);

    uint64_t max2N(uint64_t a, uint64_t &digits, ABYParty *pt, e_role role);
    uint64_t max2N(uint64_t a, ABYParty *pt, e_role role) ;
    uint64_t right_shift(uint64_t a, uint64_t digits, ABYParty *pt, e_role role);
    uint64_t right_shift_const(uint64_t a, uint64_t digits, ABYParty *pt, e_role role);
    vector<uint64_t> right_shift_const(vector<uint64_t >&a, uint64_t digits, ABYParty *pt, e_role role);
    vector<uint64_t> right_shift(vector<uint64_t >&a, uint64_t digits, ABYParty *pt, e_role role);
    uint64_t left_shift(uint64_t a, uint64_t digits, ABYParty *pt, e_role role);
    uint64_t left_shift_const(uint64_t a, uint64_t digits, ABYParty *pt, e_role role);

    uint range_sample(vector<uint64_t>&threshold, vector<uint64_t> &value, uint64_t p, ABYParty *pt, e_role role);
    // generate a 32-bit randomness
    uint random(ABYParty *pt, e_role role);

    share*build_gt_circuit(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    share*build_eq_circuit(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    share* build_product_circuit(uint64_t a, uint64_t b, ABYParty *pt, e_role role);
    share* build_product_circuit(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty*pt, e_role role);
    share* build_inner_product_circuit(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty *pt, e_role role);
    share *build_max_circuit(vector<uint64_t>&a, ABYParty *pt, e_role role);
    share *build_argmax_circuit(vector<uint64_t>&a, ABYParty *pt, e_role role);
    share *build_argmax_circuit2(vector<uint64_t>&a, ABYParty *pt, e_role role);
    share *build_argmax_vector_circuit(vector<uint64_t>&a, ABYParty *pt, e_role role);
    share* build_max2N_circuit(uint64_t a, share ** digits, ABYParty *pt, e_role role);
}


#endif
