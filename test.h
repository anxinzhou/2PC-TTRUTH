//
// Created by anxin on 8/9/20.
//

#ifndef TTRUTH_TEST_H
#define TTRUTH_TEST_H

#include <ABY/src/abycore/aby/abyparty.h>
#include <vector>
#include "abycore/aby/abyparty.h"

#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"
#include "common/mpc_util.h"
#include <algorithm>
#include <random>
using namespace std;
void test_inner_product(ABYParty *pt, e_role role);
void test_argmax(ABYParty *pt, e_role role);
void test_argmax2(ABYParty *pt, e_role role);
void test_argmax3(ABYParty *pt, e_role role);
void test_max(ABYParty *pt, e_role role);
void test_min(ABYParty *pt, e_role role);
void test_max2N(ABYParty *pt, e_role role);
void test_right_shift(ABYParty *pt, e_role role);
void test_log(ABYParty *pt, e_role role);
void test_rep_square_root(ABYParty *pt, e_role role);
void test_sigmoid(ABYParty *pt, e_role role);
void test_eq(ABYParty *pt, e_role role);
void test_right_shift2(ABYParty *pt, e_role role);
#endif //TTRUTH_TEST_H
