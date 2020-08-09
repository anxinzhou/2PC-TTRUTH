/**
 \file 		innerproduct.cpp
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

#include "mpc_util.h"
#include "abycore/sharing/sharing.h"

const uint UINT64_LEN = 64;
const uint UINT_LEN=32;
const uint FLOAT_SCALE_FACTOR =16;


namespace MPC {
    ABYParty* init_party(e_role role, const std::string& address, uint16_t port, seclvl seclvl,uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg) {
        return new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                            mt_alg);
    }

    uint random(ABYParty *pt, e_role role) {
//        uint seed = role==SERVER? 2:3;
        srand(time(NULL));
        uint l = rand();
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto r = bcirc->PutSharedINGate(l,UINT_LEN);
        r = acirc->PutB2AGate(r);
        r = acirc->PutSharedOUTGate(r);
        pt->ExecCircuit();
        uint res = r->get_clear_value<uint>();
        return res;
    }

    uint64_t argmax(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto s_out = build_argmax_circuit2(a,  pt, role);
        s_out = circ->PutSharedOUTGate(s_out);
        pt->ExecCircuit();
        auto output = s_out->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    uint64_t max2N(uint64_t a, ABYParty *pt, e_role role) {
        uint64_t v = 0;
        return max2N(a,v,pt,role);
    }

    uint64_t max2N2(uint64_t a, uint64_t &digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        sa = bcirc->PutA2YGate(sa);

//        cout<<sa->get_bitlength()<<endl;
//        sa = ycirc->PutSplitterGate(sa);
        uint64_t sum = 0;
        digits = 0;
        vector<share*>inds;
        for(int i=1;i<UINT64_LEN;i++) {
            sa->set_wire_id(i, ycirc->PutORGate(sa->get_wire_id(i),sa->get_wire_id(i-1)));
        }
        for(int i=UINT64_LEN-1;i>=0;i--) {
            auto sid = acirc->PutY2AGate(sa->get_wire_ids_as_share(i),bcirc);
            sid = acirc->PutSharedOUTGate(sid);
            inds.push_back(sid);
        }
        pt->ExecCircuit();
        for(int i=UINT64_LEN-1; i>=0; i--) {
            uint64_t v = uint64_t(1)<<uint64_t(i);
            uint64_t ind = inds[i]->get_clear_value<uint64_t>();
            sum+= ind * v;
            digits+=ind;
        }
        pt->Reset();
        if(role == SERVER) {
            sum+=1;
        }

        return sum;
    }

    uint64_t max2N(uint64_t a, uint64_t &digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        sa = ycirc->PutA2YGate(sa);
        uint64_t sum = 0;
        digits = 0;
        vector<share*>inds(UINT64_LEN);
        for(int i=UINT64_LEN-1;i>=0;i--) {
            uint64_t v = uint64_t(1)<<uint64_t(i);
            auto sv = ycirc->PutINGate(v,UINT64_LEN,SERVER);
            auto cmp_res = ycirc->PutGTGate(sa,sv);
            cmp_res = acirc->PutY2AGate(cmp_res, bcirc);
            cmp_res = acirc->PutSharedOUTGate(cmp_res);
            inds[i]= cmp_res;
        }
        pt->ExecCircuit();
        for(int i=UINT64_LEN-1; i>=0; i--) {
            uint64_t v = uint64_t(1)<<uint64_t(i);
            uint64_t ind = inds[i]->get_clear_value<uint64_t>();
            sum+= ind * v;
            digits+=ind;
        }
        pt->Reset();
        if(role == SERVER) {
            sum+=1;
        }

        return sum;
    }

    vector<uint64_t> right_shift_const(vector<uint64_t >&a, uint64_t digits, ABYParty *pt, e_role role) {
        uint dim = a.size();
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sdigit = acirc->PutCONSGate(digits, UINT64_LEN);
        auto sa = acirc->PutSharedSIMDINGate(dim, a.data(), UINT64_LEN);
        sdigit = ycirc->PutA2YGate(sdigit);
        sa = ycirc->PutA2YGate(sa);
        sa = ycirc->PutBarrelRightShifterGate(sa, sdigit);
        sa = acirc->PutY2AGate(sa, bcirc);
        sa = acirc->PutSharedOUTGate(sa);
        pt->ExecCircuit();
        uint *v;
        uint bitlen, nval;
        sa->get_clear_value_vec(&v,&bitlen,&nval);
        pt->Reset();
        vector<uint64_t>output(v,v+dim);
        delete v;
        return output;
    }

    vector<uint64_t> right_shift(vector<uint64_t >&a, uint64_t digits, ABYParty *pt, e_role role) {
        uint dim = a.size();
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sdigit = acirc->PutSharedINGate(digits, UINT64_LEN);
        auto sa = acirc->PutSharedSIMDINGate(dim, a.data(), UINT64_LEN);
        sdigit = ycirc->PutA2YGate(sdigit);
        sa = ycirc->PutA2YGate(sa);
        sa = ycirc->PutBarrelRightShifterGate(sa, sdigit);
        sa = acirc->PutY2AGate(sa, bcirc);
        sa = acirc->PutSharedOUTGate(sa);
        pt->ExecCircuit();
        uint *v;
        uint bitlen, nval;
        sa->get_clear_value_vec(&v,&bitlen,&nval);
        pt->Reset();
        vector<uint64_t>output(v,v+dim);
        delete v;
        return output;
    }

    uint64_t right_shift_const(uint64_t a, uint64_t digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sdigit = acirc->PutCONSGate(digits, UINT64_LEN);
        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        sdigit = ycirc->PutA2YGate(sdigit);
        sa = ycirc->PutA2YGate(sa);
        sa = ycirc->PutBarrelRightShifterGate(sa, sdigit);
        sa = acirc->PutY2AGate(sa, bcirc);
        sa = acirc->PutSharedOUTGate(sa);
        pt->ExecCircuit();
        uint64_t v = sa->get_clear_value<uint64_t>();
        pt->Reset();
        return v;
    }

    uint64_t right_shift(uint64_t a, uint64_t digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sdigit = acirc->PutSharedINGate(digits, UINT64_LEN);
        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        sdigit = ycirc->PutA2YGate(sdigit);
        sa = ycirc->PutA2YGate(sa);
        sa = ycirc->PutBarrelRightShifterGate(sa, sdigit);
        sa = acirc->PutY2AGate(sa, bcirc);
        sa = acirc->PutSharedOUTGate(sa);
        pt->ExecCircuit();
        uint64_t v = sa->get_clear_value<uint64_t>();
        pt->Reset();
        return v;
    }

    uint64_t left_shift(uint64_t a, uint64_t digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sdigit = acirc->PutSharedINGate(digits, UINT64_LEN);
        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        sdigit = ycirc->PutA2YGate(sdigit);
        sa = ycirc->PutA2YGate(sa);
        sa = ycirc->PutBarrelLeftShifterGate(sa, sdigit);
        sa = acirc->PutY2AGate(sa, bcirc);
        sa = acirc->PutSharedOUTGate(sa);
        pt->ExecCircuit();
        uint64_t v = sa->get_clear_value<uint64_t>();
        pt->Reset();
        return v;
    }

    uint64_t left_shift_const(uint64_t a, uint64_t digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sdigit = acirc->PutCONSGate(digits, UINT64_LEN);
        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        sdigit = ycirc->PutA2YGate(sdigit);
        sa = ycirc->PutA2YGate(sa);
        sa = ycirc->PutBarrelLeftShifterGate(sa, sdigit);
        sa = acirc->PutY2AGate(sa, bcirc);
        sa = acirc->PutSharedOUTGate(sa);
        pt->ExecCircuit();
        uint64_t v = sa->get_clear_value<uint64_t>();
        pt->Reset();
        return v;
    }

    // the result will scale scale_factor (1<<scale_factor)
    uint64_t log(uint64_t a, uint64_t scale_factor, uint64_t already_scaled_factor, ABYParty *pt, e_role role) {
        uint64_t digits;
        uint64_t m2N = max2N(a, digits,pt, role);

        uint64_t threshold = uint64_t(0.85 * (1<<FLOAT_SCALE_FACTOR));
        uint64_t alpha1 = 1.86511335 * (1<<FLOAT_SCALE_FACTOR) ;
        uint64_t beta1 = - uint64_t (1.8154986 * (1<<scale_factor));
        uint64_t alpha2 = 1.5617682 * (1<<FLOAT_SCALE_FACTOR);
        uint64_t beta2 = - uint64_t (1.55872625 * (1<<scale_factor));

        uint64_t tmp = m2N*threshold;
        tmp = right_shift_const(tmp, FLOAT_SCALE_FACTOR, pt, role);
        uint64_t cmp_res = gt(tmp, a, pt, role);
        uint64_t cmp_inv = -cmp_res;
        if(role==SERVER) {
            cmp_inv += 1;
        }

        uint64_t res1 = alpha1 * a;
        res1 = right_shift_const(res1, FLOAT_SCALE_FACTOR, pt, role);
        res1 = left_shift_const(res1,scale_factor,pt,role);
        res1 = right_shift(res1, digits,pt,role);

        if(role==SERVER) {
            res1 = res1 + beta1;
        }

        uint64_t res2 = alpha2 * a;
        res2 = right_shift_const(res2, FLOAT_SCALE_FACTOR, pt, role);
        res2 = left_shift_const(res2,scale_factor,pt,role);
        res2 = right_shift(res2, digits,pt,role);
        if(role==SERVER){
            res2=res2 + beta2;
        }

        uint64_t p1 = product(cmp_res, res1, pt, role);
        uint64_t p2 = product(cmp_inv, res2, pt, role);

        uint64_t res = p1 + p2 + digits * (1<<scale_factor);
        if(role == SERVER) {
            res-=already_scaled_factor * (1<<scale_factor);
        }
        return res;
    }

    uint64_t sigmoid(uint64_t a, uint64_t scale_factor, uint64_t already_scaled_factor, ABYParty *pt, e_role role) {

        uint64_t  padded_value = (uint64_t(1) << uint(32)) << (already_scaled_factor);

        vector<uint64_t> alpha(6);
        vector<uint64_t> beta(6);
        alpha[0] = uint64_t (0.07113105 * (1<<FLOAT_SCALE_FACTOR));
        alpha[1] = uint64_t (0.14951886 * (1<<FLOAT_SCALE_FACTOR));
        alpha[2] = uint64_t (0.23810934 * (1<<FLOAT_SCALE_FACTOR));
        alpha[3] = uint64_t (0.14951886 * (1<<FLOAT_SCALE_FACTOR));
        alpha[4] = uint64_t (0.07113105 * (1<<FLOAT_SCALE_FACTOR));
        alpha[5] = uint64_t (0.01930584 * (1<<FLOAT_SCALE_FACTOR));

        beta[0] = uint64_t (0.25621852 * (1<<scale_factor));
        beta[1] = uint64_t (0.41069012 * (1<<scale_factor));
        beta[2] = uint64_t (0.5 * (1<<scale_factor));
        beta[3] = uint64_t (0.58930988 * (1<<scale_factor));
        beta[4] = uint64_t (0.74378148 * (1<<scale_factor));
        beta[5] = uint64_t (0.90177876 * (1<<scale_factor));

//        {-3,-2,-1,1,2,3,5}
        vector<uint64_t> threshold(7);
        threshold[0] = -uint64_t (3 * (1<<already_scaled_factor)) ;
        threshold[1] = -uint64_t (2 * (1<<already_scaled_factor));
        threshold[2] = -uint64_t (1 * (1<<already_scaled_factor));
        threshold[3] = uint64_t (1 * (1<<already_scaled_factor));
        threshold[4] = uint64_t (2 * (1<<already_scaled_factor));
        threshold[5] = uint64_t (3 * (1<<already_scaled_factor));
        threshold[6] = uint64_t (5 * (1<<already_scaled_factor));

        vector<uint64_t>res(8);
        vector<uint64_t>cmp(7);
        res[0] = 0;
        for(int i=0;i<threshold.size();i++) {
            uint64_t tmp = role == SERVER? a+padded_value: a;
            uint64_t cmp_res = share_gt_const(tmp, threshold[i]+padded_value, pt, role);
            cmp[i] = cmp_res;

            if(i==threshold.size()-1) {
                res[i] = role ==SERVER? 1:0;
                break;
            }

            uint64_t val = alpha[i] * a;
            val = right_shift_const(val, FLOAT_SCALE_FACTOR, pt, role);
            val = left_shift_const(val,scale_factor,pt,role);
            val = right_shift_const(val, already_scaled_factor,pt,role);

            if(role==SERVER) {
                val += beta[i];
            }

            res[i+1] = val;
        }

        uint64_t sum=res[0];
        for(int i=1;i<8;i++) {
            sum+= product(cmp[i-1], res[i]-res[i-1], pt, role);
        }

        for(int i=0;i<7;i++) {cout<<cmp[i]<<endl;}

        return sum;
    }

    // the result will be scaled scale_factor   * (1<<scale_factor)
    uint64_t rep_square_root(uint64_t a, uint64_t scale_factor,uint64_t already_scaled_factor,ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        uint64_t digits;
        uint64_t m2N = max2N(a, digits,pt, role);

        uint64_t alpha = -uint64_t(0.8099868542 * (1<<FLOAT_SCALE_FACTOR));
        uint64_t beta = 1.787727479 * (1<<scale_factor);

        uint64_t res = -alpha*a;

        res = right_shift_const(res, FLOAT_SCALE_FACTOR, pt, role);
        res = left_shift_const(res,scale_factor,pt,role);
        res = right_shift(res, digits,pt,role);
        res = -res;
        if(role == SERVER) {
            res += beta;
        }

        uint64_t digits_even_half = right_shift_const(digits,1,pt,role);
        uint64_t even = eq(digits, digits_even_half*2, pt, role);
        uint64_t odd = -even;
        if(role == SERVER) {
            odd+=1;
        }

        uint64_t digits_odd_half = digits_even_half;
        if (role == SERVER) {
            digits_odd_half +=1;
        }

        uint64_t even_res = right_shift(res, digits_even_half, pt, role);
        even_res = product(even,even_res,pt,role);


        uint64_t odd_res = right_shift(res, digits_odd_half, pt, role);
//        cout<<odd_res<<endl;
        odd_res*= uint64_t(std::sqrt(2)* (1<<FLOAT_SCALE_FACTOR));
        odd_res = right_shift_const(odd_res, FLOAT_SCALE_FACTOR, pt, role );
//        cout<<odd_res<<endl;

        odd_res = product(odd,odd_res,pt,role);

        res = even_res+odd_res;
        if(already_scaled_factor!=0) {
            if(already_scaled_factor %2==0) {
                res = left_shift_const(res,already_scaled_factor/2,pt,role);
            } else {
                res = left_shift_const(res,(already_scaled_factor-1)/2,pt,role);

                res *= uint64_t(std::sqrt(2)* (1<<FLOAT_SCALE_FACTOR));
                res = right_shift_const(res, FLOAT_SCALE_FACTOR, pt, role );
            }
        }
        return res;
    }

    uint range_sample(vector<uint64_t>&threshold, vector<uint64_t> &value, uint64_t p, ABYParty *pt, e_role role) {
        assert(value.size() == threshold.size()+1);
        uint64_t sum = value[0];
        for(int i=0;i<threshold.size();i++) {
            uint cmp = gt(p, threshold[i],pt, role);
            sum+=product(cmp, value[i+1]-value[i], pt, role);
        }
        return sum;
    }

    uint64_t min(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto sa = acirc->PutSharedINGate(a,UINT64_LEN);
        auto sb = acirc->PutSharedINGate(b,UINT64_LEN);
        sa = ycirc->PutA2YGate(sa);
        sb = ycirc->PutA2YGate(sb);
        share**val;
        val = (share **)malloc(sizeof(share *)*2);
        val[0] = sa;
        val[1] = sb;
        auto soutput = ycirc->PutMinGate(val,2);
        soutput = acirc->PutY2AGate(soutput, bcirc);
        soutput = acirc->PutSharedOUTGate(soutput);
        pt->ExecCircuit();
        uint64_t out = soutput->get_clear_value<uint64_t>();
        pt->Reset();
        return out;
    }

    vector<uint64_t> argmax_vector(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        uint64_t dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
//        ArithmeticCircuit*
        auto boolcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto yaocirc = (ArithmeticCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        share ** val, ** id, *maxval, *maxindex;
        val = (share **)malloc(sizeof(share *)*dim);
        id = (share **)malloc(sizeof(share *)*dim);
        for(uint64_t i=0; i<a.size(); i++) {
            val[i] = arithcirc->PutSharedINGate(a[i],UINT64_LEN);
            id[i] = arithcirc->PutCONSGate(i, UINT64_LEN);
        }

        for(uint64_t i=0; i<a.size(); i++) {
            val[i] = boolcirc->PutA2BGate(val[i],yaocirc);
            id[i] = boolcirc->PutA2BGate(id[i],yaocirc);
        }
        boolcirc->PutMaxIdxGate(val, id, dim, &maxval, &maxindex);

        vector<uint64_t> indexs(dim);
        vector<share*> sindex;
        for(uint64_t i=0;i<dim;i++) {
            auto tmp_index = boolcirc->PutINGate(i,UINT64_LEN,SERVER);
            auto cmp_res = boolcirc->PutEQGate(maxindex,tmp_index);
            cmp_res = arithcirc->PutB2AGate(cmp_res);
            cmp_res = arithcirc->PutOUTGate(cmp_res,ALL);
            sindex.push_back(cmp_res);
        }
        pt->ExecCircuit();
        for(int i=0;i<dim;i++) {
            uint64_t index = sindex[i]->get_clear_value<uint64_t>();
            indexs[i] = index;
            cout<<index<<endl;
        }
        pt->Reset();

        delete val;
        delete id;
        return indexs;
    }

    uint64_t share_gt_const(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        auto sb = acirc->PutCONSGate(b, UINT64_LEN);
        auto ba = ycirc->PutA2YGate(sa);
        auto bb = ycirc->PutA2YGate(sb);
        auto cmp = ycirc->PutGTGate(ba,bb);
        auto cmp_result = acirc->PutY2AGate(cmp, bcirc);
        cmp_result = acirc->PutSharedOUTGate(cmp_result);
        pt->ExecCircuit();
        auto output = cmp_result->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    uint64_t gt(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto s_out = build_gt_circuit(a, b, pt, role);
        s_out = circ->PutSharedOUTGate(s_out);
        pt->ExecCircuit();
        auto output = s_out->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    uint64_t share_eq_const(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        auto sb = acirc->PutCONSGate(b, UINT64_LEN);
        auto ba = ycirc->PutA2YGate(sa);
        auto bb = ycirc->PutA2YGate(sb);
        auto cmp = ycirc->PutEQGate(ba,bb);
        auto cmp_result = acirc->PutY2AGate(cmp, bcirc);

        cmp_result = acirc->PutSharedOUTGate(cmp_result);
        pt->ExecCircuit();
        auto output = cmp_result->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    uint64_t eq(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto s_out = build_eq_circuit(a, b, pt, role);
        s_out = circ->PutSharedOUTGate(s_out);
        pt->ExecCircuit();
        auto output = s_out->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    vector<uint64_t> product(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty*pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto s_out = build_product_circuit(a, b, pt, role);
        s_out = circ->PutSharedOUTGate(s_out);
        pt->ExecCircuit();
        uint *v;
        uint bitlen, nval;
        s_out->get_clear_value_vec(&v, &bitlen, &nval);
        vector<uint64_t> output (v,v+a.size());
        pt->Reset();
        return output;
    }

    uint64_t product(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto s_out = build_product_circuit(a, b, pt, role);
        s_out = circ->PutSharedOUTGate(s_out);
        pt->ExecCircuit();
        auto output = s_out->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    uint64_t inner_product(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto s_out = build_inner_product_circuit(a, b, pt, role);
        s_out = circ->PutSharedOUTGate(s_out);
        pt->ExecCircuit();
        auto output = s_out->get_clear_value<uint64_t>();
        pt->Reset();
        return output;
    }

    uint64_t argmax2(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        uint64_t dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();

        uint64_t v_max = a[0];
        uint64_t v_index = 0;

        for(uint64_t i=1; i<a.size();i++) {
            uint64_t cmp_result = gt(a[i], v_max, pt, role);
            uint64_t inv_cmp = -cmp_result;
            if(role == SERVER) {
                inv_cmp +=1;
            }

            v_max = product(cmp_result, a[i], pt, role) +
                    product(inv_cmp, v_max, pt, role);
            cout<<v_max<<endl;

            uint64_t index = 0;
            if(i%2==0) {
                index = i/2;
            } else {
                index = role==SERVER? (i+1)/2:(i-1)/2;
            }
            v_index = product(cmp_result, index, pt, role) +
                      product(inv_cmp, v_index, pt, role);
        }
        return v_index;
    }

    share* build_product_circuit(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty*pt, e_role role){
        assert(a.size() == b.size());
        uint64_t dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
        ArithmeticCircuit* circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();

        share * s_a_vec = circ->PutSharedSIMDINGate(dim, a.data(), UINT64_LEN);
        share * s_b_vec = circ->PutSharedSIMDINGate(dim, b.data(), UINT64_LEN);

        // generate multiplication triplet
        MulTripletVec triplet(dim, role);
        auto trip_a = circ->PutSharedSIMDINGate(dim, triplet.getA().data(), UINT64_LEN);
        auto trip_b = circ->PutSharedSIMDINGate(dim, triplet.getB().data(), UINT64_LEN);

        share* su = circ->PutSUBGate(s_a_vec, trip_a);
        share* sv = circ->PutSUBGate(s_b_vec, trip_b);

        su = circ->PutOUTGate(su,ALL);
        sv = circ->PutOUTGate(sv, ALL);

        pt->ExecCircuit();

        uint out_bitlen , out_nvals;
        uint64_t *u, *v;
        su->get_clear_value_vec(&u, &out_bitlen, &out_nvals);
        sv->get_clear_value_vec(&v, &out_bitlen, &out_nvals);
        vector<uint64_t>uv_div2(dim);
        for(uint64_t i=0; i<dim;i++) {
            uv_div2[i] = u[i] * v[i];
            if(uv_div2[i] %2 !=0) {
                uv_div2[i] = role == SERVER? (uv_div2[i]-1)/2: (uv_div2[i]+1)/2;
            } else {
                uv_div2[i]/=2;
            }
        }

        pt->Reset();

        auto va = triplet.getA().data();
        auto ub = triplet.getB().data();
        for(uint64_t i = 0; i< dim; i++) {
            va[i] *= v[i];
            ub[i] *= u[i];
        }

        auto s_va = circ->PutSharedSIMDINGate(dim, va, UINT64_LEN);
        auto s_ub = circ->PutSharedSIMDINGate(dim, ub, UINT64_LEN);
        auto s_c = circ->PutSharedSIMDINGate(dim, triplet.getC().data(), UINT64_LEN);

        share* c_uv_div2 = circ->PutSharedSIMDINGate(dim,uv_div2.data(), UINT64_LEN);

        s_va = circ->PutADDGate(s_va, s_ub);
        s_va = circ->PutADDGate(s_va, s_c);
        s_va = circ->PutADDGate(s_va, c_uv_div2);


        delete u;
        delete v;
        return s_va;
    }

    share* build_product_circuit(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        std::vector<Sharing*>& sharings = pt->GetSharings();
        auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();

        share * sa = circ->PutSharedINGate(a, UINT64_LEN);
        share * sb = circ->PutSharedINGate(b, UINT64_LEN);

        MulTriplet triplet(role);
        auto trip_a = circ->PutSharedINGate(triplet.getA(), UINT64_LEN);
        auto trip_b = circ->PutSharedINGate(triplet.getB(), UINT64_LEN);

        share* su = circ->PutSUBGate(sa, trip_a);
        share* sv = circ->PutSUBGate(sb, trip_b);

        su = circ->PutOUTGate(su,ALL);
        sv = circ->PutOUTGate(sv, ALL);

        pt->ExecCircuit();

        uint64_t u, v;
        u = su->get_clear_value<uint64_t>();
        v = sv->get_clear_value<uint64_t>();
        uint64_t uv_div2 = u*v;
        if(uv_div2 %2 !=0) {
            uv_div2 = role == SERVER? (uv_div2-1)/2: (uv_div2+1)/2;
        } else {
            uv_div2/=2;
        }

        pt->Reset();

        auto va = triplet.getA();
        auto ub = triplet.getB();
        va*=v;
        ub*=u;


        auto s_va = circ->PutSharedINGate(va, UINT64_LEN);
        auto s_ub = circ->PutSharedINGate(ub, UINT64_LEN);
        auto s_c = circ->PutSharedINGate(triplet.getC(), UINT64_LEN);

        share* c_uv_div2 = circ->PutSharedINGate(uv_div2, UINT64_LEN);

        s_va = circ->PutADDGate(s_va, s_ub);
        s_va = circ->PutADDGate(s_va, s_c);
        s_va = circ->PutADDGate(s_va, c_uv_div2);

        return s_va;
    }

    share* build_inner_product_circuit(vector<uint64_t> &a, vector<uint64_t> &b, ABYParty *pt, e_role role) {
        assert(a.size() == b.size());
        uint64_t dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
        ArithmeticCircuit* circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();

        share * s_a_vec = circ->PutSharedSIMDINGate(dim, a.data(), UINT64_LEN);
        share * s_b_vec = circ->PutSharedSIMDINGate(dim, b.data(), UINT64_LEN);

        // generate multiplication triplet
        MulTripletVec triplet(dim, role);
        auto trip_a = circ->PutSharedSIMDINGate(dim, triplet.getA().data(), UINT64_LEN);
        auto trip_b = circ->PutSharedSIMDINGate(dim, triplet.getB().data(), UINT64_LEN);

        share* su = circ->PutSUBGate(s_a_vec, trip_a);
        share* sv = circ->PutSUBGate(s_b_vec, trip_b);

        su = circ->PutOUTGate(su,ALL);
        sv = circ->PutOUTGate(sv, ALL);

        pt->ExecCircuit();

        uint out_bitlen , out_nvals;
        uint64_t *u, *v;
        su->get_clear_value_vec(&u, &out_bitlen, &out_nvals);
        sv->get_clear_value_vec(&v, &out_bitlen, &out_nvals);
        uint64_t uv_div2 = 0;
        for(uint64_t i=0; i<dim;i++) {
            uv_div2 += u[i] * v[i];
        }
        if(uv_div2 %2 !=0) {
            uv_div2 = role == SERVER? (uv_div2-1)/2: (uv_div2+1)/2;
        } else {
            uv_div2/=2;
        }

        pt->Reset();

        auto va = triplet.getA().data();
        auto ub = triplet.getB().data();
        for(uint64_t i = 0; i< dim; i++) {
            va[i] *= v[i];
            ub[i] *= u[i];
        }

        auto s_va = circ->PutSharedSIMDINGate(dim, va, UINT64_LEN);
        auto s_ub = circ->PutSharedSIMDINGate(dim, ub, UINT64_LEN);
        auto s_c = circ->PutSharedSIMDINGate(dim, triplet.getC().data(), UINT64_LEN);

        share* c_uv_div2 = circ->PutSharedINGate(uv_div2, UINT64_LEN);

        s_va = circ->PutADDGate(s_va, s_ub);
        s_va = circ->PutADDGate(s_va, s_c);
        s_va = circ->PutSplitterGate(s_va);
        for (uint64_t i = 1; i < dim; i++) {
            s_va->set_wire_id(0, circ->PutADDGate(s_va->get_wire_id(0), s_va->get_wire_id(i)));
        }
        s_va->set_bitlength(1);
        s_va = circ->PutADDGate(s_va, c_uv_div2);


        delete u;
        delete v;
        return s_va;
    }

    share* build_logarithm_circuit(vector<uint64_t>&a, ABYParty *party) {
        share *s_a_vec;
        return s_a_vec;
    }

    share *build_argmax_circuit(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        uint64_t dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
//        ArithmeticCircuit*
        auto boolcirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        share ** val, ** id, *maxval, *maxindex;
        val = (share **)malloc(sizeof(share *)*dim);
        id = (share **)malloc(sizeof(share *)*dim);
        for(uint64_t i=0; i<a.size(); i++) {
            val[i] = arithcirc->PutSharedINGate(a[i],UINT64_LEN);
            id[i] = arithcirc->PutCONSGate(i, UINT64_LEN);
        }

        for(uint64_t i=0; i<a.size(); i++) {
            val[i] = boolcirc->PutA2YGate(val[i]);
            id[i] = boolcirc->PutA2YGate(id[i]);
        }
        boolcirc->PutMaxIdxGate(val, id, dim, &maxval, &maxindex);

        delete val;
        delete id;
        return maxindex;
    }

    share *build_argmax_circuit2(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        uint64_t dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
//        ArithmeticCircuit*
        auto boolcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto yaocirc = (ArithmeticCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        share ** val, ** id, *maxval, *maxindex;
        val = (share **)malloc(sizeof(share *)*dim);
        id = (share **)malloc(sizeof(share *)*dim);
        for(uint64_t i=0; i<a.size(); i++) {
            val[i] = arithcirc->PutSharedINGate(a[i],UINT64_LEN);
            id[i] = arithcirc->PutCONSGate(i, UINT64_LEN);
        }

        for(uint64_t i=0; i<a.size(); i++) {
            val[i] = boolcirc->PutA2BGate(val[i],yaocirc);
            id[i] = boolcirc->PutA2BGate(id[i],yaocirc);
        }

        delete val;
        delete id;
        return maxindex;
    }

    share *build_argmax_vector_circuit(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        uint dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
//        ArithmeticCircuit*
        auto boolcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto yaocirc = (ArithmeticCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        share ** val, ** id, *maxval, *maxindex;
        val = (share **)malloc(sizeof(share *)*dim);
        id = (share **)malloc(sizeof(share *)*dim);
        for(uint i=0; i<a.size(); i++) {
            val[i] = arithcirc->PutSharedINGate(a[i],UINT64_LEN);
            id[i] = arithcirc->PutCONSGate(i, UINT64_LEN);
        }

        for(uint i=0; i<a.size(); i++) {
            val[i] = boolcirc->PutA2BGate(val[i],yaocirc);
            id[i] = boolcirc->PutA2BGate(id[i],yaocirc);
        }
        boolcirc->PutMaxIdxGate(val, id, dim, &maxval, &maxindex);
        vector<uint> indexs(dim);
        vector<share*> sindex;
        for(uint i=0;i<dim;i++) {
            auto tmp_index = boolcirc->PutCONSGate(i,UINT64_LEN);
            auto cmp_res = boolcirc->PutEQGate(maxindex,tmp_index);
            cmp_res = arithcirc->PutB2AGate(cmp_res);
            cmp_res = arithcirc->PutOUTGate(cmp_res,ALL);
            sindex.push_back(cmp_res);
        }
        pt->ExecCircuit();
        for(uint i=0;i<dim;i++) {
            uint index = sindex[i]->get_clear_value<uint>();
            indexs[i] = index;
            cout<<index<<endl;
        }
        pt->Reset();

        delete val;
        delete id;
        return maxindex;
    }

    share*build_gt_circuit(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        std::vector<Sharing*>& sharings = pt->GetSharings();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto sa = arithcirc->PutSharedINGate(a, UINT64_LEN);
        auto sb = arithcirc->PutSharedINGate(b, UINT64_LEN);
        auto ba = ycirc->PutA2YGate(sa);
        auto bb = ycirc->PutA2YGate(sb);
        auto cmp = ycirc->PutGTGate(ba,bb);
        auto cmp_result = arithcirc->PutY2AGate(cmp, bcirc);
        return cmp_result;
    }

    share*build_eq_circuit(uint64_t a, uint64_t b, ABYParty *pt, e_role role) {
        std::vector<Sharing*>& sharings = pt->GetSharings();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        auto sa = arithcirc->PutSharedINGate(a, UINT64_LEN);
        auto sb = arithcirc->PutSharedINGate(b, UINT64_LEN);
        auto ba = ycirc->PutA2YGate(sa);
        auto bb = ycirc->PutA2YGate(sb);
        auto cmp = ycirc->PutEQGate(ba,bb);
        auto cmp_result = arithcirc->PutY2AGate(cmp, bcirc);
        return cmp_result;
    }

    share *build_max_circuit(vector<uint64_t>&a, ABYParty *pt, e_role role) {
        uint dim = a.size();
        std::vector<Sharing*>& sharings = pt->GetSharings();
//        ArithmeticCircuit*
        auto boolcirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto arithcirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        share ** val,  *maxval;
        val = (share **)malloc(sizeof(share *)*dim);
        for(uint i=0; i<a.size(); i++) {
            val[i] = arithcirc->PutSharedINGate(a[i],UINT64_LEN);
        }

        for(uint i=0; i<a.size(); i++) {
            val[i] = boolcirc->PutA2YGate(val[i]);
        }
        auto out = boolcirc->PutMaxGate(val, dim);


        delete val;
        return out;
    }

    share* build_max2N_circuit(share* sa, share ** digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
        sa = ycirc->PutA2YGate(sa);
        uint sum = 0;
        digits = 0;
        vector<share*>inds(UINT64_LEN);
        for(int i=UINT64_LEN-1;i>=0;i--) {
            uint v = uint(1)<<uint(i);
            auto sv = ycirc->PutINGate(v,UINT64_LEN,SERVER);
            auto cmp_res = ycirc->PutGTGate(sa,sv);
            cmp_res = acirc->PutY2AGate(cmp_res, bcirc);
            inds[i]= cmp_res;
        }
        vector<share*> svs;
        for(int i=UINT64_LEN-1; i>=0; i--) {
            uint v = uint(1)<<uint(i);
            auto sv = acirc->PutCONSGate(v, UINT64_LEN);
            sv = acirc->PutMULCONSTGate(sv, inds[i]);
            svs.push_back(sv);
        }

        for(int i=1;i<UINT64_LEN;i++) {
            svs[0] = acirc->PutADDGate(svs[0],svs[i]);
            inds[0] = acirc->PutADDGate(inds[0], inds[1]);
        }
        *digits = inds[0];
        return svs[0];
    }

    share* build_max2N_circuit(uint64_t a, share ** digits, ABYParty *pt, e_role role) {
        auto sharings = pt->GetSharings();
        auto acirc = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
        auto ycirc = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
        auto bcirc = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();

        auto sa = acirc->PutSharedINGate(a, UINT64_LEN);
        return build_max2N_circuit(sa, digits, pt, role);
    }
}