////
//// Created by anxin on 8/9/20.
////
//
//#include "test.h"
//
//
//
//void test_inner_product(ABYParty *pt, e_role role) {
//    int len = 100;
//    srand(1);
//    vector<uint64_t> a(len);
//    vector<uint64_t> b(len);
//    uint64_t v_sum = 0;
//    for(int i=0;i<len;i++) {
//        a[i] = rand()%100;
//        b[i] = rand()%100;
//        v_sum += a[i] * b[i];
//    }
//
//    MPC::VecShare a_share (a,role);
//    MPC::VecShare b_share (b,role);
//
//    std::vector<Sharing*>& sharings = pt->GetSharings();
//    ArithmeticCircuit* circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
//
//
//    auto s_out = MPC::build_inner_product_circuit(a_share.val(),b_share.val(),pt, role);
//    s_out = circ->PutOUTGate(s_out, ALL);
//
//    pt->ExecCircuit();
//
//    auto output = s_out->get_clear_value<uint64_t>();
//    pt->Reset();
////    assert(output == v_sum);
////
//    std::cout << "\nCircuit Result: " << output;
//    std::cout << "\nVerification Result: " << v_sum << std::endl;
//}
//
//void test_argmax(ABYParty *pt, e_role role) {
//    int len = 20;
//    srand(1);
//    vector<uint64_t> a(len);
//
//    for(int i=0;i<len;i++) {
//        a[i] = i+1;
//    }
//    shuffle(a.begin(), a.end(),std::default_random_engine(1));
//    int v_max = 0;
//    int v_index = 0;
//    for(int i=0;i<len;i++) {
//        if(v_max<a[i]) {
//            v_max = a[i];
//            v_index = i;
//        }
//    }
//
//
//    MPC::VecShare a_share (a,role);
//    std::vector<Sharing*>& sharings = pt->GetSharings();
//    BooleanCircuit* circ = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
////    ArithmeticCircuit *arc = (ArithmeticCircuit *) sharings[S_ARITH]->GetCircuitBuildRoutine();
////
////    len = 50;
////    vector<uint64_t> t1(len);
////    vector<uint64_t> t2(len);
////
////    auto x = circ->PutSIMDINGate(len,t1.data(),UINT64_LEN, SERVER);
////    auto y =  circ->PutSIMDINGate(len,t2.data(),UINT64_LEN,CLIENT);
////    auto output = circ->PutGTGate(x,y);
////    output = circ->PutOUTGate(output, ALL);
//
////    pt->ExecCircuit();
//
//    auto s_out = MPC::build_argmax_circuit(a_share.val(), pt, role);
//    auto s_out2 = MPC::build_argmax_circuit(a_share.val(), pt, role);
//
//    s_out = circ->PutOUTGate(s_out, ALL);
//    s_out2 = circ->PutOUTGate(s_out2, ALL);
//
//    pt->ExecCircuit();
//
//    auto output = s_out->get_clear_value<uint64_t>();
//
//    pt->Reset();
//
//    std::cout << "\nCircuit Result: " << output;
//    std::cout << "\nVerification Result: " << v_index << std::endl;
//    cout<<a[v_index]<<endl;
//}
//
//void test_argmax2(ABYParty *pt, e_role role) {
//    int len = 10;
//    srand(1);
//    vector<uint64_t> a(len);
//
//    for(int i=0;i<len;i++) {
//        a[i] = i+1;
//    }
//    shuffle(a.begin(), a.end(),std::default_random_engine(1));
//    int v_max = 0;
//    int v_index = 0;
//    for(int i=0;i<len;i++) {
//        if(v_max<a[i]) {
//            v_max = a[i];
//            v_index = i;
//        }
//    }
//
//
//    MPC::VecShare a_share (a,role);
//    std::vector<Sharing*>& sharings = pt->GetSharings();
//    BooleanCircuit* circ = (BooleanCircuit*) sharings[S_BOOL]->GetCircuitBuildRoutine();
//
//    auto s_out = MPC::build_argmax_circuit2(a_share.val(), pt, role);
//
//    s_out = circ->PutOUTGate(s_out, ALL);
//
//    pt->ExecCircuit();
//
//    auto output = s_out->get_clear_value<uint64_t>();
//
//    pt->Reset();
//
//    std::cout << "\nCircuit Result: " << output;
//    std::cout << "\nVerification Result: " << v_index << std::endl;
//    cout<<a[v_index]<<endl;
//}
//
//void test_argmax3(ABYParty *pt, e_role role) {
//    int len = 10;
//    srand(1);
//    vector<uint64_t> a(len);
//
//    for(int i=0;i<len;i++) {
//        a[i] = i+1;
//    }
//    shuffle(a.begin(), a.end(),std::default_random_engine(1));
//    int v_max = 0;
//    int v_index = 0;
//    for(int i=0;i<len;i++) {
//        if(v_max<a[i]) {
//            v_max = a[i];
//            v_index = i;
//        }
//    }
//
//
//    MPC::VecShare a_share (a,role);
//    std::vector<Sharing*>& sharings = pt->GetSharings();
//    auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
//
//    auto s_out = MPC::argmax(a_share.val(), pt, role);
//
//
//    std::cout << "\nCircuit Result: " << s_out;
//    std::cout << "\nVerification Result: " << v_index << std::endl;
//    cout<<a[v_index]<<endl;
//}
//
//void test_argmax4(ABYParty *pt, e_role role) {
//    int len = 10;
//    srand(1);
//    vector<uint64_t> a(len);
//
//    for(int i=0;i<len;i++) {
//        a[i] = i+1;
//    }
//    shuffle(a.begin(), a.end(),std::default_random_engine(1));
//    int v_max = 0;
//    int v_index = 0;
//    for(int i=0;i<len;i++) {
//        if(v_max<a[i]) {
//            v_max = a[i];
//            v_index = i;
//        }
//    }
//
//
//    MPC::VecShare a_share (a,role);
//    std::vector<Sharing*>& sharings = pt->GetSharings();
//    auto circ = (ArithmeticCircuit*) sharings[S_ARITH]->GetCircuitBuildRoutine();
//
//    auto s_out = MPC::argmax_vector(a_share.val(), pt, role);
//
//    for(uint i=0;i<len;i++) {
//        if(s_out[i] == 1) {
//            std::cout << "\nCircuit Result: " << i<<endl;
//
//        }
//    }
//    std::cout << "\nVerification Result: " << v_index << std::endl;
//    cout<<a[v_index]<<endl;
//}
//
//void test_max(ABYParty *pt, e_role role) {
//    int len = 2;
//    srand(1);
//    vector<uint64_t> a(len);
//
//    for(int i=0;i<len;i++) {
//        a[i] = i+1;
//    }
//    shuffle(a.begin(), a.end(),std::default_random_engine(1));
//    int v_max = 0;
//    int v_index = 0;
//    for(int i=0;i<len;i++) {
//        if(v_max<a[i]) {
//            v_max = a[i];
//            v_index = i;
//        }
//    }
//
//
//    MPC::VecShare a_share (a,role);
//    std::vector<Sharing*>& sharings = pt->GetSharings();
//    BooleanCircuit* circ = (BooleanCircuit*) sharings[S_YAO]->GetCircuitBuildRoutine();
//
//    auto s_out = MPC::build_max_circuit(a_share.val(), pt, role);
//
//    s_out = circ->PutOUTGate(s_out, ALL);
//
//    pt->ExecCircuit();
//
//    auto output = s_out->get_clear_value<uint64_t>();
//
//    pt->Reset();
//
//    std::cout << "\nCircuit Result: " << output;
//    std::cout << "\nVerification Result: " << v_max << std::endl;
//    cout<<a[v_index]<<endl;
//}
//
//void test_min(ABYParty *pt, e_role role) {
//    uint64_t a =5;
//    uint64_t b= 10;
//    auto output = MPC::min(a,b,pt,role);
//    cout<<output<<endl;
//}
//
//void test_max2N(ABYParty *pt, e_role role) {
//    uint64_t a=1024;
//    cout<< MPC::max2N(a,pt,role) <<endl;
//}
//
//void test_right_shift(ABYParty *pt, e_role role) {
//    srand(1);
//    uint64_t a = 2048;
//    uint64_t tmp = rand();
//    if (role == SERVER) {
//        a = tmp;
//    } else {
//        a = a-tmp;
//    }
//    auto res = MPC::right_shift(a,3,pt,role);
////    cout<<a<<endl;
//    cout<<res<<endl;
//}
//
//void test_log(ABYParty *pt, e_role role) {
//    clock_t start, end;
//    start = clock();
//    for(int i=0; i<1;i++) {
//        uint64_t a = 21<<FLOAT_SCALE_FACTOR;
//        srand(2);
//        uint64_t x = rand();
//        a = role==SERVER? a-x:x;
//
//        uint64_t res = MPC::log(a,FLOAT_SCALE_FACTOR,FLOAT_SCALE_FACTOR,pt,role);
//        cout<<res<<endl;
//    }
//    end = clock();
//    cout<<"Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
////    cout<<res<<endl;
//}
//
//void test_rep_square_root(ABYParty *pt, e_role role) {
//    clock_t start, end;
//    start = clock();
//    for(int i=0; i<1;i++) {
//        uint64_t a = (0)<<20;
//        srand(2);
//        uint64_t x = rand();
//        a = role==SERVER? a-x:x;
//        auto res = MPC::rep_square_root(a, 20,20, pt, role);
//        cout<<res<<endl;
//    }
//    end = clock();
//    cout<<"Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
////    cout<<res<<endl;
//}
//
//void test_sigmoid(ABYParty *pt, e_role role) {
//    clock_t start, end;
//    start = clock();
//    for(int i=0; i<1;i++) {
//        uint64_t a = (3)<<16;
//        srand(2);
//        uint64_t x = rand();
//        a = role==SERVER? a-x:x;
//        auto res = MPC::sigmoid(a, 16,16, pt, role);
//        cout<<res<<endl;
//    }
//    end = clock();
//    cout<<"Run time: "<<(double)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
////    cout<<res<<endl;
//}
//
//void test_eq(ABYParty *pt, e_role role) {
//    auto res = MPC::eq(3,2,pt,role);
//    cout<<res<<endl;
//}
//
//void test_right_shift2(ABYParty *pt, e_role role) {
//    uint64_t a = role==SERVER? 4:5;
//    auto res = MPC::right_shift(a,1,pt,role);
//    cout<<res<<endl;
//}
