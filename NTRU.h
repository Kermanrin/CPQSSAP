#pragma once
#ifndef NTRU_H
#define NTRU_H

#include "Polynomial.h"
#include <iostream>
#include <unordered_map>

class NTRU {
public:
    
    NTRU(std::vector<int> vec, int sl) : N(security_level_table[sl % 4][0]), p(security_level_table[sl % 4][1]), q(security_level_table[sl % 4][2]), cntf_1(security_level_table[sl % 4][3])
        , cntf_neg1(security_level_table[sl % 4][4]), cntg_1(security_level_table[sl % 4][5]), cntg_neg1(security_level_table[sl % 4][6]), h(KeyGen_h()), f(KeyGen_f()), g(KeyGen_g())
        , f_q(KeyGen_f_q()), f_p(KeyGen_f_p()) {
        int len = vec.size();
        for (int i = 0; i < len; i++) {
            vec[i] = vec[i] % p - 1;
        }
        m = Polynomial(N, q, p, vec);
    }
    Polynomial KeyGen_f();//生成私钥f
    Polynomial KeyGen_g();//生成私钥g
    Polynomial KeyGen_h();//生成公钥h,需要调用Polynomial类的求多项式逆函数
    Polynomial KeyGen_f_q();
    Polynomial KeyGen_f_p();
    Polynomial KeyGen_f_mod(int mod);//计算多项式模mod的逆，需要调用Polynomial类的求多项式逆函数，还没实现

    Polynomial get_sk_f() { 
        printf("The sk_f is : ");
        return f; 
    };
    Polynomial get_sk_g() {
        printf("The sk_g is : ");
        return g; 
    };
    Polynomial get_pk_h() {
        printf("The pk_h is : ");
        return h; 
    };
    Polynomial get_sk_f_q() {
        printf("The sk_f_q is : ");
        return f_q; 
    };
    Polynomial get_sk_f_p() {
        printf("The sk_f_p is : ");
        return f_p; 
    };

    Polynomial get_message() {
        printf("The plaintext is : ");
        return m;
    };


private:
    friend class Polynomial;
    int security_level_table[4][7] = { {3, 3, 100, 2, 0, 1, 1},{251, 3, 257, 72, 71, 72, 72},{768, 3, 10000001, 10000, 9999}, {1024, 3, 10000001, 10000, 9999} };
    int security_level;
    unsigned int const N;
    unsigned int const p;//通常为3,不特属于某个成员的变量，必须在类外初始化
    unsigned int const q;//不特属于某个成员的statec变量
    unsigned int const cntf_1;
    unsigned int const cntf_neg1;
    unsigned int const cntg_1;
    unsigned int const cntg_neg1;
    Polynomial f;//私钥
    Polynomial g;//私钥
    Polynomial f_q;//f mod q的逆
    Polynomial f_p;// f mod p 的逆
    Polynomial h;//公钥
    Polynomial m;
    
};


#endif

