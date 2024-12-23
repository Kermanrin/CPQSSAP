#pragma once
#ifndef NTRU_H
#define NTRU_H

#include "Polynomial.h"
#include <iostream>
#include <unordered_map>

class NTRU {
public:
    
    NTRU(int sl) : N(security_level_table[sl % 4][0]), p(security_level_table[sl % 4][1]), q(security_level_table[sl % 4][2]), cntf_1(security_level_table[sl % 4][3])
        , cntf_neg1(security_level_table[sl % 4][4]), cntg_1(security_level_table[sl % 4][5]), cntg_neg1(security_level_table[sl % 4][6]), f(KeyGen_f()), g(KeyGen_g()) {};

    Polynomial KeyGen_f();//生成私钥f
    Polynomial KeyGen_g();//生成私钥g
    Polynomial KeyGen_h();//生成公钥h,需要调用Polynomial类的求多项式逆函数，还没实现
    Polynomial KeyGen_f_mod(int mod);//计算多项式模mod的逆，需要调用Polynomial类的求多项式逆函数，还没实现

    Polynomial get_sk_f() { return f; };
    Polynomial get_sk_g() { return g; };
    Polynomial get_pk_h() { return h; };


private:
    friend class Polynomial;
    int security_level_table[4][7] = { {3, 3, 100, 2, 0, 1, 1},{512, 3, 10000001, 10000, 10000, 10000, 9999},{768, 3, 10000001, 10000, 9999}, {1024, 3, 10000001, 10000, 9999} };
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
    Polynomial h;//公钥
    Polynomial f_q;//f mod q的逆
    Polynomial f_p;// f mod p 的逆
    
};


#endif

